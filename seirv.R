
seirv_model <- function(
    N, 
    beta, 
    delta, 
    gamma, 
    nu = function(t, mat) {
      (t >= 50) * (1 - exp(-0.02)) * (1 - mat[, 5]) *
        (mat[, 1] == 1 | mat[, 2] == 1)
    }, 
    steps = 365
  ) {
  
  # initialize compartment assignments
  S <- c(0, rep(1, N - 1))
  E <- c(rep(1, 10), rep(0, N - 10))
  I <- rep(0, N)
  R <- rep(0, N)
  V <- rep(0, N)
  
  # initialize matrix
  mat <- cbind(S, E, I, R, V)
  
  # initialize sim list
  sims <- list()
  
  for (t in 1:steps) {
    # calculate compartment populations
    pop <- colSums(mat)
      
    # transition probabilities
    p_SE <- 1 - exp(-(beta * pop[2] / N) * mat[, 1]) # S to E
    p_EI <- 1 - exp(-delta * mat[, 2]) # E to I 
    p_IR <- 1 - exp(-gamma * mat[, 3]) # I to R
    
    p_V <- nu(t, mat) 
    
    # simulate transitions
    sims[[t]] <- make_transitions(mat, p_SE, p_EI, p_IR, p_V)
    
    mat <- sims[[t]]
  }
  
  # bind all sims together and create time and id cols
  sims <- rbindlist(lapply(sims, as.data.table), idcol = "step")
  sims$id <- rep(1:N, steps)
    
  return(sims)
}

make_transitions <- function(mat, p_SE, p_EI, p_IR, p_V) {
  # copy matrix
  mat_new <- mat

  # draw transitions
  SE <- rbinom(nrow(mat), size = 1, p_SE)
  EI <- rbinom(nrow(mat), size = 1, p_EI)
  IR <- rbinom(nrow(mat), size = 1, p_IR)
  
  # draw vaccinations
  V <- rbinom(nrow(mat), size = 1, p_V)
  
  # assign new states
  mat_new[SE == 1, 1] <- 0
  mat_new[SE == 1, 2] <- 1
  
  mat_new[EI == 1, 2] <- 0
  mat_new[EI == 1, 3] <- 1
  
  mat_new[IR == 1, 3] <- 0
  mat_new[IR == 1, 4] <- 1
  
  mat_new[V == 1, 5] <- 1
  
  return(mat_new)
}

observational_study <- function(sims, study_start = 50, study_end = 100) {
  
  # filter post campaign data
  post <- filter(sims, step >= study_start)
  
  # figure out who has already been infected
  immune_ids <- post$id[which(post$step == study_start & (post$R == 1 | post$I == 1))]
  
  # figure out who is ever vaccinated
  vacc_ids <- post$id[which(post$step == study_end & post$V == 1)]
  
  # filter out those who are immune
  post <- filter(post, !id %in% immune_ids)
  
  # make outcome
  post$case <- post$R + post$I
  
  # calculate event times
  post <-
    post |>
    group_by(id) |>
    mutate(
      case_time = step[cumsum(case) == 1 | (cumsum(case) == 0 & step == study_end)],
      case_time = case_time + runif(n()),
      vacc_time = step[cumsum(V) == 1 | (cumsum(V) == 0 & step == study_end)],
      vaccinated = max(V),
      event = max(case),
      event = replace(event, vaccinated & step < vacc_time, 0),
      cox_keep = ifelse(
        step == study_start | (vaccinated == 1 & step == vacc_time), 
        1, 
        0
      ),
      stop = ifelse(
        vaccinated == 1 & step < vacc_time,
        vacc_time,
        case_time
      ),
      start = 0,
      diff = stop - step
    ) |>
    ungroup()

  # naive logistic regression
  logit <- glm(
    formula = case ~ V,
    family = poisson(link = "log"),
    data = filter(post, step == study_end)
  )
  
  # time-dependent cox regression
  cox <- coxph(
    formula = Surv(step, stop, event) ~ V,
    data = filter(post, cox_keep == 1),
    id = id
  )
  
  # time-dependent cox reset
  cox_reset <- coxph(
    formula = Surv(start, diff, event) ~ V,
    data = filter(post, cox_keep == 1),
    id = id
  )
  
  # target trial
  tt <- coxph(
    formula = Surv(diff, event) ~ V + strata(step),
    data = filter(post, step <= case_time & step <= vacc_time),
    id = id
  )
  
  #return(list(logit, cox, filter(post, step <= case_time & step <= vacc_time) |> arrange(id, step)))
  #return(list(logit, cox, filter(post, cox_keep == 1) |> arrange(id, step)))
  #return(list(filter(post, cox_keep == 1) |> arrange(id, step), filter(post, step <= case_time & step <= vacc_time) |> arrange(id, step)))
  return(list(logit, cox, cox_reset, tt))
}
