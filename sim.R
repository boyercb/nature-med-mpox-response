# load packages
library(data.table)
library(tidyverse)
library(patchwork)
library(survival)
library(splines)


# load user-defined functions
source("seirv.R")


# run simulation ----------------------------------------------------------

# number of monte carlo simulations
SIMS <- 100

set.seed(23415125)

# grid of simulation parameters
sim_grid <- expand_grid(
  nu = list(
    "all together" = function(t, mat) {
      (t == 50) * 0.30 *
        (mat[, 1] == 1 | mat[, 2] == 1)
    },
    "constant rate" = function(t, mat) {
      (t >= 50) * (1 - exp(-0.015)) * (1 - mat[, 5]) *
        (mat[, 1] == 1 | mat[, 2] == 1)
    },
    "accelerating" = function(t, mat) {
      (t >= 50) * (1 - exp(-8e-04 * t + 0.04)) * (1 - mat[, 5]) *
        (mat[, 1] == 1 | mat[, 2] == 1)
    }
  ),
  sim = 1:SIMS
)

# run simulation
with_progress({
  p <- progressor(steps = nrow(sim_grid))
  
  res <- pmap(
    sim_grid,
    function(sim, nu, p) {
      s <- seirv_model(
        N = 10000, 
        beta = 0.35,
        delta = 0.20,
        gamma = 0.20,
        nu = nu,
        steps = 100
      )
      p()
      observational_study(s)
    },
    p = p
  )
})


# make plots --------------------------------------------------------------

examples <- pmap_dfr(
  filter(sim_grid, sim == 1),
  function(sim, nu, p) {
    set.seed(1235)
    s <- seirv_model(
      N = 10000, 
      beta = 0.35,
      delta = 0.20,
      gamma = 0.20,
      nu = nu,
      steps = 100
    )
    
    s[, .(
      S = sum(S), 
      E = sum(E), 
      I = sum(I), 
      R = sum(R), 
      V = sum(V), 
      IV = sum(V * I)
    ), by = step]
    
  },
  .id = "scenario",
  p = p
)

examples <-
  examples |>
  mutate(
    scenario = factor(
      scenario,
      levels = c(
        "all together",
        "constant rate",
        "accelerating"
      )
    )
  )

p1 <- 
  ggplot(
    pivot_longer(examples, c(I, IV)),
    aes(x = step, y = value, group = name, color = name, fill = name)
  ) + 
  facet_grid(~scenario) +
  geom_col(linewidth = 0.75) + 
  geom_vline(xintercept = 50, linetype = "dashed") +
  scale_fill_brewer(
    palette = "Set2",
    name = "",
    labels = c("unvaccinated", "vaccinated")
  ) +
  scale_color_brewer(
    palette = "Set2",
    name = "",
    labels = c("unvaccinated", "vaccinated")
  ) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(
    title = "A",
    x = "time",
    y = "infections\n"
  )

p2 <- 
  ggplot(
    pivot_longer(examples, c(V)),
    aes(x = step, y = value, group = name, color = name, fill = name)) + 
  facet_grid(~scenario) +
  geom_col(linewidth = 0.75, fill = 'grey80', color = 'grey80') + 
  geom_vline(xintercept = 50, linetype = "dashed") +
  theme_minimal() +
  labs(
    title = "B",
    x = "time",
    y = "vaccinated\n"
  )

df <- tibble(
  sim = rep(1:SIMS, 3),
  scenario = rep(
    c("all together",
      "constant rate",
      "accelerating"),
    each = SIMS
  ),
  cox = sapply(res, function(x) coef(x[[2]])),
  cox_reset = sapply(res, function(x) coef(x[[3]])),
  tt = sapply(res, function(x) coef(x[[4]])),
  logit = sapply(res, function(x) coef(x[[1]])[2])
)

df_long <- pivot_longer(df, -c(sim, scenario)) |>
  mutate(
    name = factor(
      name,
      levels = c("logit", "cox_reset", "cox", "tt"),
      labels = c(
        "naive",
        "cox\n(reset T0)",
        "cox",
        "target\ntrial"
      )
    ),
    scenario = factor(
      scenario,
      levels = c(
        "all together",
        "constant rate",
        "accelerating"
      )
    )
  )


p3 <- ggplot(df_long, aes(x = name, y = value)) +
  facet_grid(~scenario) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(fill = 'lightblue') +
  theme_minimal() + 
  labs(
    title = "C",
    x = NULL,
    y = "estimate"
  )



ggsave(
  filename = "figure1.png",
  plot = (p1 / p2 / p3),
  width = 9,
  height = 6.5,
  units = "in",
  dpi = 400
)

