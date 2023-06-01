# load packages
library(ggplot2)
library(tidyverse)
library(patchwork)
library(flextable)


km <- read_csv("KM.csv")

# number of study week observations in KM data
len <- length(km$u_risk) 

# infer number vaxed each week
new_vax_calwk <- diff(rev(c(km$v_risk, 0))) 

# calculate cumulative number vaxed 
cum_vax_calwk <- cumsum(new_vax_calwk) 

# vax_at and cum_vax are indexed by calendar time, not by time since vaccination as the other vaccinated variables.

# total study size calculated from the number ever vaccinated and the number always unvax
N <- km$v_risk[1] + km$u_risk[1] 

# incident infections by calendar week in unvax
U_incinf_calwk <- c(diff(c(0, km$u_event)))

# incidence rate in unvaxed
U_incrate_calwk <-
  U_incinf_calwk / (km$u_risk + km$v_risk[1] - cum_vax_calwk) 

# incident cases each week with placebo vaccine
V_inc_null_calwk <-
  U_incrate_calwk * cum_vax_calwk 

# total incident cases by each week with placebo vax
V_cum_null_calwk <-
  cumsum(V_inc_null_calwk) 

RR_crude_recalc <- km$v_event[len] / V_cum_null_calwk[len]

VE_crude_recalc <- 1 - RR_crude_recalc
VE_crude_recalc

# standardize to size of unvaxxed
V_inc_null_calwk_std <- 
  V_inc_null_calwk * km$u_risk[1] / km$v_risk[1]

V_cum_null_calwk_std <-
  cumsum(V_inc_null_calwk_std)

# incident infections by calendar week in vax
V_incinf_calwk <- c(diff(c(0, km$v_event)))

# incidence rate in unvaxed
V_person_days <- sum(cum_vax_calwk * 5)
U_person_days <- sum((km$u_risk + km$v_risk[1] - cum_vax_calwk) * 5)

1 - (km$v_event[len] / V_person_days) / (km$u_event[len] / U_person_days)

# make plot ---------------------------------------------------------------

# helper function
specd <- function(x, k) trimws(format(round(x, k), nsmall=k))

# create data.frame
df <- data.frame(
  time = seq(as.Date('2022/07/31'), as.Date('2022/09/09'), by = "5 days"),
  null = V_cum_null_calwk,
  vaccinated = km$v_event,
  unvaccinated = km$u_event
)

# reshape to long
df_long <- 
  pivot_longer(df, -time) %>%
  mutate(
    name = factor(
      name, 
      levels = c("unvaccinated", "vaccinated", "null")
    )
  )

# plot of cumulative cases over time
p1 <- ggplot(
  df_long %>% filter(name != "vaccinated"), 
  aes(x = time, y = value, color = name)) +
  geom_point(
    size = 2,
    data = df_long %>% 
      filter(name == "vaccinated" & time == "2022-09-09")
  ) +
  geom_step(aes(linetype = name == "null"), size = 1.2) +
  geom_errorbar(
    aes(
      x = time + 4,
      # xend = time + 1,
      y = NULL,
      ymin = vaccinated,
      ymax = null
    ),
    data = df %>% filter(time == "2022-09-09"),
    color = "grey70"
  ) +
  geom_errorbar(
    aes(
      x = time + 8,
      # xend = time + 1,
      y = NULL,
      ymin = vaccinated,
      ymax = unvaccinated
    ),
    data = df %>% filter(time == "2022-09-09"),
    color = "grey70"
  ) +
  annotate(
    geom = "text",
    x = as.Date('2022/09/09') + 4,
    y = 2,
    label = bquote(VE %~~% 72*"%"),
    hjust = "center",
    color = "grey70"
    # parse = TRUE
  ) +
  annotate(
    geom = "text",
    x = as.Date('2022/09/09') + 8,
    y = 16,
    label = bquote(VE == 84*"%"),
    hjust = "center",
    color = "grey70"
    # parse = TRUE
  ) +
  theme_classic(base_size = 14) +
  scale_color_manual(
    name = "",
    values = c("#FF0000", "#ADD8E6", "#333333"),
    labels = c("unvaccinated", "vaccinated", "null (unvaccinated rate x vaccinated enrollment)"),
    drop = FALSE
  ) +
  scale_linetype_discrete(guide = "none") +
  #scale_x_continuous(breaks = seq(0, 7 * nrow(foo) - 1, 7)) +
  scale_x_date(
    limits = c(as.Date('2022/07/31'), as.Date('2022/09/09') + 10),
    breaks = seq(as.Date('2022/07/31'), as.Date('2022/09/09'), by = "5 days"),
    labels = function(x) str_wrap(format(x, format = "%b %d"), width = 4)
  ) +
  labs(
    x = NULL,
    y = "Cumulative cases of mpox"
  ) +
  theme(legend.position = "top")


# numbers at risk
p3 <- ggplot(
  pivot_longer(df, -time) %>% filter(name != "vaccinated"), 
  aes(
    x = time, 
    y = factor(name, 
               levels = c("vaccinated", "null", "unvaccinated"),
               labels = c("vaccinated (calc.)", "null (calc.)", "unvaccinated")
    ),
    label = format(round(value, 1), big.mark = ",")
  )
) +
  geom_text(hjust = 'center', size = 3.5) +
  scale_y_discrete(labels = function(x) str_wrap(x, 30)) +
  scale_x_date(
    limits = c(as.Date('2022/07/31'), as.Date('2022/09/09') + 10),
    breaks = seq(as.Date('2022/07/31'), as.Date('2022/09/09'), by = "5 days"),
    labels = function(x) str_wrap(format(x, format = "%b %d"), width = 4)
  ) +
  # scale_x_continuous(breaks = seq(0, 7 * nrow(foo) - 1, 7)) +
  theme_classic(base_size = 14) +
  labs(
    x = NULL,
    y = NULL,
    title = "Number at risk"
  ) +
  theme(axis.text.x = element_blank())

# combine plots in single layout using patchwork
p <- p1 / p3 + plot_layout(heights = c(3, 1))
p[[1]] <- p[[1]] + theme(axis.title.y = element_text(vjust = -15))

ggsave(
  filename = "figure2.png",
  plot = p,
  width = 9,
  height = 6.5,
  units = "in",
  dpi = 400
)



# make supplementary table ------------------------------------------------

tab <- data.frame(
  date = seq(as.Date('2022/07/31'), as.Date('2022/09/09'), by = "5 days"),
  unvaccinated = km$u_risk + km$v_risk[1] - cum_vax_calwk,
  events = km$u_event,
  ir = U_incrate_calwk,
  cumvax = cum_vax_calwk,
  vaccinc0 = V_inc_null_calwk
)

tab %>%
  mutate(
    date = format(date, "%b %d"),
    ir = ir * 1000,
    across(-date, round, 1)
  ) %>%
  flextable() %>%
  save_as_docx(path = "tableS1.docx")
