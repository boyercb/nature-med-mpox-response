Replication files for:

Boyer, C. \& Lipsitch, M. (2023) "Comment on: Real-world effectiveness of a single dose of mpox vaccine in males."

## Abstract
Wolff Sagy et al. report that a single dose of MVA-BN mpox vaccine was 86% effective in preventing incident mpox. In their cohort study, all eligible men were enrolled at the start of a vaccination campaign on July 31, 2022 and followed until they tested positive for mpox or December 25, 2022, whichever was earlier. To avoid the potential for immortal time bias, participants vaccinated after July 31 were considered unvaccinated from July 31 through their date of vaccination and then considered vaccinated from their date of vaccination, which was set as their new Day 0. Incidence of mpox was compared using time-dependent Cox regression on this time scale, with typically different calendar days being Day 0 for vaccinated and unvaccinated persons. In this communication, we show that while their method eliminates immortal time, it does not account for confounding by calendar time which is particularly important during an outbreak where incidence rates vary drastically over time. We evaluate the potential magnitude of this bias through simulation and a re-analysis of initial data reported by Wolff Sagy et al. and conclude with a discussion of how it can be resolved in future studies.

## How to run
- `sim.R` reproduces the simulation of a vaccination campaign during an outbreak in Figure 1.
- `reanalysis.R` reproduces the re-analysis of Wolff Sagy et al. adjusting for calendar time in Figure 2 and Table S1.

## Inputs
- `KM.csv` extracted data from Kaplan Meier table in Figure 1 of Arbel, R., Sagy, Y. W., Zucker, R., Arieh, N. G., Markovits, H., Abu-Ahmad, W., Battat, E., Ramot, N., Carmeli, G., Mark-Amir, A., Wagner-Kolasko, G., Duskin-Bitan, H., Yaron, S., Peretz, A., Hammerman, A., Lavie, G., & Netzer, D. (2022). Effectiveness of a single-dose Modified Vaccinia Ankara in Human Monkeypox: An observational study [Preprint]. Research Square. https://doi.org/[https://doi.org/10.21203/rs.3.rs-1976861/v2.
- `seirv.R` helper functions defining the SEIR model used in the simulation.

## Outputs
- `figure1.png`
- `figure2.png`
- `TableS1.docx`

