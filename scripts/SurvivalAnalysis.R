library(tidyverse)
library(survival)
library(survminer)

# set working directory
setwd(glue::glue("{dirname(rstudioapi::getActiveDocumentContext()$path)}/.."))

# get example of Surv object from survival package
ovarian <- survival::ovarian

# look at structure of ovarian object
str(ovarian)

# learn about it from survival package
?survival::ovarian

# run through example survival analysis
fit1 <- survfit(Surv(futime, fustat) ~ resid.ds, data=ovarian)
print(fit1, rmean= 730)
summary(fit1, times= (0:4)*182.5, scale=365)

fit1_p <- survminer::ggsurvplot(fit1, legend="right", conf.int = T)

# open our example data, set the path to your file path if needed
n2dat <- readr::read_csv("data/lifespan_example_data.csv") %>%
  dplyr::mutate(censor = dplyr::case_when(!is.na(type_of_event) ~ 1,
                                          TRUE ~ 0),
                total_n = sum(num_of_events))

n2dead <- n2dat %>%
  dplyr::filter(censor == 0) %>%
  dplyr::mutate(n_dead = num_of_events) %>%
  dplyr::distinct(day, .keep_all = T)

n2censored <- n2dat %>%
  dplyr::filter(censor == 1) %>%
  dplyr::group_by(day) %>%
  dplyr::mutate(n_censored = sum(num_of_events)) %>%
  dplyr::ungroup()

n2full <- n2dat %>%
  dplyr::distinct(day) %>%
  dplyr::left_join(n2dead) %>%
  dplyr::left_join(dplyr::select(n2censored, day, n_censored)) %>%
  dplyr::distinct(day, .keep_all = T) %>%
  dplyr::select(day, strain, total_n:n_censored) %>%
  dplyr::mutate(n_censored = ifelse(is.na(n_censored), 0, n_censored),
                cum_dead = cumsum(n_dead),
                cum_censor = cumsum(n_censored),
                n_alive = total_n - cum_dead)

# plot the example data
ggplot(n2full) +
  geom_line(aes(x=day, y=n_alive))

# Use our data structure to try to run survival analysis
toydat <- readr::read_csv("data/Tim_example_lifespan_data.csv")
toyfit <- survfit(Surv(day, event2) ~ 1, data=toydat)
summary(toyfit)

# look at plot. Why is survival prob 75% at day 11?
toyfit_p <- survminer::ggsurvplot(toyfit, legend="right", conf.int = T)
toyfit_p
