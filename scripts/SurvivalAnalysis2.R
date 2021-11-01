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

# plot the data
ggplot(n2full) +
  geom_line(aes(x=day, y=n_alive))

##################################################
## Reshape Iris data
##################################################
dat_iris <- read_csv("data/iris_lifespan_data.csv")
glimpse(dat_iris)

# need data and fraction alive for each day
proc_dat_iris <- dat_iris %>%
  dplyr::group_by(day, strain) %>%
  dplyr::mutate(a = sum(num_alive_IS),
                tot = 70,
                c = sum(num_censor_IS)) %>%
  dplyr::ungroup() %>%
  dplyr::distinct(day, strain, .keep_all = T) %>%
  dplyr::group_by(strain) %>%
  dplyr::arrange(day) %>%
  dplyr::mutate(cum_c = cumsum(c)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(frac_c = cum_c/tot,
                frac_retained = (1-frac_c)) %>%
  tibble::add_row(day = c(0,0,0), strain = c("CX11264", "QX1793", "QX1794"), frac_retained = c(1,1,1)) 

# plot cumulative censoring
p <- ggplot(proc_dat_iris) +
  geom_step(aes(x=day, y=frac_retained, color = strain)) +
  theme_bw() +
  xlim(0,15) +
  ylim(0,1) +
  labs(x = "days post L4", y = "Fraction not censored")

cowplot::ggsave2(p, filename = "plots/20210810_fraction_not_censored.png", width = 5, height = 5)
########################################
# Example data from https://stackoverflow.com/questions/63178033/how-do-i-group-a-kaplan-meier-and-a-geom-line-under-one-legend-with-linetype
########################################
treatment <- data.frame(
  treatment = c(1, 1, 1, 1, 1, 1), 
  t = c(5.525, 1.9493, 4.9473, 5.9466, 1.5797, 0.5038), 
  event = c(1, 1, 1, 0, 1, 1)
)

tempsurv1 <- c(1.0000000, 0.9129731, 0.8337045, 0.7614860, 0.6956758, 0.6356917, 0.50, 0.43, 0.37)
tempsurv2 <- c(1.0000000, 0.9324888, 0.8671987, 0.8042297, 0.7436717, 0.6856045, 0.6300962, 0.5772029, 0.5269681)
x<- c(0:8)
temp <- data.frame(x, tempsurv1, tempsurv2)

temp<- temp %>% 
  gather(key, value, -c(x))

f2 <- survfit(Surv(t, event)~1, data=treatment)
f2 <- ggsurvplot(f2, legend="right") 
f2 <- f2$plot +  geom_line(data = temp, aes(x=x, y=value, group=key, linetype=key))+
  theme(legend.position="bottom");f2
