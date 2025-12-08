# probit regression models run with brms
# for manuscript: "Biological, freshwater, and marine drivers of age at maturity in wild Chinook Salmon"
# authors: Jennifer L. Gosselin, Benjamin P. Sandford, Caitlin S. O'Brien and Eric R. Buhle
# last update: 2025-12-07

# packages
library(here)
library(brms)
library(tidyverse)
library(GGally)
library(scales)
library(shinystan)

# Age of return data with covariate data
DATA <- read_csv(here("data", "adult_age_covariate_data.csv"))


# filter, arrange DATA
data.age.length <- DATA %>%
  #bin adult_age for model
  mutate(age_group = case_when(
    adult_age == 3 ~ 3,
    adult_age == 4 ~ 4,
    between(adult_age, 5,6) ~ 5
  )) %>%
  #sort age
  arrange(age_group) %>%
  # set latent variable as ordered factor
  mutate(age_group = factor(age_group, ordered = TRUE)) %>%
  #set reference categories for categorical variables
  mutate(
    pass_type_T_R = factor(pass_type_T_R, levels = c("ROR", "T")),
    rel_age = factor(rel_age, levels = c("ParrAboveLGR", "SmoltAboveLGR", "SmoltAtLGR")),
    #basin = factor(basin, levels = c("Lower Snake", "Salmon", "Clearwater")),
    MPG = factor(MPG, levels = c("Grande Ronde/Imnaha", "Lower Snake River",
                                 "Clearwater River","South Fork Salmon River",
                                 "Lower Salmon River","Middle Fork Salmon River","Upper Salmon River")
    ),
    MPG.bin = factor(MPG.bin, levels = c("Grande Ronde/Imnaha", "Lower Snake","Clearwater and Salmon")
    )
  ) %>%
  mutate(SYage1 = as.character(SYage1)) %>%
  mutate(SYage1_DOY = as.numeric(SYage1_DOY))


# scale data
data.scaled <- data.age.length %>%
  mutate_at(c(16:22), ~ (scale(.) %>% as.vector())) # adjust as needed


# Consider excluding combinations of covariates with correlation coefficients > 0.4
cor(data.scaled[,16:22], use="pairwise.complete.obs")



## fit_DOY_npgo3: main model
fit_DOY_npgo3 <- brm(age_group ~ length * rel_age + MPG + SYage1_DOY + LGR.flow.7d + pass_type_T_R + NPGO.ONDJFM.T.AVG3 + (1 | SYage1),
                     data = data.scaled,
                     family = cumulative("probit"),
                     prior = c(prior(normal(0, 5), class = Intercept),
                               prior(normal(0, 5), class = b)),
                     chains = 3,
                     cores = 3,
                     iter = 4000, warmup = 1000, thin=1, seed=1234,
                     control = list(adapt_delta = 0.995, max_treedepth = 15),
                     init = 0,
                     # file = here("results/models.local", "fit_DOY_npgo3")
)
summary(fit_DOY_npgo3, prob=c(0.90))
# launch_shinystan(fit_DOY_npgo3)
# prior_summary(fit_DOY_npgo3)
# plot(fit_DOY_npgo3)
# conditional_effects(fit_DOY_npgo3, categorical = T, prob = .95)
# posterior_summary(fit_DOY_npgo3)


ranef.df <- ranef(fit_DOY_npgo3, probs = c(0.05, 0.95))
ranef.df <- round(as.data.frame(ranef.df),4)
write.csv(ranef.df, file=here("results/models.local", "random_effects.csv"))



## fit_TEMP_npgo3: same model with temp instead of DOY
fit_TEMP_npgo3 <- brm(age_group ~ length * rel_age + MPG + LGR.temp.7d + LGR.flow.7d + pass_type_T_R + NPGO.ONDJFM.T.AVG3 + (1 | SYage1),
                      data = data.scaled,
                      family = cumulative("probit"),
                      chains = 3,
                      cores = 3,
                      iter = 4000, warmup = 1000, thin=1, seed=1234,
                      control = list(adapt_delta = 0.995, max_treedepth = 15),
                      init = 0,
                      # file = here("results/models.local", "fit_TEMP_npgo3")
)
summary(fit_TEMP_npgo3, prob=c(0.90))
# launch_shinystan(fit_TEMP_npgo3)



# Subset model by grouping covariates:

## fit_fish
fit_fish <- brm(age_group ~ length*rel_age + MPG + SYage1_DOY+ (1 | SYage1),
                data = data.scaled,
                family = cumulative("probit"),
                chains = 3,
                cores = 3,
                iter = 4000, warmup = 1000, thin=1, seed=1234,
                control = list(adapt_delta = 0.995, max_treedepth = 15),
                init = 0,
                # file = here("results/models.local", "fit_fish")
)
summary(fit_fish, prob=c(0.90))



## fit_freshwater
fit_freshwater1<- brm(age_group ~ LGR.temp.7d + LGR.flow.7d + pass_type_T_R + (1 | SYage1),
                      data = data.scaled,
                      family = cumulative("probit"),
                      chains = 3,
                      cores = 3,
                      iter = 4000, warmup = 1000, thin=1, seed=1234,
                      control = list(adapt_delta = 0.995, max_treedepth = 15),
                      init = 0,
                      # file = here("results/models.local", "fit_freshwater")
)
summary(fit_freshwater1, prob=c(0.90))

fit_freshwater2<- brm(age_group ~ LGR.flow.7d + pass_type_T_R + (1 | SYage1),
                      data = data.scaled,
                      family = cumulative("probit"),
                      chains = 3,
                      cores = 3,
                      iter = 4000, warmup = 1000, thin=1, seed=1234,
                      control = list(adapt_delta = 0.995, max_treedepth = 15),
                      init = 0,
                      # file = here("results/models.local", "fit_freshwater2")
)
summary(fit_freshwater2, prob=c(0.90))



## fit_marine
fit_marine<- brm(age_group ~ NPGO.ONDJFM.T.AVG3 + (1 | SYage1),
                 data = data.scaled,
                 family = cumulative("probit"),
                 chains = 3,
                 cores = 3,
                 iter = 4000, warmup = 1000, thin=1, seed=1234,
                 control = list(adapt_delta = 0.995, max_treedepth = 15),
                 init = 0,
                 # file = here("results/models.local", "fit_marine")
)
summary(fit_marine, prob=c(0.90))



## compare
loo.results<-loo(fit_DOY_npgo3, fit_TEMP_npgo3, fit_fish, fit_freshwater1, fit_freshwater2, fit_marine)


