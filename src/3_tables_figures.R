# Figures and Tables
# for manuscript: "Biological, freshwater, and marine drivers of age at maturity in wild Chinook Salmon"
# authors: Jennifer L. Gosselin, Benjamin P. Sandford, Caitlin S. O'Brien and Eric R. Buhle
# last update: 2025-12-07

# Packages
library(here)
library(tidyverse)
library(ggplot2)
library(ggridges) # for forklength joyplot
library(pals) # check color contrasts
library(patchwork) # align plots and inset figures
library(brms) # import model
library(bayesplot) # for model fit and ppc_check
library(data.table) # fread() for covariate data
library(tidybayes) # gathering draws and epred
library(ggdist) # plotting stat_half_eye color ramp
library(flextable) # use for all table design word or pdf- https://ardata-fr.github.io/flextable-book/index.html#create-a-flextable
library(ftExtra) # use to split headers better then separate_headers() https://cran.r-project.org/web/packages/ftExtra/vignettes/transform-headers.html
library(modelr) # seq_range
library(readxl)
library(ggstream)
library(showtext)
library(ggtext)


# color palettes
# https://jfly.uni-koeln.de/color/

pal_all <- c("#000000", "#d1ae41", "#56B4E9", "#D55E00", "#0072B2", "#CC79A7", "#009E73", "#F0D449", "#E69F00")

pal_age <- c("#444444", "#83dcdd", "#e04555", "#000000")

pal_cov <- c("#d1ae41", "#009E73", "#2b92b4") # gray: #555d66

# check contrasts
# pal.safe(pal_cov)


# set function for SYage1 shown on graphs
every_nth <- function(n) {
  return(function(x) {
    x[c(TRUE, rep(FALSE, n - 1))]
  })
}


# import data
DATA <- read_csv(here("data", "adult_age_covariate_data.csv"))

# scale data
data.scaled <- DATA %>%
  mutate_at(c(16:21), ~ (scale(.) %>% as.vector())) %>%
  # add each OA for figure 2 & 3
  mutate(OA.true = case_when(
    adult_age == 3 ~ "1-Ocean",
    adult_age == 4 ~ "2-Ocean",
    adult_age == 5 ~ "3-Ocean",
    adult_age == 6 ~ "4-Ocean",
    .default = "other"
  )) %>%
  mutate(age_group = case_when(
    adult_age == 3 ~ 3,
    adult_age == 4 ~ 4,
    between(adult_age, 5, 6) ~ 5
  )) %>%
  # sort age
  arrange(age_group) %>%
  # set latent variable as ordered factor
  mutate(age_group = factor(age_group, ordered = TRUE))


# import model
m <- readRDS(here("results/models.local", "fit_DOY_npgo3.rds"))



## Figure 2. Counts and proportion by age
d.return <- data.scaled %>%
  mutate_if(sapply(data.scaled, is.character), as.factor) %>%
  group_by(pass_type_T_R, SYage1, OceanAge, OA.true) %>%
  summarise(n = n())

d.return.yr <- d.return %>%
  group_by(SYage1, OceanAge) %>%
  summarize(count = sum(n))


p.line.count <- ggplot(d.return.yr, aes(x = SYage1, y = count, col=OceanAge)) +
  geom_line() +
  geom_point(size = 3, shape = 15) +
  labs(
    x = "", y = "Count of adult returns",
    color = NULL,
    fill = NULL,
    title = NULL,
    caption = ""
  ) +
  scale_color_manual(values = c(pal_age)) +
  scale_x_continuous(breaks = 1998:2020, limits = c(1998,2020), expand = c(0.02,0.02)) +
  ggtitle("(a) Count by age") +
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.line.y = element_line(linewidth = 0.2),
    legend.position = "top",
    legend.key.size = unit(0.5, "cm"),
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(margin = margin(t = 0, r = 0, b = 0, l = 0), size = 12),
    axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0), size = 12),
    plot.title = element_text(size = 12, vjust = -7),
    panel.grid.minor = element_blank()
  )

p.bar.prop <- ggplot(d.return, aes(x = as.factor(SYage1), y = n, fill = forcats::fct_rev(OceanAge))) +
  geom_bar(position = position_fill(), stat = "identity", width = .7) +
  labs(
    x = "Smolt migration year", y = "Proportion",
    color = "",
    fill = "",
    title = NULL,
    subtitle = NULL,
    caption = ""
  ) +
  scale_x_discrete(breaks = 1998:2020) +
  ggtitle("(b) Proportion by age") +
  scale_fill_manual(
    values = pal_age,
    breaks = c("1-Ocean", "2-Ocean", "3-Ocean or greater")
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 12),
    axis.title.x = element_text(vjust = -1.5, size = 12),
    panel.background = element_blank()
  )

p_Fig2 <- p.line.count + p.bar.prop + plot_layout(nrow=2)

# pdf(here("results", "figures", "manuscript", "Fig2 Count by Age.pdf"), width=7, height=5, onefile=TRUE)
p_Fig2
# dev.off()



## Figure 3. hypothetical distr. and tao cutpoints
# load(here("results/models.local", "fit_DOY_npgo3.rds"))
m <- fit_DOY_npgo3

draws <- as_draws_df(m) %>%
  select(.draw, `b_Intercept[1]`:`b_Intercept[2]`)

# summarize draws (int 1 and 2 only)
taos <- m %>%
  gather_draws(`b_Intercept[1]`, `b_Intercept[2]`) %>%
  summarise_draws()
taos

# only select medians with 5 to 95 CI
med_taos <- taos %>%
  select(1, 4, 7:8) %>%
  pivot_wider(names_from = .variable, values_from = c("median", "q5", "q95"))
med_taos

# plot with verticle lines to represent median, with 5 to 95 CI
tibble(x = seq(from = -3, to = 3, by = .01)) %>%
  mutate(d = dnorm(x)) %>%
  ggplot(aes(x = x, ymin = 0, ymax = d)) +
  geom_ribbon(fill = "darkgrey") +
  geom_vline(xintercept = as.numeric(med_taos[1, 1:2]), color = "white", linetype = 1, linewidth = 1) +
  geom_vline(xintercept = as.numeric(med_taos[1, 3:6]), color = "white", linetype = 2, linewidth = .3) +
  labs(y = "Probability") +
  scale_x_continuous(NULL,
                     breaks = c(-3, -2, -1, 0, 1, 2, 3),
                     labels = c(-3, -2, -1, 0, 1, 2, 3)
  ) +
  annotate(geom="text", x=-1.38, y=.18, label=expression(paste(tau[1]))) +
  annotate(geom="text", x=.935, y=.29, label=expression(paste(tau[2]))) +
  # scale_x_continuous(NULL,
  #                    breaks = c(as.numeric(med[1, 1:2])),
  #                    labels = c(parse(text = str_c("tao[", 1:2, "]")))
  # ) +
  # ggtitle("Standard normal distribution underlying the ordinal Y data:",
  #         subtitle = "The solid vertical lines mark the posterior means for the thresholds. \nThe dashed vertical lines represent the 95 CI") +
  theme_minimal() +
  theme(
    strip.text = element_blank(),
    axis.line.x = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.line.y = element_line(linewidth = 0.3),
    axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
    axis.title.y = element_text(margin = margin(t = 0, r = 8, b = 0, l = 0))
  ) +
  coord_cartesian(xlim = c(-3, 3))





## Figure 4. Beta effects
# pdf(here("results", "figures", "manuscript", "Fig3 Beta params.pdf"), width=7, height=9, onefile=TRUE)
m %>%
  gather_draws(
    b_MPGLowerSnakeRiver, b_MPGClearwaterRiver,
    b_MPGSouthForkSalmonRiver, b_MPGLowerSalmonRiver,
    b_MPGMiddleForkSalmonRiver, b_MPGUpperSalmonRiver,
    b_length, b_rel_ageSmoltAboveLGR, b_rel_ageSmoltAtLGR,
    `b_length:rel_ageSmoltAboveLGR`, `b_length:rel_ageSmoltAtLGR`,
    b_SYage1_DOY, b_pass_type_T_RT, b_LGR.flow.7d, b_NPGO.ONDJFM.T.AVG3
  ) %>%
  mutate(subgroup = if_else(.variable %in% c("b_NPGO.ONDJFM.T.AVG3"), "Marine",
                            if_else(.variable %in% c("b_pass_type_T_RT", "b_LGR.flow.7d"), "Freshwater", "Fish")
  )) %>%
  mutate(.variable = factor(.variable,
                            levels = c(
                              "b_NPGO.ONDJFM.T.AVG3",
                              "b_LGR.flow.7d",
                              "b_pass_type_T_RT",
                              "b_SYage1_DOY",
                              "b_MPGUpperSalmonRiver", "b_MPGMiddleForkSalmonRiver", "b_MPGLowerSalmonRiver",
                              "b_MPGSouthForkSalmonRiver", "b_MPGClearwaterRiver", "b_MPGLowerSnakeRiver",
                              "b_length:rel_ageSmoltAtLGR", "b_length:rel_ageSmoltAboveLGR",
                              "b_rel_ageSmoltAtLGR", "b_rel_ageSmoltAboveLGR",
                              "b_length"
                            ),
                            labels = c(
                              "NPGO",
                              "Flow",
                              "Transported",
                              "Day of Year",
                              "Upper Salmon",
                              "MF Salmon",
                              "Lower Salmon",
                              "SF Salmon",
                              "Clearwater",
                              "Lower Snake",
                              "Length x \nSmolt&AtLGR",
                              "Length x \nSmolt&AboveLGR",
                              "Smolt&AtLGR",
                              "Smolt&AboveLGR",
                              "Length"
                            )
  )) %>%
  ggplot(aes(y = .variable, x = .value)) + # reorder(variable, desc(.value))
  geom_vline(
    xintercept = 0, linetype = "dashed",
    color = "darkgrey"
  ) +
  stat_halfeye(
    aes(fill = subgroup, fill_ramp = after_stat(level)),
    .width = c(.50, .80, .90), scale = .7, # adjust based on CI wanted--poster includes .85 &.95 only
    # NOTE: we use position = "dodgejust" (a dodge that respects the
    # justification of intervals relative to slabs) instead of
    # position = "dodge" here because it ensures the topmost slab does
    # not extend beyond the plot limits
    position = "dodgejust"
  ) +
  # a range from 1 down to 0.2 ensures the fill goes dark to light inside-out
  # and doesn't get all the way down to white (0) on the lightest color
  scale_fill_ramp_discrete(na.translate = FALSE) +
  scale_x_continuous(breaks = seq(-1, 0.5, by = .5), limits = c(-1, 0.5)) +
  scale_y_discrete(expand = c(-.05, .9)) +
  scale_fill_manual(
    breaks = c("Fish", "Freshwater", "Marine"),
    values = pal_cov, guide = NULL
  ) +
  geom_hline(yintercept = 15.7, color = pal_cov[1], alpha = 0.3, linewidth = 3.5) +
  annotate("segment", x = -.98, y = 15.6, xend = -.98, yend = 3.83, color = pal_cov[1], alpha = .4, linewidth = 2.3) +
  annotate("text", y = 15.7, x = -0.98, label = "  Biological / behavioural covariates", color = "black", hjust = 0, size = 3.1) +

  geom_hline(yintercept = 3.7, color = pal_cov[2], alpha = 0.2, linewidth = 3.5) +
  annotate("segment", x = -.98, y = 3.6, xend = -.98, yend = 1.8, color = pal_cov[2], alpha = .3, linewidth = 2.3) +
  annotate("text", y = 3.7, x = -0.98, label = "  Freshwater covariates", color = "black", hjust = 0, size = 3.1) +

  geom_hline(yintercept = 1.7, color = pal_cov[3], alpha = 0.2, linewidth = 3.5) +
  annotate("segment", x = -.98, y = 1.6, xend = -.98, yend = 0.8, color = pal_cov[3], alpha = .3, linewidth = 2.3) +
  annotate("text", y = 1.7, x = -0.98, label = "  Marine covariates", color = "black", hjust = 0, size = 3.1) +

  annotate("segment", x = -.98, y = 15.1, xend = -.98, yend = 14.9, color = "#a37a29", alpha = .4, linewidth = 2.3) +

  annotate("text", x = -.775, y = 13.65, label = "Tag Age & Location", size = 3, color = "#a37a29") +
  annotate("text", x = -.73, y = 13.35, label = "Reference: Parr&AboveLGR", size = 2.3, color = "#a37a29") +
  annotate("segment", x = -.98, y = 14.3, xend = -.98, yend = 12.7, color = "#a37a29", alpha = .4, linewidth = 2.3) +

  annotate("text", x = -.69, y = 11.65, label = "Length x Tag Age & Location", size = 3, color = "#a37a29") +
  annotate("text", x = -.665, y = 11.35, label = "Reference: Length x Parr&AboveLGR", size = 2.3, color = "#a37a29") +
  annotate("segment", x = -.98, y = 12.3, xend = -.98, yend = 10.7, color = "#a37a29", alpha = .4, linewidth = 2.3) +

  annotate("text", x = -.9, y = 7.65, label = "Origin", size = 3, color = "#a37a29") +
  annotate("text", x = -.69, y = 7.35, label = "Reference: Grande Ronde/Imnaha", size = 2.3, color = "#a37a29") +
  annotate("segment", x = -.98, y = 10.3, xend = -.98, yend = 4.7, color = "#a37a29", alpha = .4, linewidth = 2.3) +
  annotate("segment", x = -.98, y = 4.1, xend = -.98, yend = 3.9, color = "#a37a29", alpha = .4, linewidth = 2.3) +

  annotate("text", x = -.7, y = 3.15, label = "Hydrosystem Passage Type", size = 3, color = pal_cov[2]) +
  annotate("text", x = -.8, y = 2.85, label = "Reference: In-river", size = 2.3, color = pal_cov[2]) +
  annotate("segment", x = -.98, y = 3.3, xend = -.98, yend = 2.7, color = pal_cov[2], alpha = .4, linewidth = 2.3) +
  annotate("segment", x = -.98, y = 2.1, xend = -.98, yend = 1.9, color = pal_cov[2], alpha = .4, linewidth = 2.3) +

  annotate("segment", x = -.98, y = 1.1, xend = -.98, yend = 0.9, color = pal_cov[3], alpha = .4, linewidth = 2.3) +
  labs(
    fill_ramp = "CI", # is this confidence or credible interval? double check
    fill = "Covariates",
    y = expression(bold("Covariate")),
    x = expression(bold(paste("Covariate effect (", beta, ")")))
  ) +
  theme_minimal() +
  theme(
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2, color = "gray"),
    axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
    legend.position = "right",
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 9),
    # axis.text.y = element_text(hjust = .01),
    plot.title.position = "plot"
  )
# dev.off()



## Figure 5. Conditional effects

c_eff <- conditional_effects(m, prob = .9, categorical = TRUE, method = "posterior_epred", re_formula = NA)


# length
FLz_xaxis <- (seq(40, 160, 10) - mean(DATA$length)) / sd(DATA$length)
FL_lab <- as.character(seq(40, 160, 10))

p_length <- plot(c_eff, plot = FALSE)[[1]] +
  scale_color_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater")
  ) +
  scale_fill_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater")
  ) +
  scale_shape_manual(guide = "none") +
  scale_y_continuous(breaks = seq(0, 0.8, by = .2), limits = c(0, 0.8)) +
  labs(x = "Length (mm)", y = "Probability", fill = NULL, color = NULL, shape = NULL, y = NULL) +
  ggtitle("(a) Length") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.line.y = element_line(linewidth = 0.2),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9.5),
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(size = 9),
    plot.margin = unit(c(0.3, 0.3, 0.6, 0.3), "cm")
  ) +
  scale_x_continuous(breaks = FLz_xaxis, labels = FL_lab)


# Tag_Age&Location
p_rel_age <- plot(c_eff, plot = FALSE)[[2]] +
  scale_color_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater"),
    guide = "none"
  ) +
  scale_fill_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater"),
    guide = "none"
  ) +
  scale_shape_manual(
    values = 1,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater"),
    guide = "none"
  ) +
  scale_x_discrete(
    breaks = c("ParrAboveLGR", "SmoltAboveLGR", "SmoltAtLGR"),
    labels = c("Parr&\nAboveLGR", "Smolt&\nAboveLGR", "Smolt&\nAtLGR")
  ) +
  scale_y_continuous(breaks = seq(0, 0.8, by = .2), limits = c(0, 0.8)) +
  labs(x = "Tag age & location", fill = NULL, color = NULL, shape = NULL, y = NULL) + # coord_capped_cart(bottom='both', left = "both", xlim= c(-3,3)) +
  ggtitle("(b) Tag age & location") +
  theme_minimal() +
  theme(
    legend.position = "right",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.line.y = element_line(linewidth = 0.2),
    # axis.text.x = element_text(angle=10),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9.5),
    legend.key.size = unit(0.2, "cm"),
    legend.text = element_text(size = 9),
    plot.margin = unit(c(0.3, 0.3, 0.6, 0.3), "cm")
  )


# length x tag_age&location

c_eff_intrxn_1 <- conditional_effects(
  m,
  effects = "length",
  conditions = data.frame(rel_age = "ParrAboveLGR"),
  categorical = TRUE
)

p_intrxn_1 <- plot(c_eff_intrxn_1, plot = FALSE)[[1]] +
  scale_color_manual(values = c(
    "3" = pal_age[1],
    "4" = pal_age[2],
    "5" = pal_age[3]
  )) +
  scale_fill_manual(values = c(
    "3" = pal_age[1],
    "4" = pal_age[2],
    "5" = pal_age[3]
  )) +
  scale_y_continuous(breaks = seq(0, 0.8, by = .2), limits = c(0, 0.8)) +
  labs(x = "Length (mm)", y = "Probability", fill = NULL, color = NULL, shape = NULL) +
  ylab("Probability") +
  ggtitle("(c) incl. Length x ParrAboveLGR") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.line.y = element_line(linewidth = 0.2),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9.5),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")
  ) +
  scale_x_continuous(breaks = FLz_xaxis, labels = FL_lab)


c_eff_intrxn_2 <- conditional_effects(
  m,
  effects = "length",
  conditions = data.frame(rel_age = "SmoltAboveLGR"),
  categorical = TRUE
)

p_intrxn_2 <- plot(c_eff_intrxn_2, plot = FALSE)[[1]] +
  scale_color_manual(values = c(
    "3" = pal_age[1],
    "4" = pal_age[2],
    "5" = pal_age[3]
  )) +
  scale_fill_manual(values = c(
    "3" = pal_age[1],
    "4" = pal_age[2],
    "5" = pal_age[3]
  )) +
  scale_y_continuous(breaks = seq(0, 0.8, by = .2), limits = c(0, 0.8)) +
  labs(x = "Length (mm)", y=NULL, fill = NULL, color = NULL, shape = NULL) +
  ggtitle("(d) incl. Length x SmoltAboveLGR") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.line.y = element_line(linewidth = 0.2),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9.5),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")
  ) +
  scale_x_continuous(breaks = FLz_xaxis, labels = FL_lab)



c_eff_intrxn_3 <- conditional_effects(
  m,
  effects = "length",
  conditions = data.frame(rel_age = "SmoltAtLGR"),
  categorical = TRUE
)

p_intrxn_3 <- plot(c_eff_intrxn_3, plot = FALSE)[[1]] +
  scale_color_manual(values = c(
    "3" = pal_age[1],
    "4" = pal_age[2],
    "5" = pal_age[3]
  )) +
  scale_fill_manual(values = c(
    "3" = pal_age[1],
    "4" = pal_age[2],
    "5" = pal_age[3]
  )) +
  scale_y_continuous(breaks = seq(0, 0.8, by = .2), limits = c(0, 0.8)) +
  labs(x = "Length (mm)", y=NULL, fill = NULL, color = NULL, shape = NULL) +
  ggtitle("(e) incl. Length x SmoltAtLGR") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.line.y = element_line(linewidth = 0.2),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9.5),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")
  ) +
  scale_x_continuous(breaks = FLz_xaxis, labels = FL_lab)



# Origin
p_origin <- plot(c_eff, plot = FALSE)[[3]] +
  scale_color_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater"),
    guide = "none"
  ) +
  scale_fill_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater"),
    guide = "none"
  ) +
  scale_shape_manual(
    values = 1,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater"),
    guide = "none"
  ) +
  labs(x = "Origin", y = "Probability", fill = NULL, color = NULL, shape = NULL, y = NULL) +
  scale_x_discrete(
    breaks = c(
      "Grande Ronde/Imnaha", "Lower Snake River",
      "Clearwater River", "South Fork Salmon River",
      "Lower Salmon River", "Middle Fork Salmon River",
      "Upper Salmon River"
    ),
    labels = c(
      "Grande Ronde/\nImnaha", "L Snake",
      "Clearwater", "SF Salmon",
      "L Salmon", "MF Salmon",
      "U Salmon"
    )
  ) +
  scale_y_continuous(breaks = seq(0, 0.8, by = .2), limits = c(0, 0.8)) +
  ggtitle("(f) Fish tag & release origin") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.line.y = element_line(linewidth = 0.2),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9.5),
    plot.margin = unit(c(0.3, 0.3, 0.6, 0.3), "cm")
  )
# axis.text.x = element_text(angle = 270,vjust=.5, hjust = 1)


# DOY
DOYz_xaxis <- (seq(80, 160, 10) - mean(DATA$SYage1_DOY)) / sd(DATA$SYage1_DOY)
DOY_lab <- as.character(seq(80, 160, 10))

p_DOY <- plot(c_eff, plot = FALSE)[[4]] +
  scale_color_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater")
  ) +
  scale_fill_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater")
  ) +
  scale_shape_manual(guide = "none") +
  labs(x = "Day of year", y = "Probability", fill = NULL, color = NULL, shape = NULL, y = NULL) +
  scale_y_continuous(breaks = seq(0, 0.8, by = .2), limits = c(0, 0.8)) +
  # scale_x_continuous(breaks = seq(-3,3,1))+
  ggtitle("(g) Day of year at LGR") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.line.y = element_line(linewidth = 0.2),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9.5),
    plot.margin = unit(c(0.3, 0.3, 0.6, 0.3), "cm")
  ) +
  scale_x_continuous(breaks = DOYz_xaxis, labels = DOY_lab)


# pass_type
p_pass_type <- plot(c_eff, plot = FALSE)[[6]] +
  scale_color_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater"),
    guide = "none"
  ) +
  scale_fill_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater"),
    guide = "none"
  ) +
  scale_shape_manual(
    values = 1,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater"),
    guide = "none"
  ) +
  scale_x_discrete(labels = c("In-river", "Transported")) +
  scale_y_continuous(breaks = seq(0, 0.8, by = .2), limits = c(0, 0.8)) +
  labs(x = "Passage type", fill = NULL, color = NULL, shape = NULL, y = NULL) +
  ggtitle("(h) Hydrosystem passage type") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.line.y = element_line(linewidth = 0.2),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9.5),
    plot.margin = unit(c(0.3, 0.3, 0.6, 0.3), "cm")
  )



# flow
flowz_xaxis <- (seq(0, 300, 20) - mean(DATA$LGR.flow.7d)) / sd(DATA$LGR.flow.7d)
flow_lab <- as.character(seq(0, 300, 20))

p_flow <- plot(c_eff, plot = FALSE)[[5]] +
  scale_color_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater")
  ) +
  scale_fill_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater")
  ) +
  scale_shape_manual(guide = "none") +
  scale_y_continuous(breaks = seq(0, 0.8, by = .2), limits = c(0, 0.8)) +
  labs(x = "Flow (kcfs)", y = "Probability", fill = NULL, color = NULL, shape = NULL, y = NULL) +
  ggtitle("(i) Flow at LGR") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.line.y = element_line(linewidth = 0.2),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9.5),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")
  ) +
  scale_x_continuous(breaks = flowz_xaxis, labels = flow_lab)


# NPGO
npgoz_xaxis <- (seq(-2.5, 2.5, 0.5) - mean(DATA$NPGO.ONDJFM.T)) / sd(DATA$NPGO.ONDJFM.T)
npgo_lab <- as.character(seq(-2.5, 2.5, 0.5))

p_NPGO <- plot(c_eff, plot = FALSE)[[7]] +
  scale_color_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater")
  ) +
  scale_fill_manual(
    values = pal_age,
    breaks = c("3", "4", "5"),
    labels = c("1-Ocean", "2-Ocean", "3-Ocean or greater")
  ) +
  scale_shape_manual(guide = "none") +
  scale_y_continuous(breaks = seq(0, 0.8, by = .2), limits = c(0, 0.8)) +
  labs(x = "NPGO index", fill = NULL, color = NULL, shape = NULL, y = NULL) +
  ggtitle("(j) NPGO index, 3-year ONDJFM mean") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.line.y = element_line(linewidth = 0.2),
    plot.title = element_text(size = 10),
    axis.title = element_text(size = 9.5),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm")
  ) +
  scale_x_continuous(breaks = npgoz_xaxis, labels = npgo_lab)

# Fig 5 a-e
patch1 <- (p_length + p_rel_age) / (p_intrxn_1 + p_intrxn_2 + p_intrxn_3)
# pdf(here("results", "figures", "manuscript", "Fig5a_e length effects.pdf"), width=8, height=7, onefile=TRUE)
patch1
# dev.off()

# Fig 5 f-j
patch2 <- p_origin / (p_DOY + p_pass_type ) / (p_flow + p_NPGO)
# pdf(here("results", "figures", "manuscript", "Fig5f_j cond effects.pdf"), width=7, height=9, onefile=TRUE)
patch2
# dev.off()



# Figure 6a.
# pdf(here("results", "figures", "manuscript", "Fig6a Obs vs Pred proportions no RE.pdf"), width=7, height=9, onefile=TRUE)
y_int <- as.integer(data.scaled$age_group)

# set yrep
yrep_int <- posterior_predict(m, newdata = NULL, re_formula = NA)

color_scheme_set(scheme = "mix-darkgray-gray")


ppc_bars_grouped(
  y = y_int,
  yrep = yrep_int,
  group = data.scaled$SYage1,
  freq = FALSE, # adjusts from TRUE = count and FALSE = proportion
  prob = 0.95,
  fatten = 1,
  size = 1.5,
  facet_args = list(
    ncol = 3,
    scales = "fixed"
  ) # "fixed" or "free_y"
) + scale_x_continuous(
  breaks = c(1, 2, 3), # check back on this
  labels = c("1-Ocean", "2-Ocean", "3-Ocean")
) +
  labs(x = "Ocean Age") + theme(text = element_text(size = 14)) +
  theme_tidybayes() +
  theme(
    legend.position = c(.85, .05),
    legend.key.size = unit(.5, "cm"),
    legend.spacing.y = unit(.02, "cm"),
    legend.background = element_blank()
  )
# dev.off()



## Figure 6b.
# pdf(here("results", "figures", "manuscript", "Fig6b Obs vs Pred counts no RE.pdf"), width=7, height=9, onefile=TRUE)

ppc_bars_grouped(
  y = y_int,
  yrep = yrep_int,
  group = data.scaled$SYage1,
  freq = TRUE, # adjusts from TRUE = count and FALSE = proportion
  prob = 0.9,
  fatten = 1,
  size = 1.5,
  facet_args = list(
    ncol = 3,
    scales = "free_y"
  ) # "fixed" or "free_y"
) + scale_x_continuous(
  breaks = c(1, 2, 3),
  labels = c("1-Ocean", "2-Ocean", "3-Ocean")
) +
  labs(x = "Ocean Age", y = "Count of adult returns") + theme(text = element_text(size = 14)) +
  theme_tidybayes() +
  theme(
    legend.position = "none",
    strip.text.x = element_text(margin = margin(0.05, 0, 0.05, 0, "cm")),
    strip.text = element_text(size = 8),
    panel.spacing = unit(1, "lines")
  )
# dev.off()




# Supplementary Tables


# Table S1.
## Number of adult returns by origin and age
tb1 <- DATA %>%
  group_by(MPG, adult_age) %>%
  summarise(count = n()) %>%
  pivot_wider(
    names_from = adult_age, values_from = count,
    values_fill = 0, names_sort = TRUE,
    names_expand = TRUE
  )

tb2 <- DATA %>%
  group_by(MPG) %>%
  summarise(n = sum(n()))

tj <- full_join(tb1, tb2, by = c("MPG")) %>%
  mutate(MPG = factor(MPG,
                      levels = c(
                        "Grande Ronde/Imnaha", "Lower Snake River", "Clearwater River",
                        "Lower Salmon River", "South Fork Salmon River", "Middle Fork Salmon River",
                        "Upper Salmon River"
                      ),
                      labels = c(
                        "Grande Ronde Imnaha River", "Lower Snake River", "Clearwater River",
                        "Lower Salmon River", "South Fork Salmon River", "Middle Fork Salmon River",
                        "Upper Salmon River"
                      )
  ))

f.tj <- flextable(tj)

f.adj <- f.tj %>%
  labelizor(
    part = "header",
    labels = c("MPG" = "Origin")
  ) %>%
  add_header_row(
    values = c("", "Age", ""),
    colwidths = c(1, 4, 1), top = TRUE
  ) %>%
  bold(j = c(6), bold = TRUE, part = "all") %>%
  align(align = "left", part = "all") %>%
  font(fontname = "Calibri") %>%
  hline(part = "header") %>%
  autofit()

f.final <- set_table_properties(
  f.adj,
  layout = "autofit",
  width = 0,
  align = "center",
  word_title = NULL,
  word_description = NULL
)
f.final



# Table S2.
## Number of adult returns by tag_age&location and ocean age per smolt migration year
tb1 <- DATA %>%
  group_by(SYage1, rel_age, adult_age) %>%
  summarise(count = n()) %>%
  pivot_wider(
    names_from = c(rel_age, adult_age), values_from = count,
    values_fill = 0, names_sort = TRUE,
    names_expand = TRUE
  )

tb2 <- DATA %>%
  group_by(SYage1, rel_age) %>%
  summarise(sum = sum(n())) %>%
  pivot_wider(
    names_from = rel_age, values_from = sum,
    values_fill = 0, names_glue = "{rel_age}_n",
    names_sort = TRUE, names_expand = TRUE
  )

tj <- full_join(tb1, tb2, by = "SYage1") %>%
  select(SYage1, contains("parr"), contains("SmoltAbove"), contains("SmoltAt"))



f.tj <- flextable(tj) %>%
  span_header()


f.adj <- labelizor(
  x = f.tj,
  part = "header",
  labels = c(
    "SYage1" = "Juvenile\nmigration\nyear",
    "ParrAboveLGR" = "Parr&AboveLGR",
    "SmoltAboveLGR" = "Smolt&AboveLGR",
    "SmoltAtLGR" = "Smolt&AtLGR"
  )
)

f.adj <- colformat_num(f.adj,
                       j = 1,
                       big.mark = ""
) %>%
  bold(j = c(6, 11, 16), bold = TRUE, part = "all") %>%
  align(align = "left", part = "all") %>%
  font(fontname = "Calibri") %>%
  hline(part = "header") %>%
  autofit()


f.final <- set_table_properties(
  f.adj,
  layout = "autofit",
  width = 0,
  align = "center",
  word_title = NULL,
  word_description = NULL
)
f.final




# Table S3.
## Number of adult returns by passage type and age
tb1 <- DATA %>%
  group_by(SYage1, pass_type_T_R, adult_age) %>%
  summarise(count = n()) %>%
  pivot_wider(names_from = c(pass_type_T_R, adult_age), values_from = count) %>%
  select(SYage1, "ROR_3", "ROR_4", "ROR_5", "ROR_6", "T_3", "T_4", "T_5", "T_6") %>%
  replace(is.na(.), 0)

tb2 <- DATA %>%
  group_by(SYage1, pass_type_T_R) %>%
  summarise(sum = sum(n())) %>%
  pivot_wider(names_from = pass_type_T_R, values_from = sum)

tj <- full_join(tb1, tb2, by = "SYage1") %>%
  select(SYage1, "ROR_3", "ROR_4", "ROR_5", "ROR_6", "ROR", "T_3", "T_4", "T_5", "T_6", "T") %>%
  rename(
    "Juvenile\nmigration\nyear" = SYage1,
    "ROR_n" = ROR,
    "T_n" = T
  )

f.tj <- flextable(tj) %>%
  span_header() %>%
  labelizor(
    part = "header",
    labels = c(
      "ROR" = "In-River",
      "T" = "Transported"
    )
  ) %>%
  align(align = "left", part = "all") %>%
  font(fontname = "Calibri") %>%
  colformat_num(
    j = "Juvenile\nmigration\nyear",
    big.mark = ""
  ) %>%
  bold(j = c(6, 11), bold = TRUE, part = "all") %>%
  autofit()

f.tj


# Supplementary Figures

# Figure S1
## Passage type
df.sum <- data.scaled %>%
  group_by(SYage1) %>%
  summarise(n = n())

df.t <- data.scaled %>%
  group_by(SYage1, pass_type_T_R) %>%
  filter(pass_type_T_R == "T") %>%
  summarise(nT = n())

df.pctT <- inner_join(df.t, df.sum, by = "SYage1") %>%
  mutate(pctT = scales::percent(nT / n)) %>%
  mutate(pctR = scales::percent(1-(nT/n)))

# pdf(here("results", "figures", "manuscript", "FigS1 Passtype.pdf"), width=8, height=4, onefile=TRUE)
ggplot(d.return, aes(x = as.factor(SYage1))) +
  geom_bar(aes(fill = pass_type_T_R, y = n), position = "stack", stat = "identity", width = .6 ) +
  geom_text(data = df.pctT, aes(label = pctT, y = n, group = SYage1), position = "identity", vjust = -.5, hjust = .3, size = 3) +
  scale_fill_grey(
    breaks = c("ROR", "T"),
    labels = c("In-river", "Transported"),
    start = 0.8,
    end = 0.2
  ) +
  labs(
    x = "Smolt migration year", y = "Number of adult returns",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = c(.9, .8),
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)
  )
# dev.off()



# Figure S2
## Number of adult returns by passage type and smolt age
# pdf(here("results", "figures", "manuscript", "FigS2 Passtype and tag age & location.pdf"), width=8, height=4, onefile=TRUE)
DATA %>%
  group_by(SYage1, OceanAge, pass_type_T_R, rel_age) %>%
  summarise(n = n()) %>%
  ggplot(aes(x = as.factor(SYage1), y = n, fill = rel_age)) +
  geom_bar(position = "stack", stat = "identity", width = .6) +
  scale_fill_grey(
    breaks = c("ParrAboveLGR", "SmoltAboveLGR", "SmoltAtLGR"),
    labels = c("Parr/aboveLGR", "Smolt/aboveLGR", "Smolt/atLGR")
  ) +
  labs(
    x = "Smolt migration year", y = "Number of adult returns",
    color = NULL,
    fill = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)
  ) +
  facet_wrap(~pass_type_T_R, labeller = labeller(pass_type_T_R = p.labs))
# dev.off()



# Figure S3
## proportion of origin grouping: stacked barplot
# pdf(here("results", "figures", "manuscript", "FigS3 Origin.pdf"), width=7, height=4, onefile=TRUE)
DATA %>%
  group_by(SYage1, MPG) %>%
  summarise(n = n()) %>%
  mutate(MPG = factor(MPG, levels = c(
    "Grande Ronde/Imnaha", "Lower Snake River", # look into ordering data with fill order
    "Clearwater River", "South Fork Salmon River",
    "Lower Salmon River", "Middle Fork Salmon River",
    "Upper Salmon River"
  ))) %>%
  ggplot(aes(x = SYage1, y = n, fill = MPG, label = MPG)) +
  geom_bar(position = position_fill(), stat = "identity", width = .8) +
  labs(
    fill = "Origin",
    x = "Smolt migration year",
    y = "Proportion"
  ) +
  scale_x_continuous(breaks = seq(1998, 2020, by = 1)) +
  scale_fill_manual(
    breaks = c(
      "Grande Ronde/Imnaha", "Lower Snake River", # look into ordering data with fill order
      "Clearwater River", "South Fork Salmon River",
      "Lower Salmon River", "Middle Fork Salmon River",
      "Upper Salmon River"
    ),
    labels = c(
      "Grande Ronde/Imnaha", "L Snake",
      "Clearwater", "SF Salmon",
      "L Salmon", "MF Salmon",
      "U Salmon"
    ),
    values = pal_all
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
# dev.off()




# Figure S4
## Number by origin grouping and passage type
# pdf(here("results", "figures", "manuscript", "FigS4 Passtype and Origin.pdf"), width=8, height=5, onefile=TRUE)
DATA %>%
  group_by(SYage1, OceanAge, pass_type_T_R, MPG) %>%
  summarise(n = n()) %>%
  mutate(MPG = factor(MPG, levels = c(
    "Grande Ronde/Imnaha", "Lower Snake River",
    "Clearwater River", "South Fork Salmon River",
    "Lower Salmon River", "Middle Fork Salmon River",
    "Upper Salmon River"
  ))) %>%
  ggplot(aes(x = as.factor(SYage1), y = n, fill = MPG)) +
  geom_bar(position = "stack", stat = "identity", width = .6) +
  scale_fill_manual(
    breaks = c(
      "Grande Ronde/Imnaha", "Lower Snake River",
      "Clearwater River", "South Fork Salmon River",
      "Lower Salmon River", "Middle Fork Salmon River",
      "Upper Salmon River"
    ),
    labels = c(
      "Grande Ronde/Imnaha", "L Snake",
      "Clearwater", "SF Salmon",
      "L Salmon", "MF Salmon",
      "U Salmon"
    ),
    values = pal_all
  ) +
  labs(
    x = "Smolt migration year", y = "Number of adult returns",
    color = NULL,
    fill = NULL
  ) +
  theme_classic() +
  # scale_x_discrete(breaks = every_nth(n = 2))+
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1)
  ) +
  guides(fill = guide_legend(nrow = 2)) +
  facet_wrap(~pass_type_T_R, labeller = labeller(pass_type_T_R = p.labs))
# dev.off()



# Figure S5
## Fork length
# pdf(here("results", "figures", "manuscript", "FigS5 Fork Length.pdf"), width=4, height=8, onefile=TRUE)
FLz_xaxis <- (seq(40, 160, 10) - mean(DATA$length)) / sd(DATA$length)
FL_lab <- as.character(seq(40, 160, 10))

ggplot(data.scaled, aes(x = length, y = as.character(SYage1), fill = factor(after_stat(quantile)))) +
  stat_density_ridges(
    rel_min_height = 0.03, scale = 0.9,
    quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975),
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = "o", point_size = 0.7, point_alpha = 0.7, alpha = 0.5, size = .25
  ) +
  scale_fill_manual(
    name = "Probability", values = alpha(c("#333333", "#CCCCCC", "#CCCCCC", "#333333"), 0.8),
    labels = c("(0, 0.025]", "(0.025, 0.5]", "(0.5, 0.975]", "(0.975, 1]")
  ) +
  labs(x = "Fork Length (mm)", y = "Smolt Migration Year") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.line.y = element_line(linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = FLz_xaxis, labels = FL_lab)
# dev.off()



# Figure S6
## DOY of LGR passage
# pdf(here("results", "figures", "manuscript", "FigS6 LGR passage DOY.pdf"), width=4, height=8, onefile=TRUE)
DOYz_xaxis <- (seq(80, 160, 10) - mean(DATA$SYage1_DOY)) / sd(DATA$SYage1_DOY)
DOY_lab <- as.character(seq(80, 160, 10))

ggplot(data.scaled, aes(x = SYage1_DOY, y = as.character(SYage1), fill = factor(after_stat(quantile)))) +
  stat_density_ridges(
    rel_min_height = 0.03, scale = 0.9,
    quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975),
    geom = "density_ridges_gradient",
    calc_ecdf = TRUE,
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = "o", point_size = 0.7, point_alpha = 0.7, alpha = 0.5, size = .25
  ) +
  scale_fill_manual(
    name = "Probability", values = alpha(c("#333333", "#CCCCCC", "#CCCCCC", "#333333"), 0.8),
    labels = c("(0, 0.025]", "(0.025, 0.5]", "(0.5, 0.975]", "(0.975, 1]")
  ) +
  labs(x = "Day of Year", y = "Smolt Migration Year") +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.line.y = element_line(linewidth = 0.2)
  ) +
  scale_x_continuous(breaks = DOYz_xaxis, labels = DOY_lab)
# dev.off()



# Figure S7
## LGR River Temperature
# pdf(here("results", "figures", "manuscript", "FigS7 LGR temperature.pdf"), width=7, height=8, onefile=TRUE)
ggplot(DATA, aes(x = SYage1_DOY, y = LGR.temp.7d, group=as.factor(SYage1))) +
  xlim(80, 160) +
  ylim(2, 18) +
  geom_point(size = 0.25, color = pal_cov[2] ) +
  geom_line(colour = pal_cov[2]) +
  gghighlight::gghighlight( )+
  labs(x = "Day of year", y = expression(paste("River temperature (", degree, "C)"))) +
  facet_wrap(~SYage1) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.line.y = element_line(linewidth = 0.2)
  )
# dev.off()



# Figure S8
# LGR River flow
# pdf(here("results", "figures", "manuscript", "FigS8 LGR flow.pdf"), width=7, height=8, onefile=TRUE)
ggplot(DATA, aes(x = SYage1_DOY, y = LGR.flow.7d, group=as.factor(SYage1))) +
  xlim(80, 160) +
  geom_point(size = 0.25, color = pal_cov[2] ) +
  geom_line(colour = pal_cov[2]) +
  gghighlight::gghighlight( )+
  labs(x = "Day of year", y = "River flow (kcfs)") +
  facet_wrap(~SYage1) +
  theme_bw() +
  theme(
    legend.position = "none",
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.line.y = element_line(linewidth = 0.2)
  )
# dev.off()



# Figure S9.
## Marine indices

### NPGO raw data
NPGO.raw.data <- fread("https://www.o3d.org/npgo/data/NPGO.txt",
                       skip = 30,
                       header = FALSE,
                       fill=TRUE
)
NPGO.raw.data <- NPGO.raw.data[1:889,]

NPGO.ONDJFM.T <- NPGO.raw.data %>%
  rename(year = 1, month = 2, NPGO = 3) %>%
  mutate(CY = as.numeric(year)) %>%
  mutate(BYage0 = ifelse(month %in% c(1:9), CY - 1, CY)) %>%
  mutate(
    SYage1 = BYage0,
    AYageT = BYage0 + 1
  ) %>%
  #  mutate(across(!NPGO, as.character)) %>%
  filter(month %in% c(10:12, 1:3),
         between (SYage1, 1998, 2023)) %>%
  select(SYage1, month, NPGO)


# 3-year NPGO
NPGO.raw.data <- NPGO.raw.data %>%
  rename(year = 1, month = 2, NPGO = 3) %>%
  mutate(CY = as.numeric(year)) %>%
  mutate(
    BYage0 = ifelse(month %in% c(1:9), CY - 1, CY),
    WY = ifelse(month %in% c(1:9), CY, (CY + 1))
  ) %>%
  mutate(
    SYage1 = BYage0,
    AYageT = BYage0 +1
  ) %>%
  mutate(across(!NPGO, as.character))

NPGO.ONDJFM.T.3 <- NPGO.raw.data %>%
  filter(month %in% c(10:12, 1:3)) %>%
  summarise(.by = SYage1, NPGO.ONDJFM.T = round(mean(NPGO), 3))
NPGO.ONDJFM.T.3$NPGO.ONDJFM.T.AVG3 <- rollmean(NPGO.ONDJFM.T.3$NPGO.ONDJFM.T,
                                               k=3, align="left", fill=c(NA,NA,NA))
NPGO.ONDJFM.T.3 <- NPGO.ONDJFM.T.3[,-2]


# NPGO plots
p_raw_npgo <- ggplot(NPGO.ONDJFM.T, aes(x = SYage1, y = NPGO)) +
  geom_point(shape = rep(c("O", "N", "D", "J", "F", "M"), 26), size = 2) +
  stat_summary(
    geom = "line", group = 1, fun = "mean",
    size = 1, colour = pal_cov[3], alpha = 0.5
  ) +
  stat_summary(
    geom = "point", group = 1, fun = "mean",
    shape = 19, size = 3, col = pal_cov[3], alpha = 0.7
  ) +
  scale_x_continuous(breaks = seq(1998,2023, 1))+
  labs(x = "", y = "NPGO index") +
  ggtitle("(a)") +
  theme_minimal() +
  theme(
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.line.y = element_line(linewidth = 0.2)
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
    plot.title = element_text(size = 10)
  )

# 3-year NPGO
NPGO.ONDJFM.T.3.ggplot <- NPGO.ONDJFM.T.3[NPGO.ONDJFM.T.3$SYage1 %in% 1998:2023,]
NPGO.ONDJFM.T.3.ggplot[,1] <- as.numeric(NPGO.ONDJFM.T.3.ggplot[,1])

p_npgo3 <- ggplot(NPGO.ONDJFM.T.3.ggplot, aes(x = SYage1, y = NPGO.ONDJFM.T.AVG3)) +
  geom_line(size = 1, colour = pal_cov[3], alpha = 0.5) +
  geom_point(shape = 19, size = 3, col = pal_cov[3], alpha = 0.7) +
  scale_x_continuous(breaks = seq(1998,2023, 1))+
  labs(x = expression("Year" ~ italic(T) ),
       y = "3-year mean of NPGO index") +
  ggtitle("(b)") +
  theme_minimal() +
  theme(
    axis.line.x = element_line(linewidth = 0.2),
    axis.ticks = element_line(linewidth = 0.2),
    axis.line.y = element_line(linewidth = 0.2)
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 8, r = 0, b = 0, l = 0)),
    plot.title = element_text(size = 10)
  )

patch3 <- p_raw_npgo / p_npgo3

# pdf(here("results", "figures", "manuscript", "FigS9 NPGO.pdf"), width=7, height=7, onefile=TRUE)
patch3
# dev.off()

