library(rstan)
library(ggplot2)
library(rstanarm)
library(bayesplot)
library(tidyr) #per pivot_longer
library(dplyr) #per %>%
library(truncnorm) # per rnorm con valori minimi e massimi
library(loo)
library(lubridate) # progettato per semplificare il lavoro con le date

# Create custom palette because I wanted to distinguish well between countries
custom_palette <- c("#9F002E", "#BC72F0", "#FF50FF", "#FF7F00", "#FFB900", "#66E5FF", "#5FCEBE", "#00FF00", "#4DAF4A")

# Load the dataset
#install.packages("mlmRev")
library(mlmRev)
data("Mmmec")

# Check for missing values
colSums(is.na(Mmmec))

# Summary statistics of the deaths overall
summary(Mmmec$deaths)
# 25% of the observations has less than 8 deaths
# median(deaths) = 14.50
# mean(deaths) = 27.83
# 75% of the observations has less than 31 deaths
# max(deaths) = 313
# mean(deaths) = mean(expected deaths) = 27.8
summary(Mmmec$uvb)

# Explorative analysis

# Some histograms of the distribution of deaths and expected deaths

# Only deaths
ggplot(Mmmec, aes(x = deaths)) +
  geom_histogram(binwidth = 1, fill = "red", color = "black") +
  ggtitle("Distribution of Deaths")
ggsave(file="images/distribution_deaths.pdf", width=12,height=7)

# Only expected deaths
ggplot(Mmmec, aes(x = expected)) +
  geom_histogram(binwidth = 1, fill = "yellow", color = "black") +
  ggtitle("Distribution of Expected Deaths")
ggsave(file="images/distribution_expecteddeaths.pdf", width=12,height=7)

# Both deaths and expected deaths
# Gather the data into a long format using pivot_longer
Mmmec_long <- Mmmec %>%
  pivot_longer(cols = c(deaths, expected), names_to = "type", values_to = "value")
# Create an overlapped (position="identity") histogram with transparency (alpha=0.6)
ggplot(Mmmec_long, aes(x = value, fill = type)) +
  geom_histogram(binwidth = 1, alpha = 0.6, position = "identity", color = "black") +
  ggtitle("Distribution of Deaths and Expected Deaths") +
  xlab("Deaths / Expected Deaths") +
  scale_fill_manual(values = c("deaths" = "red", "expected" = "yellow"), labels = c("Deaths", "Expected Deaths"))
ggsave(file="images/distribution_deaths_expected.pdf", width=12,height=7)

# Some scatter plots of deaths VS UVB exposure

# With a single regression line
ggplot(Mmmec, aes(x = uvb, y = deaths, color = nation)) +
  geom_smooth(method = "lm", color = "blue") +
  scale_colour_manual(values = custom_palette) +
  ggtitle("Deaths VS UVB Exposure") +
  geom_jitter()
ggsave(file="images/scatter_deaths_vs_uvb.pdf", width=9.3,height=7)

# With a regression line for each nation
ggplot(Mmmec, aes(x = uvb, y = deaths, color = nation)) +
  geom_smooth(method = "lm", se = FALSE) +  # Adds a regression line for each nation
  scale_colour_manual(values = custom_palette) +
  ggtitle("Deaths VS UVB Exposure") +
  geom_jitter()
ggsave(file="images/scatter_deaths_vs_uvb_regression_by_nation.pdf", width=9.3,height=7)

# Or dividing by nation using facet_wrap
ggplot(Mmmec, aes(x = uvb, y = deaths, color = nation)) +
  geom_point() +
  geom_smooth(method = "lm", color = "blue") +
  facet_wrap(~ nation) +
  scale_colour_manual(values = custom_palette) +
  ggtitle("Deaths VS UVB Exposure by Nation") +
  theme(legend.position = "none")
ggsave(file="images/scatter_deaths_vs_uvb_by_nation.pdf", width=8,height=7)

# Boxplots of deaths by nation
ggplot(Mmmec, aes(x = nation, y = deaths)) +
  geom_boxplot(fill = custom_palette)
ggsave(file="images/boxplot_deaths_by_nation.pdf", width=8,height=7)


# Generalized linear mixed model with frequentist approach (in particular ML approach)
M1 <- glmer(deaths ~ uvb + (1 | region) + (1 | nation), Mmmec, poisson, offset = log(expected))
# (1|region) includes varying intercepts for each region
# log(expected) is an offset term to adjust the model to account for the expected number of deaths
summary(M1)


# Bayes approach relying on HMC sampling from the posterior distribution
M1.rstanarm <- stan_glmer(deaths ~ uvb + (1 | region) + (1 | nation), Mmmec, poisson, offset = log(expected))
print(M1.rstanarm)
summary(M1.rstanarm)

# Compute AIC and BIC for comparing with M1 results
loo_M1rs <- loo(M1.rstanarm) # extract LOO
print(loo_M1rs)
elpd_loo_M1rs <- loo_M1rs$estimates["elpd_loo", "Estimate"]
p <- length(M1.rstanarm$coefficients) # number of estimated parameters
n <- length(M1.rstanarm$y) # number of observations = 354
aic_M1rs <- - 2 * elpd_loo_M1rs + 2 * p
print(aic_M1rs) # AIC
bic_M1rs <- - 2 * elpd_loo_M1rs + log(n) * p
print(bic_M1rs) # BIC


# Posterior predictive check of reproduced data
pdf("images/densitiesM1rstanarm.pdf", width=8, height=7)
pp_check(M1.rstanarm)
dev.off()

# Plot the credible intervals
beta_names <- c(paste0("beta^", c("uvb")), "gl.intercept")
alpha_names<-c()
for (i in 1:51){
  alpha_names[i] <- paste0(expression(alpha), "[", i,"]")
}
posterior_M1 <- as.matrix(M1.rstanarm)
pdf("images/logistic_credible_intervals.pdf", width=10, height=11)
mcmc_intervals(posterior_M1, regex_pars=c( "uvb",
                                           "(Intercept)", "b"))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.4))+
  scale_y_discrete(labels = ((parse(text= c(beta_names, alpha_names)))))
dev.off()

# Plot the posterior marginal densities along with 50% intervals for the ‘fixed-effects’ uvb
mcmc_areas(posterior_M1, regex_pars = c("uvb"))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.4))
ggsave(file = "images/logistic_fixed_effects.pdf", width=8, height=7)

# Plot the random effects (posterior mean +- s.e.)
pdf("images/random_effects_log.pdf", height =7, width =11)
int_ord <- sort(coef(M1.rstanarm)$region[,1], index.return=TRUE)$x
ord <- sort(coef(M1.rstanarm)$region[,1], index.return=TRUE)$ix
region.abbr <- as.character(1:79)
region.abbr.ord <- region.abbr[ord]
se_ord <- M1.rstanarm$ses[ord]
par(xaxt="n", mfrow=c(1,1), mar = c(5,2,2,1))
plot(1:length(int_ord), int_ord, ylim=c(-1.4,1.4), pch=19, bg=2, xlab="Regions", 
      ylab="Intercepts",  cex.main=1.9, cex.lab=1.9)
for (h in 1:length(int_ord)){
  segments(h, int_ord[h]-se_ord[h], h, int_ord[h]+se_ord[h], col="red")
  is.wholenumber <-
    function(x, tol = .Machine$double.eps^0.5) 
      abs(x - round(x)) < tol
  if (is.wholenumber(h/2)){
    text(h, int_ord[h]+se_ord[h]+0.1, region.abbr.ord[h], cex=1.1)}else{
      text(h, int_ord[h]-se_ord[h]-0.1, region.abbr.ord[h], cex=1.1)
    }
}
dev.off()




## Poisson regression

# Create UVB index as in Table III of the paper
uvbi_params <- data.frame(
  nation = c(unique(Mmmec$nation)),
  mean_UVBI = c(12.70, 12.79, 9.96, 17.18, 10.91, 21.45, 10.54, 13.26, 11.40),
  sd_UVBI = c(0.29, 1.35, 0.38, 2.59, 1.50, 3.51, 0.60, 0.05, 0.47),
  min_UVBI = c(12.17, 10.45, 9.47, 12.92, 6.69, 16.83, 9.64, 13.22, 10.62),
  max_UVBI = c(13.10, 15.15, 10.49, 23.24, 13.46, 28.95, 11.70, 13.31, 11.94)
)
# Join UVBI data to the dataframe Mmmec
Mmmec2 <- left_join(Mmmec, uvbi_params, by = "nation")
# From mean_UVBI and sd_UVBI generate UVBI values for each county, respecting min_UVBI and max_UVBI
Mmmec2 <- Mmmec2 %>%
  mutate(UVBI = rtruncnorm(n(), a = min_UVBI, b = max_UVBI, mean = mean_UVBI, sd = sd_UVBI))

stan_data <- list(
  N = nrow(Mmmec2),
  deaths = Mmmec2$deaths,
  expected = Mmmec2$expected,
  K = length(unique(Mmmec2$nation)),
  J = length(unique(Mmmec2$region)),
  k = as.integer(as.factor(Mmmec2$nation)),
  j = as.integer(as.factor(Mmmec2$region)),
  UVBI = Mmmec2$UVBI
)
comp_model_A <- stan_model('poisson_regression.stan')
fit_model_A <- sampling(comp_model_A, data = stan_data, seed = 123)
print(fit_model_A, pars = c('beta0','beta1','sigma_s','sigma_u','sigma_e'))
y_rep <- as.matrix(fit_model_A, pars = "y_rep")


# Posterior predictive checking

pdf(file="images/traceplot.pdf", width=8, height=7)
traceplot(fit_model_A, pars = c('beta0', 'beta1', 'sigma_s', 'sigma_u', 'sigma_e'))
dev.off()

# densities 
pdf(file="images/densities.pdf", width=8, height=7)
ppc_dens_overlay(y = stan_data$deaths, y_rep[1:200,])+
  xaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
dev.off()

# empirical distribution function
pdf(file="images/ecdf.pdf", width=8, height=7)
ppc_ecdf_overlay(y = stan_data$deaths, y_rep[1:200,])+
  xaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
dev.off()

# proportion of zero
pdf(file="images/proportion_zero.pdf", width =8, height =7)
prop_zero <- function(x) mean(x == 0)
ppc_stat(y = stan_data$deaths, yrep = y_rep, stat = "prop_zero")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
dev.off()

# statistics
pdf(file="images/mean.pdf", width=5, height=5)
ppc_stat(y = stan_data$deaths, yrep = y_rep, stat="mean")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
dev.off()
pdf(file="images/sd.pdf", width=5, height=5)
ppc_stat(y = stan_data$deaths, yrep = y_rep, stat="sd")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
dev.off()
pdf(file="images/median.pdf", width=5, height=5)
ppc_stat(y = stan_data$deaths, yrep = y_rep, stat="median")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
dev.off()
pdf(file="images/max.pdf", width=5, height=5)
ppc_stat(y = stan_data$deaths, yrep = y_rep, stat="max")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
dev.off()

# standardized residuals
pdf(file="images/standardized_residuals.pdf", width =8, height =7)
mean_y_rep <- colMeans(y_rep)
std_resid <- (stan_data$deaths - mean_y_rep) / sqrt(mean_y_rep)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)+
  labs(x="Mean of y_rep", y= "Standardized residuals")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16))
dev.off()

# predictive intervals
pdf(file="images/predictive_intervals.pdf", width =8, height =7)
ppc_intervals(
  y = stan_data$deaths, 
  yrep = y_rep,
  x = stan_data$uvb
) + 
  labs(x = "Deaths", y = "Expected deaths")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))
dev.off()

# Compute AIC, BIC and log-likelihood for comparing with previous results
log_lik_A <- extract_log_lik(fit_model_A)
loo_A <- loo(log_lik_A)
print(loo_A)
elpd_loo_A <- loo_A$estimates["elpd_loo", "Estimate"]
p <- 5 + length(unique(Mmmec2$nation)) + length(unique(Mmmec2$region)) + nrow(Mmmec2) # number of estimated parameters
aic_A <- - 2 * elpd_loo_A + 2 * p
print(aic_A)
bic_A <- - 2 * elpd_loo_A + log(n) * p
print(bic_A)
