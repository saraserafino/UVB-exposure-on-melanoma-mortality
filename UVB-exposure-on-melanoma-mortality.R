library(rstan)
library(ggplot2)
library(rstanarm)
library(bayesplot)
library(tidyr) # for pivot_longer
library(dplyr) # for %>%
library(truncnorm) # for rnorm with minimum and maximum values
library(loo)
library(lubridate) # for simplifying working with dates
library(posterior) # for checking warnings of stan model

# Create custom palette because I want to distinguish well between nations and don't like the existing options
custom_palette <- c("#9F002E", "#B23AEE", "#FF50FF", "#FF7F00", "#FFB900", "#00EEEE", "#4EEE94","#458B00", "#4876FF")

# Load the dataset
#install.packages("mlmRev")
library(mlmRev)
data("Mmmec")

# Check for missing values
colSums(is.na(Mmmec))

# Summarize statistics about deaths and uvb overall
summary(Mmmec$deaths)
# 1st quantile = 8 means that 25% of the observations have less than 8 deaths
# 3rd quantile = 31 means that 75% of the observations have less than 31 deaths
deaths_by_country <- aggregate(deaths ~ nation, data = Mmmec, sum)
print(deaths_by_country)
summary(Mmmec$uvb)

# Explorative analysis

# Some histograms of the distribution of deaths and expected deaths

# Only deaths
ggplot(Mmmec, aes(x = deaths)) +
  geom_histogram(binwidth = 1, fill = "red", color = "black") +
  ggtitle("Distribution of Deaths")
ggsave(file="images/distribution_deaths.pdf", width=18,height=10.5)

# Only expected deaths
ggplot(Mmmec, aes(x = expected)) +
  geom_histogram(binwidth = 1, fill = "yellow", color = "black") +
  ggtitle("Distribution of Expected Deaths")
ggsave(file="images/distribution_expecteddeaths.pdf", width=18,height=10.5)

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
ggsave(file="images/distribution_deaths_expected.pdf", width=18,height=10.5)

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

# Or dividing by nation using facet_wrap, so that it is less chaotic
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
  geom_boxplot(fill = custom_palette) +
  ggtitle("Deaths across counties in nations")
ggsave(file="images/boxplot_deaths_by_nation.pdf", width=8,height=7)


# Generalized linear mixed model with frequentist approach (in particular hierarchical approach)
M1 <- glmer(deaths ~ uvb + (1 | region) + (1 | nation), Mmmec, poisson, offset = log(expected))
# (1 | k) includes varying intercepts for each k
# log(expected) is an offset term to adjust the model to account for the expected number of deaths
summary(M1)

# Bayes approach relying on HMC sampling from the posterior distribution
M1.rstanarm <- stan_glmer(deaths ~ uvb + (1 | region) + (1 | nation), Mmmec, poisson, offset = log(expected))
print(M1.rstanarm)
summary(M1.rstanarm)

# Extract Leave-One-Out Cross-Validation
loo_M1rs <- loo(M1.rstanarm)
print(loo_M1rs)
# Using 10-fold CV as suggested by warning due to some observations with too high pareto_k,
# there is a warning about 1 divergent transition after warmup
#kfold_result <- kfold(M1.rstanarm, K = 10)
# But since Rhat=1 and ESS (n_eff)>1000, the resulting posterior is often good enough to move forward


# Posterior predictive checks
y_rep <- posterior_predict(M1.rstanarm)

# densities 
pp_check(M1.rstanarm)
ggsave(file="images/densities_M1rstanarm.pdf", width=8,height=7)

# Plot the credible intervals
beta_names <- c(paste0("beta^", c("uvb")), "gl.intercept")
alpha_names<-c()
for (i in 1:25){
  alpha_names[i] <- paste0("region[", i,"]")
} # region 26 is missing
for (i in 27:79){
  alpha_names[i-1] <- paste0("region[", i,"]")
}
alpha_names[79] <- paste0("Belgium")
alpha_names[80] <- paste0("WG")
alpha_names[81] <- paste0("Denmark")
alpha_names[82] <- paste0("France")
alpha_names[83] <- paste0("UK")
alpha_names[84] <- paste0("Italy")
alpha_names[85] <- paste0("Ireland")
alpha_names[86] <- paste0("Luxembourg")
alpha_names[87] <- paste0("Netherlands")
alpha_names[88] <- paste0(expression(sigma[region]))
alpha_names[89] <- paste0(expression(sigma[nation]))
posterior_M1 <- as.matrix(M1.rstanarm)
mcmc_intervals(posterior_M1, regex_pars=c( "uvb",
                                           "(Intercept)", "b"))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.4))+
  scale_y_discrete(labels = ((parse(text= c(beta_names, alpha_names)))))
ggsave(file="images/logistic_credible_intervals.pdf", width=5, height=length(alpha_names) * 0.25)

# Plot the posterior marginal densities along with 50% intervals for the ‘fixed-effects’ uvb
mcmc_areas(posterior_M1, regex_pars = c("uvb"))+
  xaxis_text(on =TRUE, size=rel(1.9))+
  yaxis_text(on =TRUE, size=rel(1.4))
ggsave(file = "images/logistic_fixed_effects.pdf", width=8, height=7)

# Plot the random effects (posterior mean +- s.e.)
int_ord <- sort(coef(M1.rstanarm)$region[,1], index.return=TRUE)$x
ord <- sort(coef(M1.rstanarm)$region[,1], index.return=TRUE)$ix
region.abbr <- levels(Mmmec$region)
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
ggsave(file="images/random_effects_log.pdf", height=7, width=length(int_ord) * 0.4)

# empirical distribution function
ppc_ecdf_overlay(y = M1.rstanarm$y, y_rep[1:200,])+
  xaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
ggsave(file="images/ecdf_M1.rstanarm.pdf", width=8, height=7)

# proportion of zero
prop_zero <- function(x) mean(x == 0)
ppc_stat(y = M1.rstanarm$y, yrep = y_rep, stat = "prop_zero")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
ggsave(file="images/proportion_zero_M1.rstanarm.pdf", width=8, height=7)

# statistics
ppc_stat(y = M1.rstanarm$y, yrep = y_rep, stat="mean")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
ggsave(file="images/mean_M1.rstanarm.pdf", width=5, height=5)

ppc_stat(y = M1.rstanarm$y, yrep = y_rep, stat="sd")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
ggsave(file="images/sd_M1.rstanarm.pdf", width=5, height=5)

ppc_stat(y = M1.rstanarm$y, yrep = y_rep, stat="median")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
ggsave(file="images/median_M1.rstanarm.pdf", width=5, height=5)

ppc_stat(y = M1.rstanarm$y, yrep = y_rep, stat="max")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
ggsave(file="images/max_M1.rstanarm.pdf", width=5, height=5)

# standardized residuals
mean_y_rep <- colMeans(y_rep)
std_resid <- (M1.rstanarm$y - mean_y_rep) / sqrt(mean_y_rep)
qplot(mean_y_rep, std_resid) + hline_at(2) + hline_at(-2)+
  labs(x="Mean of y_rep", y= "Standardized residuals")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16))
ggsave(file="images/standardized_residuals_M1.rstanarm.pdf", width =8, height =7)

# predictive intervals
ppc_intervals(
  y = M1.rstanarm$y, 
  yrep = y_rep) + 
  labs(x = "Deaths", y = "Expected deaths")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))
ggsave(file="images/predictive_intervals_M1.rstanarm.pdf", width =8, height =7)



# Stan model of model A of the paper
# It is a variance components model with UVBI included in the fixed part of the model and
# random terms s_k, u_jk, e_ijk associated respectively with the intercept at level 3, 2, 1
# s_k ~ normal(0, sigma_s)
# u_jk ~ normal(0, sigma_u)

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
# Generate UVBI values for each county from mean_UVBI and sd_UVBI, respecting min_UVBI and max_UVBI
Mmmec2 <- Mmmec2 %>%
  mutate(UVBI = rtruncnorm(n(), a = min_UVBI, b = max_UVBI, mean = mean_UVBI, sd = sd_UVBI))
# Prepare data for stan model
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
fit_model_A <- sampling(comp_model_A, data = stan_data, seed = 123, iter=4000)
print(fit_model_A, pars = c('beta0','beta1','sigma_s','sigma_u','sigma_e'))
y_rep_model_A <- as.matrix(fit_model_A, pars = "y_rep")

# Checks due to warnings

check_divergences(fit_model_A) # no divergences
check_treedepth(fit_model_A) # no saturations of max tree depths of 10
check_energy(fit_model_A)
# Bayesian Fraction of Missing Information is in fact low. This implies that the adaptation phase of the Markov
# chains did not turn out well or the posterior has thick tails that were not well explored in the simulation.
# Should reparametrize the model but then it would be different from the paper's

# The warning about bulk and tail ESS too low disappears with 8000 iterations instead of 4000
draws <- as_draws_array(fit_model_A)
ess <- summarise_draws(draws)
print(ess)
# For iter=8000: ess_bulk=547, ess_tail=927, Rhat=1.00
# But since posterior predictive checks do not improve significantly, it is unnecessary


# Posterior predictive checks

# Traceplot of the 4 chains to see if they mix well
traceplot(fit_model_A, pars = c('beta0', 'beta1', 'sigma_s', 'sigma_u', 'sigma_e'))
ggsave(file="images/traceplot.pdf", width=8, height=7)

# densities 
ppc_dens_overlay(y = stan_data$deaths, y_rep_model_A[1:200,])+
  xaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
ggsave(file="images/densities_model_A.pdf", width=8, height=7)

# empirical distribution function
ppc_ecdf_overlay(y = stan_data$deaths, y_rep_model_A[1:200,])+
  xaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
ggsave(file="images/ecdf_model_A.pdf", width=8, height=7)

# proportion of zero
prop_zero <- function(x) mean(x == 0)
ppc_stat(y = stan_data$deaths, yrep = y_rep_model_A, stat = "prop_zero")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  legend_text(size=rel(4))
ggsave(file="images/proportion_zero_model_A.pdf", width=8, height=7)

# statistics
ppc_stat(y = stan_data$deaths, yrep = y_rep_model_A, stat="mean")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
ggsave(file="images/mean_model_A.pdf", width=5, height=5)

ppc_stat(y = stan_data$deaths, yrep = y_rep_model_A, stat="sd")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
ggsave(file="images/sd_model_A.pdf", width=5, height=5)

ppc_stat(y = stan_data$deaths, yrep = y_rep_model_A, stat="median")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
ggsave(file="images/median_model_A.pdf", width=5, height=5)

ppc_stat(y = stan_data$deaths, yrep = y_rep_model_A, stat="max")+
  xaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text( size=22))+
  legend_text(size=rel(1.6))
ggsave(file="images/max_model_A.pdf", width=5, height=5)

# standardized residuals
mean_y_rep_model_A <- colMeans(y_rep_model_A)
std_resid <- (stan_data$deaths - mean_y_rep_model_A) / sqrt(mean_y_rep_model_A)
qplot(mean_y_rep_model_A, std_resid) + hline_at(2) + hline_at(-2)+
  labs(x="Mean of y_rep", y= "Standardized residuals")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text(size=16),
        axis.title.y = element_text(size=16))
ggsave(file="images/standardized_residuals_model_A.pdf", width =8, height =7)

# predictive intervals
ppc_intervals(
  y = stan_data$deaths, 
  yrep = y_rep_model_A,
  x = stan_data$uvb
) + 
  labs(x = "Deaths", y = "Expected deaths")+
  xaxis_text(on =TRUE, size=22)+
  yaxis_text(on =TRUE, size=22)+
  theme(axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20))
ggsave(file="images/predictive_intervals_model_A.pdf", width =8, height =7)

# Extract Leave-One-Out Cross-Validation
# Note that since I wrote my own stan model, I had to store the pointwise log-likelihood
log_lik_A <- extract_log_lik(fit_model_A)
loo_A <- loo(log_lik_A)
print(loo_A)