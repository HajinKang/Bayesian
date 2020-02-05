################################################################################
#### Code to accompany Baldwin and Larson -- An Introduction to Using Bayesian
#### Methods with Clinical Data.
#### This code and data will reproduce the analyses and model-based plots
#### in the Baldwin and Larson paper. We assume familiarity with R. The code
#### for the plots can get involved. Those wanting to just on the brms code
#### for fitting models can focus just on the beginning code.
#### Note that MCMC is a simulation based method. If you use the same random
#### number seed and the same version of the software as we report in the
#### paper, the results should be the same. If you use a different version
#### of Stan, brms, or R, your results may not match the results exactly.
#### However, the results should be quite close.
################################################################################

#### Load required libraries
#### Only brms is needed to fit the models. We used the other libraries for
#### creating plots and managing the data.
library(brms)
library(ggplot2)
library(dplyr)
library(reshape2)

### read data
### You will need to adjust file path depending on where you have saved the data
ern.data <- read.csv(file = "ern.data.csv", stringsAsFactors = FALSE)

### Model 1
mod1 <- brm(ern_mean ~ anxiety, data = ern.data,
            prior = c(set_prior("normal(0,3)", class = "b"),
                      set_prior("cauchy(0,2.5)", class = "sigma")),
            family = gaussian(), seed = 394589)
summary(mod1, waic = TRUE)

### Traceplot (Figure 3)
stanplot(mod1, type = "trace", inc_warmup = TRUE, ncol = 1) + scale_color_grey()

### Fitting Model 2
### The update() function allows you to add additional predictors,
### random effects, etc. to a model without needing to recompile
### the model. This can speed up model fitting significantly.
mod2 <- update(mod1, formula. = ~ . + sex, newdata = ern.data, seed = 39585)
summary(mod2, waic = TRUE)

### Fitting Model 3
mod3 <- update(mod2, formula. = ~ . + sex*anxiety, newdata = ern.data, seed = 854367)
summary(mod3, waic = TRUE)

### Comparing Models using The WAIC() and LOO() functions from brms
WAIC(mod1, mod2, mod3)
LOO(mod1, mod2, mod3)

################################################################################
### Code beyond this point is for creating the plots in the paper
################################################################################

### This code produces Figure 4 from the paper
### Instructions for how this works can be found at
### http://stackoverflow.com/questions/31215748/how-to-shade-part-of-a-density-curve-in-ggplot-with-no-y-axis-data
### Extract posterior samples from the 'mod1' brms object
mod1.samp <- posterior_samples(mod1)

### compute the 2.5 and 97.5% quantiles for the anxiety slope
mod1.ll <- quantile(mod1.samp$b_anxiety, .025)
mod1.ul <- quantile(mod1.samp$b_anxiety, .975)

p <- ggplot(data = mod1.samp, aes(x = b_anxiety)) +
  geom_density()

d <- ggplot_build(p)$data[[1]]
p <- p + geom_area(data = subset(d, x > mod1.ll & x < mod1.ul),
                   aes(x = x, y = y), fill = "gray90", color = "black") +
  geom_vline(xintercept = mean(mod1.samp$b_anxiety), size = 2) +
  xlim(c(-1, 1)) +
  ylab("Density \n") +
  xlab("\n Anxiety Slope Value") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))
p

### This code produces Figure 5 from the paper
### This is easiest to do if you create a new object
### with just the data for the lines (mod1.samp.gg)
mod1.samp.gg <- mod1.samp[1:25, 1:2]
mod1.samp.gg$iternum <- seq(1, length(mod1.samp.gg$b_Intercept))

p <- ggplot(data = ern.data, aes(x = anxiety, y = ern_mean)) +
  geom_point(size = 2) +
  geom_abline(data = mod1.samp.gg,
              aes(intercept = b_Intercept,
                  slope = b_anxiety,
                  group = iternum), alpha = 0.5) +
  xlab("\n Standardized Trait Anxiety") +
  ylab("Error Related Negativity \n") +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  theme_bw()
p

### This code produces Figure 6 from the paper
### The fitted() function from brms produces the expected values
### The predict() function from brms produces the predicted values (i.e., new data points)
pp <- fitted(mod1, summary = TRUE)
colnames(pp) <- c("Estimate", "Est.Error", "low", "high")
pp2 <- predict(mod1, summary = TRUE)
colnames(pp2) <- c("Estimate", "Est.Error", "low", "high")

dat <- as.data.frame(cbind(ern = standata(mod1)$Y, anxiety = standata(mod1)$X, pp))
dat2 <- as.data.frame(cbind(ern = standata(mod1)$Y, anxiety = standata(mod1)$X, pp2))

p <- ggplot(dat2, aes(x = anxiety, y = ern)) +
  geom_point() +
  geom_ribbon(data = dat2, aes(x = anxiety, ymin = low, ymax = high), fill = "gray80", alpha = 0.5) +
  geom_ribbon(data = dat, aes(x = anxiety, ymin = low, ymax = high), fill = "gray60", alpha = 0.5) +
  geom_ribbon(data = dat, aes(x = anxiety, ymin = Estimate - 0.01, ymax = Estimate + 0.01), fill = "black") +
  xlab("\n Standardized Trait Anxiety") +
  ylab("Error Related Negativity \n") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14))
p

### The following code creates Figure 7 from the manuscript
### It is fairly involved and requires a lot of intermediate objects
### to ensure that shading is appropriate and correct. We welcome feedback
### on simpler ways to achieve the same results with ggplot or other plotting
### functions in R.
mod3.samp <- posterior_samples(mod3)

#### Construct the gender-specific intercepts and slopes
mod3.samp$men.slope <- mod3.samp$b_anxiety
mod3.samp$women.slope <- mod3.samp$b_anxiety + mod3.samp$"b_anxiety:sex"
mod3.samp$men.int <- mod3.samp$b_Intercept
mod3.samp$women.int <- mod3.samp$b_Intercept + mod3.samp$"b_sex"

### construct the quantitles for each of the intercepts and slopes
dens.mod3 <- mod3.samp[,c("men.int", "women.int", "men.slope", "women.slope")]
mod3.int.men.ll <- quantile(dens.mod3$men.int, .025)
mod3.int.men.ul <- quantile(dens.mod3$men.int, .975)
mod3.int.women.ll <- quantile(dens.mod3$women.int, .025)
mod3.int.women.ul <- quantile(dens.mod3$women.int, .975)
mod3.slope.men.ll <- quantile(dens.mod3$men.slope, .025)
mod3.slope.men.ul <- quantile(dens.mod3$men.slope, .975)
mod3.slope.women.ll <- quantile(dens.mod3$women.slope, .025)
mod3.slope.women.ul <- quantile(dens.mod3$women.slope, .975)
colnames(dens.mod3) <- c("Intercept: Men", "Intercept: Women", "Slope: Men", "Slope: Women")

### Make a long dataset
dens.mod3.long <- melt(dens.mod3, id.vars = NULL)

### Create a new variables for faceting in ggplot
dens.mod3.long$int <- ifelse(dens.mod3.long$variable == "Intercept: Men", "Intercept",
                             ifelse(dens.mod3.long$variable == "Intercept: Women", "Intercept", "Slope"))

### Build the foundational part of the plot
p <- ggplot() +
  geom_density(data = dens.mod3.long, aes(x = value, linetype = variable)) +
  facet_grid( ~ int) +
  scale_linetype_manual(name = "Parameter", values = c("solid", "dotted", "solid", "dotted"))

### build datasets for shading
d <- ggplot_build(p)$data[[1]]
d.int.men <- filter(d, PANEL == 1 & group == 1, x > mod3.int.men.ll & x < mod3.int.men.ul)
d.int.men$int = "Intercept"
d.int.women <- filter(d, PANEL == 1 & group == 2, x > mod3.int.women.ll & x < mod3.int.women.ul)
d.int.women$int = "Intercept"
d.slope.men <- filter(d, PANEL == 2 & group == 3, x > mod3.slope.men.ll & x < mod3.slope.men.ul)
d.slope.men$int = "Slope"
d.slope.women <- filter(d, PANEL == 2 & group == 4, x > mod3.slope.women.ll & x < mod3.slope.women.ul)
d.slope.women$int = "Slope"

### add shading to the plot
p <- p +
  geom_area(data = d.int.women,
            aes(x = x, y = y), fill = "gray60", color = "black", alpha = .5, linetype = "dotted") +
  geom_area(data = d.int.men,
            aes(x = x, y = y), fill = "gray80", color = "black", alpha = .5, linetype = "solid") +
  geom_area(data = d.slope.women,
            aes(x = x, y = y), fill = "gray60", color = "black", alpha = .5, linetype = "dotted") +
  geom_area(data = d.slope.men,
            aes(x = x, y = y), fill = "gray80", color = "black", alpha = .5, linetype = "solid") +
  ylab("Density \n") +
  xlab("\n Sampled Parameter") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.text.x = element_text(size = 12))
p


### The following code creates Figure 8 in the manuscript
pred.mod3 <- fitted(mod3, summary = TRUE)
colnames(pred.mod3) <- c("Estimate", "Est.Error", "low", "high")

mod3.dat <- as.data.frame(cbind(ern = ern.data$ern_mean,
                                anxiety = ern.data$anxiety,
                                sex = ern.data$sex,
                                pred.mod3))
mod3.dat$gender <- ifelse(mod3.dat$sex == 1, "Women", "Men")


mod3.plot <- ggplot(data = mod3.dat, aes(x = anxiety, y = ern))
mod3.plot <- mod3.plot + geom_point() +
  geom_ribbon(aes(x = anxiety, ymin = low, ymax = high), fill = "gray80", alpha = 0.5) +
  geom_ribbon(aes(x = anxiety, ymin = Estimate - 0.01, ymax = Estimate + 0.01),
              fill = "black") +
  facet_wrap( ~ gender) +
  xlab("\n Standardized Trait Anxiety") +
  ylab("Error Related Negativity \n") +
  theme_bw() +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        strip.text.x = element_text(size = 12))
mod3.plot


###################################################################
##### Additional data analysis examples
###################################################################

### Logistic regression

rank.df <- read.csv("https://stats.idre.ucla.edu/stat/data/binary.csv")
head(rank.df)

logistic <- brm(admit ~ gre + gpa + factor(rank), data = rank.df,
                prior = c(set_prior("normal(0,3)", class = "b")),
                family = bernoulli(), seed = 295508)
summary(logistic, waic = TRUE)


### Poisson regression

awards <- read.csv("https://stats.idre.ucla.edu/stat/data/poisson_sim.csv")
head(awards)

poisson <- brm(num_awards ~ factor(prog) + math, data = awards,
               prior = c(set_prior("normal(0,3)", class = "b")),
               family = poisson(), seed = 3099221)
summary(poisson, waic = TRUE)


### Multilevel model

library(foreign)

math <- read.dta(file = "https://stats.idre.ucla.edu/stat/stata/examples/mlm_imm/imm10.dta")

## Upped the iterations because the Effective Sample Size was fairly
## small with the default number of iterations.
mlm <- brm(math ~ homework + (1 | schnum), data = math,
           prior = c(set_prior("normal(0,30)", class = "b", coef = "Intercept"),
                     set_prior("normal(0,3)", class = "b", coef = "homework"),
                     set_prior("cauchy(0,10)", class = "sd"),
                     set_prior("cauchy(0,10)", class = "sigma")),
           family = gaussian(), seed = 5092313,
           iter = 4000)
summary(mlm, waic = TRUE)
