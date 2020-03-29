
library(brms) # for the analysis
library(haven) # to load the SPSS .sav file
library(tidyverse) # needed for data manipulation.
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(ggmcmc)
library(ggthemes)
library(ggridges)

popular2data <- read_sav(file = "https://github.com/MultiLevelAnalysis/Datasets-third-edition-Multilevel-book/blob/master/chapter%202/popularity/SPSS/popular2.sav?raw=true")

model1 <- brm(popular ~ 1 + extrav,
							data = popular2data,
							warmup = 1000, iter = 3000,
							cores = 2, chains = 2,
							seed = 123) #to run the model

model_cens <- brm(popular | cens(censored) ~ 1 + extrav,
							data = popular2data %>% mutate(censored = sample(c(1,0), replace = T, size = n())),
							warmup = 1000, iter = 3000,
							cores = 2, chains = 2,
							seed = 123) #to run the model


bf(y | cens(censor_variable) ~ predictors)


for (n in 1:N) {
	// special treatment of censored data
	if (cens[n] == 0) {
		target += normal_lpdf(Y[n] | mu[n], sigma);
	} else if (cens[n] == 1) {
		target += normal_lccdf(Y[n] | mu[n], sigma);
	}
}


# Dirichlet

frml <- bf(proportions ~ .)

fit <- brm(frml,
					 data = R80_data,
					 family = dirichlet(link = ‘logit’, link_phi = ‘log’),
					 future = TRUE)


A = rbind(
	c(0.2, 0.3, 0.5),
	c(0.8, 0.1, 0.1)
)
df = data.frame(x = rnorm(2))
df$A = A
m = brm(A ~ 1 + x, data = df, family = dirichlet())


m_censored = brm(A | cens(censor_variable) ~ 1 + x, data = df %>% mutate(censored = sample(c(1,0), replace = T, size = n())), family = dirichlet())

# Beta

A = rbind(
	c(0.2),
	c(0.8)
)
df = data.frame(x = rnorm(2))
df$A = A

m = brm(A ~ 1 + x, data = df, family = "beta")


m_censored = brm(A | cens(censored) ~ 1 + x, data = df %>% mutate(censored = sample(c(1,0), replace = T, size = n())), family = "beta")

