
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





# Dirichlet

frml <- bf(proportions ~ .)

fit <- brm(frml,
					 data = R80_data,
					 family = dirichlet(link = 'logit', link_phi = 'log'),
					 future = TRUE)


A = rbind(
	c(0.2, 0.3, 0.5),
	c(0.8, 0.1, 0.1)
)
df = data.frame(x = rnorm(2))
df$A = A
m = brm(A ~ 1 + x, data = df, family = dirichlet())

DF = 
	res_4$proportions %>%
	select(-.variable) %>% 
	filter(alive==0) %>% 
	select(`Cell type category`, sample, PFI.time.2, .value) %>% 
	mutate(time = scale(log(PFI.time.2))) %>% 
	select(-PFI.time.2) %>% 
	group_by(sample) %>%
	mutate(.value = .value/sum(.value)) %>%
	ungroup() %>%
	spread(`Cell type category`, .value) %>%
	select(x = time, A.1 = endothelial , A.2 = epithelial, A.3 = fibroblast, A.4 = immune_cell)

my_df = DF[,1, drop=F]
my_df$A = DF[,2:5, drop=F] %>% as.matrix() 

m = brm(A ~ x, data = my_df, family = dirichlet())


m_censored = brm(A | cens(censor_variable) ~ 1 + x, data = df %>% mutate(censored = sample(c(1,0), replace = T, size = n())), family = dirichlet())

# Beta

A = rbind(
	c(0.2),
	c(0.8)
)
df = data.frame(x = rnorm(2))
df$A = A

m = brm(A ~ 1 + x, data = df, family = "beta")

data = df %>% mutate(censored = sample(c(1,0), replace = T, size = n()))

DF = 
	res_4$proportions %>%
	select(-.variable) %>% 
	select(`Cell type category`, sample, PFI.time.2, .value, alive) %>% 
	mutate(time = scale(log(PFI.time.2))) %>% 
	select(-PFI.time.2) %>% 
	group_by(sample) %>%
	mutate(.value = .value/sum(.value)) %>%
	ungroup() %>%
	mutate(.value = boot::logit(.value)) %>%
	spread(`Cell type category`, .value) %>%
	select(time = time, alive,  endothelial ,  epithelial,  fibroblast,  immune_cell)

m_censored = brm(x | cens(censored) ~ A , data = data, family = "normal", cores = 4)

