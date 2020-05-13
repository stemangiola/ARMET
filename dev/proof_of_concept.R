library(tidyverse)
library(rstan)
library(tidybayes)
library(brms)

my_theme = 	
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size=12),
		aspect.ratio=1,
		#axis.text.x = element_text(angle = 90, hjust = 1),
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
		axis.title.y  = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
	)

bind <- function(...) cbind(...)

logsumexp <- function (x) {
	y = max(x)
	y + log(sum(exp(x - y)))
}

softmax <- function (x) {
	exp(x - logsumexp(x))
}

get_alpha = function(C, slope_value){
	intercept = 1:C
	#intercept[(C-1):C] = intercept[(C-1):C] * 2
	#intercept = intercept[(C-1):C] * 2
	intercept = intercept /5
	slope = rep(0, C)
	slope[4] = slope_value
	
	alpha = matrix(intercept %>%	c(slope), ncol = 2)
}

simulate_infiltration_process = function(X, .alpha){
	
	X %*% t(.alpha) %>%
		as.data.frame() %>%
		as_tibble() %>%
		mutate(sample = 1:n()) %>%
		mutate(risk = X[,2]) %>%
		gather(cell_type, rate, -sample, -risk) 
	
	
}

C = 5
S = 100
X = matrix(rep(1, S) %>% 	c(seq(0, 1, len = S)) , ncol = 2)
alpha = get_alpha(C, 4)

set.seed(123)

# my_prop_dir_2 = 
# 	simulate_infiltration_process(X, alpha) %>%
# 	group_by(risk) %>%
# 	arrange(sample) %>%
# 	mutate(proportion = gtools::rdirichlet(1, alpha = softmax(log(rate)) * 40)) %>%
# 	ungroup()
# 
# my_prop_dir_2 %>%
# 	ggplot(aes(risk, proportion, color=cell_type)) + 
# 	geom_point()
# 
# my_prop_dir_2 %>%
# 	mutate(count = rnbinom(n(), mu = rate * 100, size = 20)) %>%
# 	ggplot(aes(risk, count, color=cell_type)) + 
# 	geom_point()

my_prop_dir_3 =
	simulate_infiltration_process(X, alpha) %>%
	mutate(count = rnbinom(n(), mu = rate * 100, size = 100)) %>%
	group_by(risk) %>%
	arrange(sample) %>%
	mutate(proportion = count / sum(count)) %>%
	ungroup()




m = rstan::stan_model("inst/stan/proof_concept.stan")
fit = 
	sampling(
		m,
		data = list(
			S = S,
			A = 2,
			C = C,
			X = X,
			prop = 
				my_prop_dir_3 %>%
				select(sample, cell_type, proportion) %>%
				spread(cell_type, proportion) %>%
				nanny::as_matrix(rownames="sample")
		),
		cores = 4
	)


# Plotting

list(
	my_prop_dir_3 %>%
		ggplot(aes(risk, rate, color=cell_type)) + 
		geom_point() +
		scale_color_brewer(palette="Set1") +
		xlab("Risk") + 
		ylab("Cell type abundance rate") +
		my_theme, 

	my_prop_dir_3 %>%
		ggplot(aes(risk, count, color=cell_type)) + 
		geom_point() +
		scale_color_brewer(palette="Set1") +
		xlab("Risk") + 
		ylab("Cell type count") +
		my_theme,
	
	my_prop_dir_3 %>%
		ggplot(aes(risk, proportion, color=cell_type)) + 
		geom_point() +
		scale_color_brewer(palette="Set1") +
		xlab("Risk") + 
		ylab("Cell type proportion") +
		my_theme,
	
	fit %>%
		gather_draws(alpha_generative[A, C]) %>% 
		filter(A ==2) %>%
		ungroup() %>%
		mutate(C = C +1) %>%
		bind_rows({
			my_sd = (.) %>% group_by(C) %>% summarise(sd(.value)) %>% pull(2) %>% mean
			
			(.) %>%
				filter(C ==2) %>%
				mutate(C = rep(1, n())) %>%
				mutate(.value = rnorm(n(), 0, my_sd))
			
		})	 %>% 
		ggplot(aes(.value, color=factor(C))) + 
		geom_density() +
		scale_color_brewer(palette="Set1") +
		xlab("Slope") + 
		ylab("Probability density") +
		my_theme,
	
	fit %>%
		gather_draws(alpha_descriptive[A, C]) %>% 
		filter(A ==2) %>%
		ungroup() %>% 
		ggplot(aes(.value, color=factor(C))) + 
		geom_density() +
		scale_color_brewer(palette="Set1") +
		xlab("Slope") + 
		ylab("Probability density") +
		my_theme
) %>%
	cowplot::plot_grid(plotlist = .,		 align = "v",		 axis = "b",	 rel_widths = 1) %>%
	ggsave(
		filename =	"dev/proof_of_concept_generative_model.pdf",
		device = "pdf",
		useDingbats=FALSE,
		units = c("mm"),
		width = 280 ,
		height = 280 
	)

# # BRMS
# fit_brms <- brm(
# 	bind(V1 ,   V2 ,   V3   , V4  ,  V5) ~ risk,
# 	my_prop_dir_3 %>%
# 		select(sample, risk, cell_type, proportion) %>%
# 		spread(cell_type, proportion),
# 	family = 'dirichlet'
# )
# fit_brms %>%  gather_draws(b_muV2_risk, b_muV3_risk, b_muV4_risk,b_muV5_risk) %>% ggplot(aes(.value, color=.variable)) + geom_density()
