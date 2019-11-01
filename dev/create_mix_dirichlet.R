# Create mix_base withn NA - imputed values

library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)
source("~/PhD/deconvolution/ARMET/R/utils.R")
source("~/PostDoc/ppcSeq/R/do_parallel.R")
n_cores = 20
S = 30

logsumexp <- function (x) {
	y = max(x)
	y + log(sum(exp(x - y)))
}

softmax <- function (x) {
	exp(x - logsumexp(x))
}

mix_base = readRDS("dev/mix_base.RDS") %>% filter(level==3)


#rgamma(mix_base %>% distinct(`Cell type category`) %>% nrow, 1.1, 0.5)
#rnorm(mix_base %>% distinct(`Cell type category`) %>% nrow, 0, 2)
sample(0:1, mix_base %>% distinct(`Cell type category`) %>% nrow,replace = T)

intercept = c(-0.9516997,  0.3748425,  0.2348878, -0.4229247, -0.9586268, -3.1402425,  0.3323535, -0.7435126,  2.9877128, -1.2978599, 4.3127303,  0.5863942,  0.3885436,  1.3248586)
intercept = intercept - sum(intercept)/length(intercept)

alpha = matrix(
	intercept %>%
		c(0 ,0 ,0, 2, 0, 2, 0, 0, 2, 0, 0, 2, 0, 2),
	ncol = 2
)

X = matrix(
	rep(1, S) %>%
	c(rep(0, S/2)) %>%
		c(rep(1, S/2)),
	ncol = 2
)

# Choose samples
foreach(r = 1:S, .combine = bind_rows) %do% {
	mix_base %>%
		distinct(`Cell type category`, sample) %>%
		group_by(`Cell type category`) %>%
		sample_n(1) %>%
		ungroup() %>%
		mutate(run = r)
} %>%

# Choose proportions
	left_join(

		X %*% t(alpha) %>%
			apply(1, softmax) %>%
			t %>%
			`*` (40) %>%
			as.data.frame() %>%
			as_tibble() %>%
			setNames( mix_base %>% distinct(`Cell type category`) %>% pull(1) ) %>%
			mutate(run = 1:n()) %>%
			gather(`Cell type category`, alpha, -run)
	) %>%
	group_by(run) %>%
	mutate(p = gtools::rdirichlet(1, alpha)) %>%
	ungroup() %>%

	# Add counts
	left_join(mix_base, by=c("Cell type category", "sample")) %>%
	write_csv("dev/mix_dirichlet_source.csv")

# make mix
read_csv("dev/mix_dirichlet_source.csv") %>%
	mutate(c = `read count normalised bayes` * p) %>%
	group_by(run, symbol) %>%
	summarise(`count mix` = c %>% sum) %>%
	ungroup %>%
	write_csv("dev/mix_dirichlet.csv")

save(X, file="dev/mix_dirichlet_X.rda")
