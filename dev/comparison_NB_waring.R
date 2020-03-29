load("dev/N52_mix.rda")
library(tidyverse)
library(ARMET)

N52_ARMET_T = ARMET_tc(
	mix = N52_mix,
	full_bayesian  = F,
	do_regression = T,
	cores = 20,
	levels = 3, iterations = 2000, sampling_iterations = 500
)

N52_ARMET_T$proportions %>%
	filter(level ==3) %>%
	filter(converged) %>%
	ggplot(aes(y = .value, x = sample)) + geom_point() + facet_wrap(~`Cell type category`)

save(N52_ARMET_T, file="dev/N52_ARMET_T.rda")


N52_ARMET_T_full = ARMET_tc(
	mix = N52_mix,
	full_bayesian  = T,
	do_regression = T,
	cores = 30,
	levels = 3
)

save(N52_ARMET_T_full, file="dev/N52_ARMET_T_full.rda")

N52_ARMET_T_full$proportions %>%
	filter(level ==3) %>%
	filter(converged) %>%
	ggplot(aes(y = .value, x = sample)) + geom_point() + facet_wrap(~`Cell type category`)


N52_ARMET_T_old = ARMET_tc(
	mix = N52_mix,
	full_bayesian  = F,
	do_regression = T,
	cores = 20,
	levels = 3
)

N52_ARMET_T_old$proportions %>%
	filter(level ==3) %>%
	filter(converged) %>%
	ggplot(aes(y = .value, x = sample)) + geom_point() + facet_wrap(~`Cell type category`)

N52_mix %>%
	gather(transcript, count, -sample) %>%
	ttBulk::deconvolve_cellularity(sample, transcript, count) %>%
	pivot_longer(names_prefix = "type: ", )
