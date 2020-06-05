# ~/unix3XX/third_party_sofware/cctools-7.1.5-x86_64-centos6/bin/makeflow  -T torque -B '-l walltime=148:60:00' --do-not-save-failed-output dev/test_noiseless_survival_standalone.makefile

# Create mix_base withn NA - imputed values

library(tidyverse)
library(magrittr)
library(purrr)
library(furrr)
library(data.tree)
library(foreach)
library(ARMET)
library(nanny)
args = commandArgs(trailingOnly=TRUE)
slope = as.numeric(args[1])

my_theme =
	theme_bw() +
	theme(
		panel.border = element_blank(),
		axis.line = element_line(),
		panel.grid.major = element_line(size = 0.2),
		panel.grid.minor = element_line(size = 0.1),
		text = element_text(size = 12),
		legend.position = "bottom",
		aspect.ratio = 1,
		strip.background = element_blank(),
		axis.title.x  = element_text(margin = margin(
			t = 10,
			r = 10,
			b = 10,
			l = 10
		)),
		axis.title.y  = element_text(margin = margin(
			t = 10,
			r = 10,
			b = 10,
			l = 10
		)),
		axis.text.x = element_text(angle = 90, hjust = 1)
	)



noiseles_test = function(mix_base, slope, which_changing, S) {
	
	cell_types =  mix_base %>% filter(level ==3) %>% pull(`Cell type category`) %>% unique
	
	alpha = ARMET:::get_alpha_test(slope, which_changing, cell_types)
	X_df = ARMET:::get_survival_X(S)
	
	mix = mix_base %>% ARMET:::generate_mixture(X_df, alpha)
	
	rr =
		mix %>%
		mutate(`count mix` = as.integer(`count mix`), run = as.character(run)) %>%
		select(-level) %>%
		ARMET_tc(
			~ censored(days, alive),
			run, symbol, `count mix`,
			prior_survival_time = X_df$real_days %>% as.numeric
		)  %>%
		ARMET_tc_continue(2) %>%
		ARMET_tc_continue(3)
	
	list(mix = mix, result = rr)
	#%>% saveRDS(sprintf("dev/test_student_noisless_%s", which_is_up_down %>% paste(collapse="_")))
}

mix_base = ARMET:::get_noiseless_harmonised()

which_is_up_down = 1:16

map(which_is_up_down,  ~ noiseles_test(mix_base, slope, .x, 30)) %>%
	saveRDS(sprintf("dev/test_noiseless_survival_regression_%s_slope.rds", slope))
