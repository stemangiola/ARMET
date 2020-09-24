#~/unix3XX/third_party_sofware/cctools-7.1.5-x86_64-centos6/bin/makeflow -T torque  --do-not-save-failed-output dev/TCGA_makeflow_pipeline/makefile_ARMET_TCGA.makeflow
#~/third_party_sofware/cctools-7.1.5-x86_64-centos6/bin/makeflow -T slurm  --do-not-save-failed-output dev/test_sim
#~/third_party_sofware/cctools-7.1.5-x86_64-centos6/bin/makeflow -T slurm  --do-not-save-failed-output dev/test_simulation_makeflow_pipeline/makefile_test_simulation.makeflow
library(tidyverse)

project_dir = "~/PhD/deconvolution/ARMET"
pipeline_dir = paste0(project_dir, "/dev/test_simulation_makeflow_pipeline")
output_dir = "dev/test_simulation" # paste0(project_dir, "/dev/test_simulation")
makelow_file = paste0(pipeline_dir, "/makefile_test_simulation.makeflow") 

prepend = function (x, values, before = 1) 
{
	n <- length(x)
	stopifnot(before > 0 && before <= n)
	if (before == 1) {
		c(values, x)
	}
	else {
		c(x[1:(before - 1)], values, x[before:n])
	}
}

expand_grid(
	slope = c(-2, -1, -.5, .5, 1, 2), 
	foreign_prop = c(0, 0.5, 0.8),
	S = c(30, 60, 90),
	which_changing = 1:16,
	run = 1:5,
	method = c("ARMET", "cibersort", "llsr", "epic")
)  %>%
	
	# Create input name
	mutate(input_file = sprintf("%s/input__slope_%s__foreignProp_%s__S_%s__whichChanging_%s__run_%s.rds", output_dir, slope, foreign_prop, S, which_changing, run)) %>%
	mutate(output_file = sprintf("%s/output__slope_%s__foreignProp_%s__S_%s__whichChanging_%s__run_%s__method_%s.rds", output_dir, slope, foreign_prop, S, which_changing, run, method)) %>%
	mutate(asses_file = sprintf("%s/asses__slope_%s__foreignProp_%s__S_%s__whichChanging_%s__run_%s__method_%s.rds", output_dir, slope, foreign_prop, S, which_changing, run, method)) %>%
	
	# Create input file
	{
		
		(.) %>%
			nest(meth = c(method, output_file, asses_file)) %>%
			mutate(command = sprintf(
				"%s:\n\tRscript %s/create_input.R %s %s %s %s %s %s", 
				input_file, pipeline_dir, slope, foreign_prop, S, which_changing, run, input_file
			)) %>% 
			pull(command) %>%
			
			# Resources
			prepend("CATEGORY=create_input\nMEMORY=8024\nCORES=2\nWALLTIME=1" ) %>%
			
			# Write to file
			write_lines(makelow_file, append = FALSE)
		
		(.)
	} %>%
	
	# Run ARMET
	{
		
		(.) %>%
			filter(method=="ARMET") %>%
			
			mutate(command = sprintf(
				"%s: %s\n\tRscript %s/infer.R %s %s %s", 
				output_file, input_file, pipeline_dir, input_file, output_file, method
			)) %>% 
			pull(command) %>%
			
			# Resources
			prepend("CATEGORY=ARMET\nMEMORY=40024\nCORES=12\nWALLTIME=1" ) %>%
			
			# Write to file
			write_lines(makelow_file, append = TRUE)
		
		(.)
	} %>%
	
	# Run third party
	{
		
		(.) %>%
			filter(method!="ARMET") %>%
			
			mutate(command = sprintf(
				"%s: %s\n\tRscript %s/infer.R %s %s %s", 
				output_file, input_file, pipeline_dir, input_file, output_file, method
			)) %>% 
			pull(command) %>%
			
			# Resources
			prepend("CATEGORY=third_party\nMEMORY=8024\nCORES=2\nWALLTIME=1" ) %>%
			
			# Write to file
			write_lines(makelow_file, append = TRUE)
		
		(.)
	} %>%
	
	# Asses
	{
		(.) %>%

			mutate(command = sprintf(
				"%s: %s\n\tRscript %s/asses_accuracy.R %s %s %s %s %s %s %s %s %s", 
				asses_file, output_file, pipeline_dir, foreign_prop ,run,S,slope,which_changing ,	method,input_file ,output_file , asses_file
			)) %>% 
			pull(command) %>%
				
			# Resources
			prepend("CATEGORY=asses\nMEMORY=8024\nCORES=2\nWALLTIME=1" ) %>%
			
			# Write to file
			write_lines(makelow_file, append = TRUE)
		
		(.)
	}

