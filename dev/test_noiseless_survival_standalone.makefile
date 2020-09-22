CORES=24
MEMORY=30024
dev/test_noiseless_survival_regression_4_slope.rds:
	dev/test_noiseless_survival_standalone.sh '4'
dev/test_noiseless_survival_regression_-4_slope.rds:
	dev/test_noiseless_survival_standalone.sh '-4'
dev/test_noiseless_survival_regression_2_slope.rds:
	dev/test_noiseless_survival_standalone.sh '2'
dev/test_noiseless_survival_regression_-2_slope.rds:
	dev/test_noiseless_survival_standalone.sh '-2'
dev/test_noiseless_survival_regression_1_slope.rds:
	dev/test_noiseless_survival_standalone.sh '1'
dev/test_noiseless_survival_regression_-1_slope.rds:
	dev/test_noiseless_survival_standalone.sh '-1'
dev/test_noiseless_survival_regression_0.5_slope.rds:
	dev/test_noiseless_survival_standalone.sh '0.5'
dev/test_noiseless_survival_regression_-0.5_slope.rds:
	dev/test_noiseless_survival_standalone.sh '-0.5'


