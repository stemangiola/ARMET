
#' convoluted_glm main
#'
#' @description The function for convoluted linear modelling takes as input a tidy table of feature count with three columns containing a sample ID, transcript ID and count, formula (continuous or discrete) and the covariate columns. The user can define a linear model with an input R formula, where the first covariate is the factor of interest. 
#'
#'
#' @param .data A tibble including a cell_group name column | sample name column | read counts column (optional depending on the input class) | covariate columns.
#' @param formula A formula. The formula describing the model for differential abundance, for example ~treatment.
#' @param .sample A column name as symbol. The sample identifier
#' @param .transcript A column name as symbol. The cell_group identifier
#' @param .abundance A column name as symbol. The cell_group abundance (read count). Used only for data frame count output. The variable in this column should be of class integer.
#' @param reference A data frame
#' @param tree A node object
#'
#' @param approximate_posterior A boolean
#' @param prior_survival_time A list
#' @param transform_time_function A function with nake survival time normally-shaped.

#' @return A nested tibble `tbl`, with the following columns
#' \itemize{
#'   \item cell_group - column including the cell groups being tested
#'   \item parameter - The parameter being estimated, from the design matrix dscribed with the input formula_composition and formula_variability
#'
#'   \item c_lower - lower (2.5%) quantile of the posterior distribution for a composition (c) parameter.
#'   \item c_effect - mean of the posterior distribution for a composition (c) parameter.
#'   \item c_upper - upper (97.5%) quantile of the posterior distribution fo a composition (c)  parameter.
#'   \item c_pH0 - Probability of the null hypothesis (no difference) for  a composition (c). This is not a p-value.
#'   \item c_FDR - False-discovery rate of the null hypothesis (no difference) for  a composition (c).
#'
#'   \item v_lower - (optional, present if variability is modelled dependent on covariates) lower (2.5%) quantile of the posterior distribution for a variability (v) parameter
#'   \item v_effect - (optional, present if variability is modelled dependent on covariates) mean of the posterior distribution for a variability (v) parameter
#'   \item v_upper - (optional, present if variability is modelled dependent on covariates) upper (97.5%) quantile of the posterior distribution for a variability (v) parameter
#'   \item v_pH0 - (optional, present if variability is modelled dependent on covariates) Probability of the null hypothesis (no difference) for a variability (v). This is not a p-value.
#'   \item v_FDR - (optional, present if variability is modelled dependent on covariates) False-discovery rate of the null hypothesis (no difference), for a variability (v).
#' }
#'
#' @examples
#'
#' data("test_mixture")
#' data("no_hierarchy_reference")
#'
#'  test_mixture |>
#'  convoluted_glm(
#'    ~ factor_of_interest,
#'    .sample = sample,
#'    .transcript = symbol,
#'    .abundance = count,
#'    reference = no_hierarchy_reference
#'   )
#'
#' @export
#'
#'
convoluted_glm = function(.data,
													.formula = ~ 1,
													.sample,
													.transcript,
													.abundance,
													reference = NULL,
													tree = NULL,
													
													# Secondary arguments
													approximate_posterior = F,
													prior_survival_time = c(),
													transform_time_function = sqrt,
													use_data = TRUE
												){
	
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	
	
	if(is.null(tree)){
		
		setup_convolved_lm_NON_hierarchical(
			.data,
			.formula = .formula,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			approximate_posterior = approximate_posterior,
			prior_survival_time = prior_survival_time,
			transform_time_function = transform_time_function,
			reference = reference
		) |>
			estimate_convoluted_lm(use_data)
	}
	else {
		
		setup_convolved_lm_hierarchical(
			.data,
			.formula = .formula,
			.sample = !!.sample,
			.transcript = !!.transcript,
			.abundance = !!.abundance,
			approximate_posterior = approximate_posterior,
			prior_survival_time = prior_survival_time,
			transform_time_function = transform_time_function,
			reference = reference,
			tree = tree
		) |>
			estimate_convoluted_lm_1() %>%
			estimate_convoluted_lm_2() %>%
			estimate_convoluted_lm_3() %>%
			estimate_convoluted_lm_4() 
	}
}