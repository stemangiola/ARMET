#' Check whether there are NA counts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#' @param .data A tibble of read counts
#' @param .abundance A character name of the read count column
#'
#' @return A tbl
#'
error_if_counts_is_na = function(.data, .abundance) {
	.abundance = enquo(.abundance)

	# Do the check
	if (.data %>% filter(!!.abundance %>% is.na) %>% nrow %>% `>` (0))
		stop("tidyBulk says: You have NA values in your counts")

	# If all good return original data frame
	.data
}

#' Check whether there are duplicated genes/transcripts
#'
#' @import dplyr
#' @import tidyr
#' @import tibble
#'
#'
#' @param .data A tibble of read counts
#' @param .sample A character name of the sample column
#' @param .transcript A character name of the transcript/gene column
#' @param .abundance A character name of the read count column
#'
#' @return A tbl
error_if_duplicated_genes <- function(.data,
																			.sample = `sample`,
																			.transcript = `transcript`,
																			.abundance = `read count`) {
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)

	duplicates <-
		distinct( .data, !!.sample,!!.transcript,!!.abundance) %>%
		count(!!.sample,!!.transcript) %>%
		filter(n > 1) %>%
		arrange(n %>% desc())

	if (duplicates %>% nrow() > 0) {
		writeLines("Those are the duplicated genes")
		duplicates %>% print()
		stop(
			"tidyBulk says: Your dataset include duplicated sample/gene pairs. Please, remove redundancies before proceeding."
		)
	}

	.data

}


check_if_data_rectangular = function(.data, .sample, .transcript, .abundance){
	
	# Parse column names
	.sample = enquo(.sample)
	.transcript = enquo(.transcript)
	.abundance = enquo(.abundance)
	
	is_rectangular =
		.data %>%
		dplyr::distinct(!!.sample, !!.transcript, !!.abundance) %>%
		count(!!.sample) %>%
		count(n, name = "nn") %>%
		nrow %>%
		equals(1)
	
	is_rectangular
	
}

