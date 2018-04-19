\name{ARMET_tc}
\alias{ARMET_tc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main script
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ARMET_tc(mix, my_design = NULL, cov_to_test = NULL, fully_bayesian = F, is_mix_microarray = F, observed_prop = NULL, ct_to_omit = c("t_CD4_naive", "adipocyte"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mix}{
%%     ~~Describe \code{mix} here~~
}
  \item{my_design}{
%%     ~~Describe \code{my_design} here~~
}
  \item{cov_to_test}{
%%     ~~Describe \code{cov_to_test} here~~
}
  \item{fully_bayesian}{
%%     ~~Describe \code{fully_bayesian} here~~
}
  \item{is_mix_microarray}{
%%     ~~Describe \code{is_mix_microarray} here~~
}
  \item{observed_prop}{
%%     ~~Describe \code{observed_prop} here~~
}
  \item{ct_to_omit}{
%%     ~~Describe \code{ct_to_omit} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (mix, my_design = NULL, cov_to_test = NULL, fully_bayesian = F, 
    is_mix_microarray = F, observed_prop = NULL, ct_to_omit = c("t_CD4_naive", 
        "adipocyte")) 
{
    if (file.exists("ARMET.ini")) {
        pars = read.ini("ARMET.ini")$ARMET
        writeLines("Importing parameters from .ini file..")
        for (n in names(pars)) {
            if (pars[n] \%in\% c("T", "F", "TRUE", "FALSE")) 
                pars[n] = as.logical(pars[n])
            assign(n, unlist(pars[n]))
        }
    }
    check_input()
    create_temp_result_directory()
    mix = as.matrix(mix)
    if (max(mix) < 50) 
        mix = exp(mix)
    if (!is_mix_microarray) 
        data(sigmatures_RNAseq)
    else data(sigmatures_array)
    data(tree_json)
    my_tree = drop_node_from_tree(tree, ct_to_omit)
    ref = ref[, colnames(ref) \%in\% get_leave_label(my_tree, last_level = 0, 
        label = "name")]
    common_genes = intersect(rownames(mix), rownames(ref))
    mix = mix[common_genes, , drop = FALSE]
    ref = ref[common_genes, , drop = FALSE]
    norm.obj = wrapper_normalize_mix_ref(mix, ref, is_mix_microarray)
    ref = norm.obj$ref
    mix = norm.obj$mix
    ggsave(sprintf("\%s/densities.png", output_dir), plot = norm.obj$plot)
    trees = lapply(1:ncol(mix), function(dummy) my_tree)
    names(trees) = colnames(mix)
    my_tree = run_coreAlg_though_tree(my_tree, list(mix = mix, 
        ref = ref, my_design = my_design, my_design.original = my_design.original, 
        cov_to_test = cov_to_test, fully_bayesian = fully_bayesian, 
        observed_prop = observed_prop, ct_to_omit = ct_to_omit, 
        my_tree = my_tree, is_mix_microarray = is_mix_microarray))
    if (!is.null(cov_to_test)) {
        data(treeYaml)
        tree.test = get_tree_hypoth_test(treeYaml, path = output_dir, 
            cov_to_test)
        print(tree.test, "significance_causal", "pvalue_causal", 
            "direction_causal", "significance_effectual", "pvalue_effectual", 
            "direction_effectual")
        save(tree.test, file = sprintf("\%s/tree_pvalues.RData", 
            output_dir))
    }
    trees = add_proportions_to_trees(trees, output_dir)
    trees = add_absolute_proportions_to_trees(trees)
    proportions = get_proportions_table(trees)
    markers = unique(do.call("c", lapply(dir(path = output_dir, 
        pattern = "model_in.RData", full.names = T), function(fn) {
        load(fn)
        rownames(model.in$x)
    })))
    signatures = list()
    signatures$orig = get_mean_signature(ref[markers, ], verbose = F)
    signatures$orig = signatures$orig[, colnames(signatures$orig) \%in\% 
        colnames(proportions)]
    signatures$predicted = NULL
    return(list(proportions = proportions, signatures = signatures, 
        mix = mix[markers, , drop = F]))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
