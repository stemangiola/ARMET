\name{ARMET_tc_coreAlg}
\alias{ARMET_tc_coreAlg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ARMET core alg}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ARMET_tc_coreAlg(obj.in, ct, real_prop_obj = NULL, is_test = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj.in}{
%%     ~~Describe \code{obj.in} here~~
}
  \item{ct}{
%%     ~~Describe \code{ct} here~~
}
  \item{real_prop_obj}{
%%     ~~Describe \code{real_prop_obj} here~~
}
  \item{is_test}{
%%     ~~Describe \code{is_test} here~~
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
function (obj.in, ct, real_prop_obj = NULL, is_test = F) 
{
    mix = obj.in$mix
    ref = obj.in$ref
    my_design = obj.in$my_design
    my_design.original = obj.in$my_design.original
    cov_to_test = obj.in$cov_to_test
    fully_bayesian = obj.in$fully_bayesian
    observed_prop = obj.in$observed_prop
    ct_to_omit = obj.in$ct_to_omit
    my_tree = obj.in$my_tree
    is_mix_microarray = obj.in$is_mix_microarray
    theta = apply(mix, 2, function(mc) length(which(mc > 0)))
    theta[theta == nrow(mix)] = theta[theta == nrow(mix)] - 1
    theta = theta/nrow(mix)
    theta = t(as.matrix(theta))
    ref.obj = node_to_balanced_sampling(my_tree, ct, ref)
    ref = cbind(do.call("cbind", lapply(ref.obj$ct_main, function(ro) ro$df)))
    markers = unique(unlist(lapply(ref.obj$ct_main, function(ro) get_genes(node_from_name(my_tree, 
        ro$name), recursive = F))))
    ref.cell_types = factor(colnames(ref))
    order_cell_types = levels(ref.cell_types)
    if (any(!markers \%in\% rownames(mix))) 
        warning(sprintf("ARMET: The following markers were not present into the mixture dataset: \%s Please reconsider the filtering process.\n", 
            paste(markers[!markers \%in\% rownames(mix)])))
    markers = markers[markers \%in\% rownames(ref) & markers \%in\% 
        rownames(mix)]
    mix = mix[markers, , drop = FALSE]
    ref = ref[markers, , drop = FALSE]
    ref.mean = get_mean_signature(ref)
    ref.mean = ref.mean[, order_cell_types]
    e.obj = prepare_input(ref, order_cell_types, markers)
    h = rev(rev(get_hierarchy(my_tree, ct))[-1])
    if (length(h) > 0) {
        for (i in 1:length(h)) {
            if (i == 1) 
                beta_bg = read.csv(sprintf("\%s_composition.csv", 
                  h[i]), row.names = 1)
            else {
                p = read.csv(sprintf("\%s_composition.csv", h[i]), 
                  row.names = 1)
                p = sweep(as.matrix(p), MARGIN = 1, beta_bg[, 
                  h[i]], `*`)
                beta_bg = beta_bg[, colnames(beta_bg) != h[i]]
                beta_bg = cbind(beta_bg, p)
            }
        }
        bg = cbind(do.call("cbind", lapply(ref.obj$ct_background, 
            function(ro) ro$df)))
        bg = bg[markers, ]
    }
    else {
        bg = matrix(rep(0, length(markers)))
        rownames(bg) = markers
        colnames(bg) = "nn"
        beta_bg = data.frame(rep(1, ncol(mix)), rep(0, ncol(mix)))
        colnames(beta_bg) = c(ct, "bg")
    }
    p_target = t(beta_bg[, ct, drop = FALSE])
    beta_bg = as.matrix(beta_bg[, !colnames(beta_bg) \%in\% ct])
    x_bg = get_mean_signature(bg)
    bg.obj = prepare_input(bg, levels(factor(colnames(bg))), 
        markers)
    y_hat_background = as.matrix(beta_bg) \%*\% t(x_bg)
    rownames(y_hat_background) = colnames(mix)
    e_map_matrix_dim = table(e.obj$e$ct_num)
    e_map_matrix_norm = median(e_map_matrix_dim)/e_map_matrix_dim
    e_i_matrix = do.call("rbind.fill", lapply(unique(sort(e.obj$e$ct_num)), 
        function(ct) data.frame(t(data.frame(rownames(e.obj$e[e.obj$e$ct_num == 
            ct, ]))))))
    e_i_matrix = apply(e_i_matrix, 2, function(mc) {
        x = as.numeric(as.character(mc))
        x[is.na(x)] = 0
        x
    })
    my_design = NULL
    my_local_design = if (is.null(my_design)) 
        matrix(rep(1, ncol(mix)), ncol = 1)
    else my_design
    model.in = list(G = length(markers), S = ncol(mix), P = length(unique(e.obj$e$ct_num)), 
        R = ncol(my_local_design), y = t(mix), X = my_local_design, 
        x = ref.mean, y_hat_background = y_hat_background, p_target = p_target, 
        theta = theta, is_mix_microarray = as.numeric(is_mix_microarray), 
        E = nrow(e.obj$e), E_MU = nrow(e.obj$e_mu), map_e_genes = e.obj$e$gene_num, 
        map_e_ct = e.obj$e$ct_num, map_e_to_mu = e.obj$e$map_to_mu, 
        e_ = e.obj$e$value, map_e_mu_gene = e.obj$e_mu$gene_num, 
        map_e_mu_ct = e.obj$e_mu$ct_num, e_i_matrix = e_i_matrix, 
        e_map_matrix_dim = e_map_matrix_dim, e_map_matrix_norm = e_map_matrix_norm, 
        B = nrow(bg.obj$e), B_MU = nrow(bg.obj$e_mu), Q = length(unique(bg.obj$e$ct_num)), 
        map_bg_genes = bg.obj$e$gene_num, map_bg_ct = bg.obj$e$ct_num, 
        map_bg_to_mu = bg.obj$e$map_to_mu, b = bg.obj$e$value, 
        map_bg_mu_gene = bg.obj$e_mu$gene_num, map_bg_mu_ct = bg.obj$e_mu$ct_num, 
        beta_bg = beta_bg)
    save(model.in, file = sprintf("\%s_model_in.RData", ct))
    model = if (fully_bayesian) 
        stanmodels$ARMET_tc_recursive
    else stanmodels$ARMET_tcFix_recursive.stan
    fit = rstan::sampling(model, data = model.in, iter = 1000, 
        control = list(adapt_delta = 0.99, stepsize = 0.01, max_treedepth = 15))
    proportions = t(parse_result_2D(fit, "beta"))
    proportions_sd = t(parse_result_2D(fit, "beta", 3))
    proportions_2.5 = t(parse_result_2D(fit, "beta", 4))
    proportions_97.5 = t(parse_result_2D(fit, "beta", 8))
    rownames(proportions) = rownames(proportions_2.5) = rownames(proportions_97.5) = rownames(proportions_sd) = colnames(mix)
    colnames(proportions) = colnames(proportions_2.5) = colnames(proportions_97.5) = colnames(proportions_sd) = order_cell_types
    save(fit, file = sprintf("\%s_fit.RData", ct))
    write.csv(proportions, sprintf("\%s_composition.csv", ct))
    proportions_4_plot = as.data.frame(proportions)
    proportions_4_plot$sample = rownames(proportions_4_plot)
    proportions_4_plot$cov = if (length(cov_to_test) > 0) 
        my_design.original[, "cov"]
    else 1
    proportions_4_plot = melt(proportions_4_plot, id.vars = c("sample", 
        "cov"))
    proportions_4_plot$perc_1 = melt(proportions_2.5)$value
    proportions_4_plot$perc_2 = melt(proportions_97.5)$value
    proportions_4_plot$minus_sd = proportions_4_plot$value - 
        melt(proportions_sd)$value
    proportions_4_plot$plus_sd = proportions_4_plot$value + melt(proportions_sd)$value
    p = ggplot(proportions_4_plot, aes(x = jitter(cov), y = value, 
        fill = factor(variable))) + geom_boxplot(coef = 6) + 
        geom_point(size = 0.1) + geom_linerange(aes(ymin = perc_1, 
        ymax = perc_2), alpha = 0.05) + geom_errorbar(aes(ymin = minus_sd, 
        ymax = plus_sd), width = 0.2, alpha = 0.2) + facet_grid(~variable) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.4, 
            hjust = 1), panel.background = element_blank())
    ggsave(sprintf("\%s_proportions.pdf", ct), useDingbats = FALSE, 
        plot = p)
    print("saved")
    if (length(cov_to_test) > 0) {
        gq <- as.matrix(fit, pars = c("beta_gen"))[1, ]
        gq = parse_summary_vector_in_2D(gq)
        colnames(gq) = colnames(mix)
        rownames(gq) = colnames(model.in$x)
        gq = as.data.frame(t(gq))
        gq$cov = my_design[, "cov"]
        gq$sample = rownames(gq)
        gq = melt(gq, id.vars = c("sample", "cov"))
        p = ggplot(gq, aes(x = factor(cov), y = value, fill = factor(variable))) + 
            geom_boxplot() + geom_jitter(position = position_dodge(width = 0.75), 
            aes(group = variable)) + facet_grid(~variable) + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.4, 
                hjust = 1))
        ggsave(sprintf("\%s_proportions_posterior.pdf", ct), useDingbats = FALSE, 
            plot = p)
    }
    print("saved")
    df_4_plot = data.frame(merge(melt(mix + 1), melt(t(proportions \%*\% 
        t(ref.mean) + 1)), by = c("X1", "X2")))
    colnames(df_4_plot)[3:4] = c("value", "predicted")
    df_4_plot$ct = sapply(df_4_plot$X1, function(x1) node_from_gene(my_tree, 
        x1)$name)
    df_4_plot$bg = apply(df_4_plot, 1, function(mr) y_hat_background[mr[2], 
        mr[1]])
    df_4_plot$bg_prop = df_4_plot$bg/df_4_plot$value
    df_4_plot$bg_prop[df_4_plot$bg_prop > 10] = 10
    df_4_plot = df_4_plot[df_4_plot$X2 \%in\% unique(df_4_plot$X2)[1:20], 
        ]
    p = ggplot(df_4_plot, aes(predicted, value, label = X1, size = bg_prop, 
        color = factor(ct))) + geom_point(alpha = 0.5) + scale_size(range = c(0, 
        5)) + expand_limits(x = 0.1, y = 0.1) + geom_abline(intercept = 0, 
        slope = 1) + scale_x_log10() + scale_y_log10() + geom_text_repel(aes(color = ct), 
        size = 2, segment.alpha = 0.2) + facet_wrap(~X2, nrow = 5) + 
        coord_fixed() + theme(panel.background = element_blank(), 
        legend.text = element_text(size = 30))
    ggsave(sprintf("\%s_predicted_values.pdf", ct), useDingbats = FALSE, 
        width = 80, height = 80, units = "cm", plot = p)
    print("saved")
    if (!is.null(observed_prop) & 0) {
        writeLines("ARMET: start plotting expected gene vallues")
        link = as_tibble(do.call("rbind", lapply(unique(colnames(observed_prop)), 
            function(op) {
                ct = get_hierarchy(my_tree, ct = op)
                c(op, ct[ct \%in\% colnames(ref.mean)])
            })))
        if (ncol(link) >= 2) {
            colnames(link) = c("ct", "ct_new")
            observed_prop = as_tibble(melt(observed_prop))
            colnames(observed_prop) = c("sample", "ct", "value")
            observed_prop = left_join(observed_prop, link, by = "ct") \%>\% 
                select(-ct) \%>\% group_by(sample, ct_new) \%>\% 
                summarise(value = sum(value, na.rm = TRUE)) \%>\% 
                spread(ct_new, value)
            for (r in colnames(ref.mean)[!colnames(ref.mean) \%in\% 
                colnames(observed_prop)]) {
                my_df = matrix(rep(0, nrow(observed_prop)))
                colnames(my_df) = r
                observed_prop = observed_prop \%>\% bind_cols(as_tibble(my_df))
            }
            observed_prop = as.data.frame(observed_prop)
            rownames(observed_prop) = observed_prop[, 1]
            observed_prop = observed_prop[, -1]
            observed_prop = as.matrix(observed_prop)
            observed_prop = observed_prop[, colnames(ref.mean)]
            print(proportions)
            df_4_plot = data.frame(merge(melt(mix + 1), melt(t(observed_prop \%*\% 
                t(ref.mean) + 1)), by = c("X1", "X2")))
            colnames(df_4_plot)[3:4] = c("value", "predicted")
            df_4_plot$ct = sapply(df_4_plot$X1, function(x1) node_from_gene(my_tree, 
                x1)$name)
            df_4_plot$bg = apply(df_4_plot, 1, function(mr) y_hat_background[mr[2], 
                mr[1]])
            df_4_plot$bg_prop = df_4_plot$bg/df_4_plot$value
            df_4_plot$bg_prop[df_4_plot$bg_prop > 10] = 10
            p = ggplot(df_4_plot[df_4_plot$X2 \%in\% unique(df_4_plot$X2)[1:20], 
                ], aes(predicted, value, label = X1, size = bg_prop)) + 
                geom_point(alpha = 0.5) + scale_size(range = c(0, 
                5)) + expand_limits(x = 0.1, y = 0.1) + geom_abline(intercept = 0, 
                slope = 1) + scale_x_log10() + scale_y_log10() + 
                geom_text_repel(aes(color = ct), size = 2, segment.alpha = 0.2) + 
                facet_wrap(~X2, nrow = 5) + coord_fixed() + theme(panel.background = element_blank(), 
                legend.text = element_text(size = 30))
            ggsave(sprintf("\%s_predicted_values_observed_prop.pdf", 
                ct), useDingbats = FALSE, width = 80, height = 80, 
                units = "cm", plot = p)
        }
    }
    list(proportions = proportions)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
