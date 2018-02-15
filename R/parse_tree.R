#' Add info to tree
add_info_to_tree = function(node, ct, label, value){
	
	if(is.na(ct)) NULL
	else if(node$name==ct) {
		node[label] = list(value)
		node
	}
	else if(length(node$children)>0){
		node$children = lapply(node$children, add_info_to_tree, ct, label, value)
		node
	}
	else node
}

#' Drop a node from the tree
drop_node_from_tree = function(node, ct){
	if( !node$name%in%ct) {
		node$children = Filter(Negate(is.null), 
													 lapply(node$children, drop_node_from_tree, ct)
		)
		node
	}
}

#' Put proportions of cell types into the tree
prop_array_to_tree = function(tree, p_array){
	for(nn in names(p_array)) tree = add_info_to_tree(tree, nn, "absolute_proportion", p_array[nn])
	add_absolute_proportions_to_tree_from_leaves(tree)
}

#' Calculate absolute proportions from relative proportions
add_absolute_proportions_to_tree_from_leaves = function(node){
	
	if(length(node$children) == 0){
		if(length(node$absolute_proportion)==0) node$absolute_proportion = NA
	}
	else {
		node$children = lapply(node$children, add_absolute_proportions_to_tree_from_leaves)
		sum_array = unlist(lapply(node$children, function(nn) nn$absolute_proportion))
		#If middle node has result
		if(length(node$absolute_proportion)==0) node$absolute_proportion = if(length(na.omit(sum_array))==0) NA else  sum(na.omit(sum_array))
		# Through error if botha node have a sum and children have sum
		else if(length(na.omit(sum_array))>0) stop(sprintf("ARMET: in calculating proportion in tree results that both children and ancestor have been given proportions which represent a contraddiction. The node is %s", node$children))
	}
	node
}

#' Put proportions on tree
add_proportions_to_trees = function(trees, path="."){
	# Add root proportion
	for(i in 1:length(trees)) trees[[i]]$relative_proportion = 1
	
	for(l in get_node_names(trees[[1]], last_level=-1)){
		#print(l)
		p = read.csv(sprintf("%s/%s_composition.csv", path, l), row.names=1)
		for(i in 1:length(trees)){
			for(ll in colnames(p)){
				trees[[i]] = add_info_to_tree(trees[[i]], ll, "relative_proportion", p[i, ll])
			}
		}
	}
	trees
}

#' Put absolute proportions on tree
add_absolute_proportions_to_tree = function(node, p_ancestor = 1){
	node$absolute_proportion = node$relative_proportion * p_ancestor
	if(length(node$children)>0){
		node$children = lapply(node$children, add_absolute_proportions_to_tree, node$absolute_proportion)
	}
	node
}

#' Put absolute proportions on trees
add_absolute_proportions_to_trees = function(trees){
	lapply(trees, add_absolute_proportions_to_tree)
}

#' Put absolute proportions on trees
get_ancestor_ct = function(node, ct){
	rev(get_hierarchy(node, ct))[2]
}

#' Get proportions from tree
get_proportions_array = function(node, sample, label="absolute_proportion"){
	
	p = as.data.frame(get_leave_label(node, last_level = 1, label = "absolute_proportion"))
	rownames(p) = get_leave_label(node, last_level = 1, label = "name")
	colnames(p) = sample
	p
	
}

#' Get proportions table from trees
get_proportions_table = function(trees){
	t(do.call("cbind", lapply(names(trees), function(na) get_proportions_array(trees[na][[1]], na) )))
}

#' Get nome from cell type
node_from_name = function(node, ct){
	
	if(is.na(ct)) NULL
	else if(node$name==ct) node
	else if(length(node$children)>0)
		do.call("c", lapply(node$children, function(n){
			x = node_from_name(n, ct)
			if(!is.null(x)) x
			else NULL
		})
		)
	else NULL
}

#' Get nome from gene
node_from_gene = function(node, gene){
	if(is.na(gene)) NULL
	else if(any(node$markers==gene)) node
	else if(length(node$children)>0)
		do.call("c", lapply(node$children, function(n){
			x = node_from_gene(n, gene)
			if(!is.null(x)) x
			else NULL
		})
		)
	else NULL
}

#' Get information from tree recursively for every node
get_node_label_recursive = function(node, last_level = 0, label = "name"){
	if(length(node$children)>0){
		if(last_level<=0) c(
			unlist(node[label]),
			do.call("c", lapply(node$children, get_node_label_recursive, last_level, label))
		)
		else if(last_level==1) do.call("c", lapply(node$children, get_node_label_recursive, last_level, label))
	}
	else if(last_level>=0) unlist(node[label])
}

get_last_existing_leaves_with_annotation = function(node, label = "absolute_proportion"){
	if(length(node$children)>0){
		
		do.call("c", lapply(node$children,function(n){
			res = get_last_existing_leaves_with_annotation(n, label)
			if(length(res)==0) {
				res = unlist( n[label])
				if(length(res)>0) {
					names(res) = n$name
				}
			}
			res
		}))
		
	}	else {
		res = unlist( node[label])
		if(length(res)>0) {
			names(res) = node$name
			res
		} else NULL
	}
}

#' Get leaves cell type names
get_leave_names = function(node, recursive = T, last_level = 0){
	
	if(recursive){
		ll =  get_node_label_recursive(node, last_level)
		ll = ll[ll!=node$name]
		ll
	}
	else do.call("c", lapply(node$children, function(n) n$name))
}

#' Get any info from the tree
get_leave_label = function(node, recursive = T, last_level = 0, label = "name"){
	
	if(recursive){
		ll =  get_node_label_recursive(node, last_level, label)
		ll = ll[ll!=node$name]
		ll
	}
	else do.call("c", lapply(node$children, function(n) n[label]))
}

#' Get cell type names
get_node_names = function(node, recursive = T, last_level = 0){
	if(recursive) get_node_label_recursive(node, last_level)
	else do.call("c", lapply(node$children, function(n) n$name))
}

#' Get cell type names from a specific level of the tree
get_node_label_level_specfic = function(node, label = "name", recursive = T, level= 0, start_level = 0, stop_level = 0){
	if(length(node$children)>0 & recursive & (level < stop_level | stop_level == 0)){
		if(stop_level > 0) level = level + 1
		my_mark = do.call("c", lapply(node$children, get_node_label_level_specfic, label, recursive, level, start_level, stop_level))
		if(level > start_level) my_mark = c(node[label], my_mark)
		return(unlist(my_mark))
	}
	else return(unlist(node[label]))
}

#' Get genes from tree
get_genes = function(node, label = "markers", recursive = T, level= 0, start_level = 0, stop_level = 0){
	get_node_label_level_specfic(
		node, 
		label = label, 
		recursive = recursive, 
		level= level, 
		start_level = start_level, 
		stop_level = stop_level
	)
}

#' Get prop table with proportions summarising at a selected level from a tree
tree_to_leveled_prop = function(tree, level = 999){
	start_level = stop_level = level
	df = data.frame(t(get_node_label_level_specfic(tree, label="absolute_proportion", start_level = start_level, stop_level = stop_level)))
	colnames(df) = get_node_label_level_specfic(tree, start_level = start_level, stop_level = stop_level)
	as_tibble(df)
}

#' Get prop table with proportions summarising at a selected level from another table
prop_table_to_leveled_prop_table = function(tree, df, level = 999){
	cn = df[,1]
	df = df[,-1] 
	df = do.call("rbind" , apply(df , 1, function(mr)
		tree_to_leveled_prop(prop_array_to_tree(tree, mr ), level = level)
	))
	
	bind_cols(cn, df)
}

#' Get prop table with proportions summarising at a selected level from trees
prop_treeList_to_leveled_prop_table = function(treeList, level = 999){
	
	df = do.call("bind_rows" , lapply(treeList, function(tree)
		tree_to_leveled_prop(tree, level)
	))
	bind_cols(tibble(X1 = names(treeList)), df)
	
}

#' Get cell type hierarchy of a node
get_hierarchy = function(node, ct){
	if(is.na(ct)) NULL
	else if(node$name==ct) node$name
	else if(length(node$children)>0){
		child = do.call("c", lapply(node$children, function(n){
			get_hierarchy(n, ct)
		}))
		
		if(length(child)>0) c(node$name, child)
		else NULL
	}
	else NULL
}

#' Get balanced sampling of samples for every node acros all its constituents
node_to_balanced_sampling = function(tree, ct, df, keep_orig_names=F, do_replacement = F){

	ct_array = levels(df$ct)
	
	get_leave_names_recursive = function(node){
		
		#df %>% dplyr:::distinct(sample, ct) %>% dplyr:::filter(ct==node$name)
		
		my_cnt = length(which(ct_array==node$name))
		my_sampling = which(ct_array==node$name)
		
		if(length(node$children)>0) {
			children = lapply(node$children, get_leave_names_recursive)
			
			all_cnt = do.call("c", lapply(children, function(ch) ch$cnt))
			
			# Don't get teh min but the second min
			if(length(all_cnt)==1) min_cnt = all_cnt
			else min_cnt = round(mean(sort(all_cnt)[1:2]))
		
			sum_cnt = min_cnt * length(all_cnt)
			children = lapply(children, function(ch) {
				ch$cnt = min_cnt
				ch
			})
			
			my_sampling = c(
				my_sampling,
				do.call("c", lapply(children, function(ch){
					
					if(length(ch$sampling)==0) stop(sprintf("ARMET: %s not present in the data set.", ch$name))
					
					ch$sampling 
					
				}))
			)
			
			list(name = node$name, cnt = sum(c(my_cnt, sum_cnt)), sampling = my_sampling, children = children )
			
		}
		else list(name = node$name, cnt = my_cnt, sampling = my_sampling)
	}
	
	get_ct_4_model = function(ct){
		
		# get childred
		ct_main = get_leave_names(node_from_name(tree, ct), recursive=F)
		
		# get backgrounbd
		get_ct_4_model_ancestor = function(ct){
			ct_childrens = get_leave_names(node_from_name(tree, ct), recursive=F)
			ct_ancestor = rev(get_hierarchy(tree, ct))[2]
			if(length(ct_ancestor)>0) {
				ancestor_leaves = get_ct_4_model_ancestor(ct_ancestor)
				c(ct_childrens, ancestor_leaves[ancestor_leaves != ct])
			}
		}
		ct_background = get_ct_4_model_ancestor(ct)
		ct_background = ct_background[!ct_background%in%ct_main]
		
		list(
			ct_main = ct_main,
			ct_background = ct_background
		)
		
	}
	
	cts = get_ct_4_model(ct)
	cts$ct_main= lapply(cts$ct_main, function(cm){
		browser()
		sampling = get_leave_names_recursive(node_from_name(tree, cm))$sampling
		
		my_df = df[,sampling, drop=F]
		if(!keep_orig_names) colnames(my_df) = rep(cm, length(sampling))
		
		list(
			name = cm,
			sampling = sampling,
			df = my_df
		)
	})
	
	cts$ct_background= lapply(cts$ct_background, function(cm){
		sampling = get_leave_names_recursive(node_from_name(tree, cm))$sampling
		my_df = df[,sampling, drop=F]
		if(!keep_orig_names) colnames(my_df) = rep(cm, length(sampling))
		
		list(
			name = cm,
			sampling = sampling,
			df = my_df
		)
	})
	
	cts
}

#' Plot violin plot for markers in a node
plot_makers = function(tree, ct, df, do_replacement=F, markers = NULL){
	
	bg.obj =  node_to_balanced_sampling(tree, ct, df, keep_orig_names = F,do_replacement=do_replacement)$ct_background
	if(is.null(markers)) markers = get_genes(node_from_name(tree, ct), recursive=F)
	bg =        cbind(do.call("cbind", lapply(bg.obj, function(ro) ro$df)))
	
	ref.obj = node_to_balanced_sampling(tree, rev(get_hierarchy(tree, ct))[2], df, keep_orig_names = F, do_replacement=do_replacement)$ct_main
	ref =     cbind(do.call("cbind", lapply(ref.obj, function(ro) ro$df)))
	ref = ref[,colnames(ref)==ct]
	
	my_df = cbind(bg, ref)
	if(!all(markers%in%rownames(my_df))) warning(sprintf("ARMET: the gene %s is not in the ref data frame.", markers[!markers%in%rownames(my_df)]))
	markers = markers[markers%in%rownames(my_df)]
	if(length(markers) <= 1) {
		my_df = as.matrix(my_df[markers,]+1)
		my_df = as.data.frame(my_df)
		my_df$X2 = rownames(my_df)
		my_df$X1 = markers
		colnames(my_df)[1] = "value"
	}
	else {
		my_df = as.matrix(my_df[markers,]+1)
		my_df = reshape:::melt(my_df)
		colnames(my_df) = c("X1", "X2", "value")
	}
	
	ggplot2:::ggplot(my_df, ggplot2:::aes(x=factor(X2), y=value,fill=factor(X2)))+
		ggplot2:::geom_violin()+ 
		ggplot2:::geom_jitter(position=position_dodge(width=0.75), ggplot2:::aes(group=X2)) +
		ggplot2:::facet_grid(.~X1) + 
		ggplot2:::theme(axis.text.x= ggplot2:::element_text(angle=90, vjust=0.4,hjust=1)) +
		ggplot2:::scale_y_log10() + 
		ggplot2:::ggtitle(ct) +
		ggplot2:::facet_wrap(~X1,  nrow = ceiling(length(markers)/6)) 
	
}

#' Plot violin plot for arbitrary markers 
plot_custom_makers = function(df){
	
	ggplot2:::ggplot(reshape:::melt(df), aes(x=factor(X2), y=value,fill=factor(X2)))+
		geom_violin()+ geom_jitter(position=position_dodge(width=0.75),aes(group=X2)) +
		facet_grid(.~X1) + theme(axis.text.x=element_text(angle=90, vjust=0.4,hjust=1)) +
		scale_y_log10() 
	
}

#' Plot violin plot for a node downstream
plot_makers_in_tree = function(tree, df, do_replacement=F){
	lapply(get_leave_names(tree), function(l) {
		plot(plot_makers(tree, l, df, do_replacement=do_replacement))
		browser()
	})
}

#' Get top differentially transcribed genes
top_de = function(df, design){
	y <- edgeR:::DGEList(counts=round(as.matrix(df)), group=design[,2])
	#y <- calcNormFactors(y, method="upperquartile")
	y <- edgeR:::estimateGLMCommonDisp(y, design)
	y <- edgeR:::estimateGLMTagwiseDisp(y, design)
	fit <- edgeR:::glmFit(y, design)
	lrt <- edgeR:::glmLRT(fit, coef=2)
	edgeR:::topTags(lrt, n=nrow(df))$table
}

#' Test differentially transcribed gene for a node
test_genes_in_node = function(tree, df, ct, n = 20, do_replacement=F, ct_to_test_against = NULL, expression_threshold=100, include_bg = F, quantile_norm = T, na_rm = F, do_higher=T, do_calc_sd = T){
	
	ancestor = get_ancestor_ct(tree, ct)
	if(!is.null(ancestor)){
		
		if(na_rm) df = na.omit(df)
		replace_NA_for_edgeR = function(my_df){
			temp =                     data.frame(ct = colnames(my_df), t(my_df), row.names=NULL)
			my_df =                  do.call("rbind", by(temp, temp$ct, function(df) 
			{
				
				ct = df[,1]
				df = apply(df[,-1, drop=F], 2, function(mc){
					if(all(is.na(mc))) mc = 0.
					if(any(is.na(mc))) mc[is.na(mc)] = exp(rnorm(length(which(is.na(mc))), mean(log(mc+1), na.rm=T), sd(log(mc+1), na.rm=T)))
					mc
				})
				if(typeof(df)!="double") df = do.call("cbind", df)
				
				# If jst one row
				if(is.null(dim(df))) df = t(as.data.frame(df))
				rownames(df) = ct
				df
			}))
			t(my_df)
		}
		
		calc_sd = function(df){
			cn = colnames(df)
			apply(df, 1, function(mr) max(stats:::aggregate(log(mr+1), by= list(cn), FUN=sd)$x))
		}
		
		my_df = cbind(do.call("cbind", lapply(node_to_balanced_sampling(tree, ancestor, df, keep_orig_names = F, do_replacement)$ct_main, function(ro) ro$df)))
		
		if(include_bg)
			my_df = cbind(my_df, cbind(do.call("cbind", lapply(node_to_balanced_sampling(tree, ancestor, df, keep_orig_names = F, do_replacement)$ct_background, function(ro) ro$df))))
		
		# Check within
		if(length(node_from_name(tree, ct)$children) > 0){
			my_df.internal = cbind(do.call("cbind", lapply(node_to_balanced_sampling(tree, ancestor, df, keep_orig_names = T, do_replacement)$ct_main, function(ro) ro$df)))
			my_df.internal = my_df.internal[,colnames(my_df.internal)%in%get_leave_names(node_from_name(tree, ct))]
			my_df.internal = replace_NA_for_edgeR(my_df.internal)
			
			design.internal <- model.matrix(~ l, data=data.frame(l = colnames(my_df.internal)))
			fit <- limma:::lmFit(log(my_df.internal+1), design.internal)
			fit <- limma:::eBayes(fit)
			top = limma:::topTable(fit, n=999999)
			gene_to_exclude = rownames(top[top$adj.P.Val<0.05,])
			
			# if no genes remained take the best 500
			at_least = nrow(my_df)/4
			if(nrow(my_df)-length(gene_to_exclude)<at_least) gene_to_exclude = gene_to_exclude[1:(length(gene_to_exclude)-at_least)]
			
			my_df = my_df[!rownames(my_df)%in%gene_to_exclude,]
			#cl <- makeCluster(30)
			# clusterEvalQ(cl, library(reshape))
			# clusterEvalQ(cl, library(stats))
			# pvalues = parApply(cl, my_df.internal, 1, function(rn) summary(aov(value ~ Var1, data = reshape:::melt(log(as.matrix(rn)+1))))[[1]][1,5] )
			# pvalues = apply(my_df.internal, 1, function(rn) summary(aov(value ~ Var1, data = reshape:::melt(log(as.matrix(rn)+1))))[[1]][1,5] )
			# pvalues = p.adjust(pvalues, method="BH")
			#gene_to_exclude = names(pvalues[pvalues<0.05])
		}
		print(sprintf("%s genes after homogeneous filtering", nrow(my_df)))
		
		# If specific comparison
		if(!is.null(ct_to_test_against)) my_df = my_df[,colnames(my_df)%in%c(ct, ct_to_test_against)]
		print(sprintf("%s genes after threshold filtering", nrow(my_df)))
		
		# Eliminate the lowly expressed genes for markes
		my_df = my_df[
			matrixStats:::rowMedians(my_df[,colnames(my_df)==ct], na.rm=T)>expression_threshold
			& !is.na(matrixStats:::rowMedians(my_df[,colnames(my_df)==ct], na.rm=T))
			,]
		
		
		# Replace 
		my_df = replace_NA_for_edgeR(my_df)
		rownames(my_df) = sapply(rownames(my_df), function(m) gsub(".", "-", m, fixed=T))
		
		# External comparison
		design <- model.matrix(~ is_l, data=data.frame(is_l = colnames(my_df)==ct))
		top = top_de(my_df, design)
		top = top[top$FDR<0.05,]
		top = if(do_higher) top[top$logFC>0,] else top[top$logFC<0,]
		writeLines(sprintf("ARMET: %s DE genes", nrow(top)))
		
		top = head(top, n=n*4)  #n=max(nrow(top)/2, n))
		top = top[order(top$logFC, decreasing = T),]
		top = head(top, n=n*2)  #max(nrow(top)/2, n))
		if(do_calc_sd){
			top = top[order(calc_sd(my_df[rownames(top),])),]
		}
		top = head(top, n=n)
		
		markers = rownames(top)
		print(sort(markers))
		
		plot_makers(tree, ct, df, markers=markers)
		
		# df_4_plot  = reshape:::melt(my_df[markers,]+1)
		# df_4_plot$X1 = factor(df_4_plot$X1, levels=unique(df_4_plot$X1)) 
		# colnames(df_4_plot) = c("Var1", "Var2", "value")
		# p = ggplot( df_4_plot, aes(x=factor(Var2), y=value,fill=factor(Var2)))+
		# 	geom_violin()+ geom_jitter(position=position_dodge(width=0.75),aes(group=Var2)) +
		# 	facet_wrap(~Var1,  nrow = ceiling(n/20)) + theme(axis.text.x=element_text(angle=90, vjust=0.4,hjust=1)) +
		# 	scale_y_log10() + ggtitle(ct)
		# plot(p)
	}
}

#' Get top differentially transcribed in all tree
test_genes_in_tree = function(tree, df, n = 20, do_replacement=F, expression_threshold=100, na_rm = F){
	lapply(get_leave_names(tree), function(l) {
		test_genes_in_node(tree, df, l, n, expression_threshold=expression_threshold, na_rm=na_rm)
		browser()
	})
}

#' Run ARMET core algorithm recursively on the tree 
#'
#' @param node A node from the tree
#' @param obj.in A object with the input
#' @return A probability array
run_coreAlg_though_tree_recursive = function(node, obj.in, bg_tree){
	
	# library(doFuture)
	# doFuture:::registerDoFuture()
	# future:::plan(future:::multiprocess)

	
	if(length(node$children)>0){
		
		obj.in$bg_tree = bg_tree
		obj.out = ARMET_tc_coreAlg(obj.in, node$name)
		
		for(ll in colnames(obj.out$proportion)){
			node = add_info_to_tree(node, ll, "relative_proportion", obj.out$proportion[, ll, drop=F])
			bg_tree = add_info_to_tree(bg_tree, ll, "relative_proportion", obj.out$proportion[, ll, drop=F])
		}
		
		if(obj.in$multithread){
			#n_cores = floor(detectCores() / 6)
			n_cores = 4
			cl <- parallel:::makeCluster(n_cores)
			parallel:::clusterExport(cl, c("obj.in", "bg_tree"), environment())
			#clusterEvalQ(cl, library("ARMET"))
			doParallel:::registerDoParallel(cl)
			
			node$children = foreach:::foreach(cc = node$children) %do% {
				run_coreAlg_though_tree_recursive(cc, obj.in, bg_tree)
			}
			
			parallel:::stopCluster(cl)
			
		} else {
			node$children = lapply(node$children, function(cc){
				run_coreAlg_though_tree_recursive(cc, obj.in, bg_tree) 
			})
		}
	}
	node
}

run_coreAlg_though_tree = function(tree, obj.in){

	
	root_proportion = data.frame(TME=rep(1, ncol(obj.in$mix))) 
	rownames(root_proportion) = colnames(obj.in$mix)
	tree$relative_proportion  = root_proportion
	
	tree = run_coreAlg_though_tree_recursive(tree, obj.in, tree)
	
	tree
}


#' Select the proportions for one specific tree form the table of proportionsd
#'
#' @param node A node from the tree
#' @param name A char
#' @return A tree
pick_proportion_for_specific_tree = function(node, name){

	if(length(node$relative_proportion)>0)
		node$relative_proportion = node$relative_proportion[name,]

	if(length(node$children)>0)
		node$children = lapply(node$children, pick_proportion_for_specific_tree, name)
	
	node
}

#' Divide one tree with probability tables in many trees with probability reals
#'
#' @param node A node from the tree
#' @return A list of trees
divide_trees_proportion_across_many_trees = function(tree){
	my_tree_names = rownames(tree$children[[1]]$relative_proportion)
	trees = lapply(my_tree_names, function(m)  pick_proportion_for_specific_tree(tree, m))
	names(trees) = my_tree_names
	trees
}
