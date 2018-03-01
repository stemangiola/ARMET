#' Add info to tree
add_info_to_tree = function(node, ct, label, value, append = F){
	
	if(is.na(ct)) stop("ARMET: cell type not provided. In add_info_to_tree")
	
	else if(node$name==ct) {
		if(append) 
			node[[label]] =
				node[[label]] %>% 
				dplyr::mutate_if(is.factor, as.character) %>%
				dplyr::bind_rows(
					value %>% 
						dplyr::mutate_if(is.factor, as.character)
				) %>%
				dplyr::mutate_if(is.character, as.factor)
		else node[[label]] = value
		node
	}
	else if(length(node$children)>0){
		node$children = lapply(node$children, add_info_to_tree, ct, label, value, append)
		node
	}
	else node
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
			node[[label]],
			do.call("c", lapply(node$children, get_node_label_recursive, last_level, label))
		)
		else if(last_level==1) do.call("c", lapply(node$children, get_node_label_recursive, last_level, label))
	}
	else if(last_level>=0) node[[label]]
}

get_last_existing_leaves_with_annotation = function(node, what = "relative_proportion"){
	
	# If did not reach the end
	if(length(node$children)>0){
		
		do.call("rbind", lapply(node$children,function(no)	
			get_last_existing_leaves_with_annotation(no, what)
		)) %>%
		{ 
			if(nrow(.)==0) 
				(.) %>% 
				dplyr::mutate_if(is.factor, as.character) %>%
				dplyr::bind_rows(
					node[[what]] %>%
						dplyr::mutate_if(is.factor, as.character)
				) %>%
				dplyr::mutate_if(is.character, as.factor)
			else 
				.
		}
	}	
	
	# If last leaf
	else node[[what]]
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
	else return(node[[label]])
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

#' Get the right level of the tree and label the leaf
get_child_and_group_background = function(node, ct){
	
	# get childred
	ct_main = get_leave_names(node_from_name(node, ct), recursive=F)
	
	# get backgrounbd
	get_background = function(ct){
		ct_childrens = get_leave_names(node_from_name(node, ct), recursive=F)
		ct_ancestor = rev(get_hierarchy(node, ct))[2]
		if(length(ct_ancestor)>0) {
			ancestor_leaves = get_background(ct_ancestor)
			c(ct_childrens, ancestor_leaves[ancestor_leaves != ct])
		}
	}
	ct_background = get_background(ct)
	ct_background = ct_background[!ct_background%in%ct_main]
	
	list(
		ct_main = ct_main,
		ct_background = ct_background
	)
	
}

#' Create a map where every leaf is mapped with the right ancestor
get_map_foreground_background = function(node, ct){
	
	fg_bg = get_child_and_group_background(node, ct)
	fg_bg = fg_bg[lapply(fg_bg,length)>0]
	
	add_children_to_line = function(df){
		no = node_from_name(node, df$ancestor)
		children = if(length(no$children)>0) get_node_label_recursive(no) else no$name
		do.call("rbind", lapply(1:length(children), function(dummy) df)) %>% 
			dplyr::mutate(ct = children)
	}
	
	reshape::melt(fg_bg) %>% 
		tibble::as_tibble() %>%
		dplyr::rename(ancestor=value, variable=L1) %>%
		dplyr::mutate(variable = gsub("ct_", "", variable)) %>%
		dplyr::group_by(ancestor) %>%
		dplyr::do(add_children_to_line(.)) %>%
		dplyr::ungroup() %>%
		dplyr::mutate(ancestor = as.character(ancestor))
}

#' Put absolute proportions on tree
add_absolute_proportions_to_tree = function(node, p_ancestor =  rep(1, nrow( node$relative_proportion )) ){
	
	node$relative_proportion = 
		node$relative_proportion %>% 
		dplyr::mutate(
			absolute_proportion = 
				relative_proportion * 
				p_ancestor
		) 
	
	if(length(node$children)>0){
		node$children = 
			lapply(
				node$children, 
				add_absolute_proportions_to_tree, 
				node$relative_proportion %>% dplyr::pull(absolute_proportion)
			)
	}
	node
}

#' Put proportions on tree
add_data_to_tree_from_table = function(tree, proportion, label, append  ){
	
	# Update tree so it does not depend
	for(
		ll in 
		proportion %>% 
		dplyr::pull(ct) %>% 
		levels()
	){
		
		tree = 
			add_info_to_tree(
				tree, 
				ll, 
				label, 
				proportion %>% dplyr::filter(ct == ll),
				append = append
			)
		
	}
	
	tree
}


#' Run ARMET core algorithm recursively on the tree 
#'
#' @param node A node from the tree
#' @param obj.in A object with the input
#' @param bg_tree The background tree that tracks the evolution of the algoritm
#' @return A probability array
run_coreAlg_though_tree_recursive = function(node, obj.in, bg_tree, log.ARMET){
	
	# library(doFuture)
	# doFuture::registerDoFuture()
	# future::plan(future::multiprocess)
	
	
	if(length(node$children)>0){
		
		# Initialize pipe
		`%>%` <- magrittr::`%>%`
		
		obj.out = ARMET_tc_coreAlg(obj.in, node)
		node = obj.out$node
		
		write(node$name,file=log.ARMET,append=TRUE)
		
		# Update results for the background tree
		bg_tree = add_data_to_tree_from_table(bg_tree, obj.out$proportion, label = "relative_proportion", append = T)
		
		# Set up background trees
		bg_tree = add_absolute_proportions_to_tree(bg_tree)
		
		obj.in$my_tree = bg_tree
		
		if(obj.in$multithread & !obj.in$do_debug){
			n_cores = length(node$children)
			
			cl <- parallel::makeCluster(n_cores)
			parallel::clusterExport(cl, c("obj.in", "bg_tree", "%>%", "log.ARMET"), environment())
			doParallel::registerDoParallel(cl)
		}
		
		`%my_do%` = ifelse(obj.in$multithread & !obj.in$do_debug, `%dopar%`, `%do%`)
		verbose = ifelse(obj.in$multithread & !obj.in$do_debug, F, T)
		
		#################################
		## Debug
		#################################
		if(obj.in$do_debug)
		node$children = lapply (node$children, function(cc)  {
			
			run_coreAlg_though_tree_recursive(cc, obj.in, bg_tree, log.ARMET)
		})
		else 
			node$children = foreach::foreach(cc = node$children, .verbose = verbose) %my_do% {
	
				run_coreAlg_though_tree_recursive(cc, obj.in, bg_tree, log.ARMET)
			}
		
		if(obj.in$multithread & !obj.in$do_debug)	parallel::stopCluster(cl)
		
		node
	}
	
	node
}

run_coreAlg_though_tree = function(node, obj.in){
	
	future::plan(future::multiprocess)
	
	
	writeLines("ARMET: Starting deconvolution")
	
	log.ARMET = sprintf("%s%s", tempfile(), Sys.getpid())
	if (file.exists(log.ARMET)) file.remove(log.ARMET)
	file.create(log.ARMET)
	
	exec_hide_std_out = function(node, obj.in, log.ARMET){
		
		if(obj.in$multithread & !obj.in$do_debug){
			cl <- parallel:::makeCluster(2)
			parallel:::clusterExport(cl, c("obj.in", "node", "%>%", "log.ARMET"), environment())
			doParallel:::registerDoParallel(cl)
		}
		
		`%my_do%` = ifelse(obj.in$multithread & !obj.in$do_debug, `%dopar%`, `%do%`)
		verbose = ifelse(obj.in$multithread & !obj.in$do_debug, F, T)
		
		node.filled = foreach:::foreach(dummy = 1, .verbose = verbose) %my_do% {
		
			run_coreAlg_though_tree_recursive(node, obj.in, node, log.ARMET)
		}
		
		if(obj.in$multithread & !obj.in$do_debug)	parallel:::stopCluster(cl)
		
		return(node.filled[[1]])
	}
	
	if(!obj.in$do_debug)
	{
		node.filled = future::future(		exec_hide_std_out(node, obj.in, log.ARMET) 		)
		
		
		
		log.array = c()
		done = F
		
		while(1){
			
			Sys.sleep(2)
			
			temp = readLines(log.ARMET)
			
			new = setdiff(temp, log.array)
			
			if(length(new)>0) {
				writeLines(sprintf("ARMET: %s deconvolution completed", new))
				log.array = c(log.array, new)
			}
			
			
			if(all( get_leave_label(node, last_level = -1) %in% log.array)) break
		}
		
		file.remove(log.ARMET)
		
		return(future::value(node.filled))
		
	} else {
		return( run_coreAlg_though_tree_recursive(node, obj.in, node, log.ARMET) )
	}
	
}


# drop_node_from_tree = function(node, ct){
# 	
# 	if( !node$name%in%ct){
# 		node$children = lapply(node$children, drop_node_from_tree, ct)
# 		node
# 	}
# }

#' Drop a node from the tree
drop_node_from_tree = function(node, ct){
	if( !node$name%in%ct) {
		node$children = Filter(Negate(is.null), 
													 lapply(node$children, drop_node_from_tree, ct)
		)
		node
	}
}

format_tree = function(tree, mix, ct_to_omit){
	
	
	
	# Drop unwanted nodes
	tree = drop_node_from_tree(tree, ct_to_omit)
	
	
	# Add empty proportion table to all nodes
	for(ct in get_leave_label(tree)) {
		tree = 
			add_info_to_tree(
				tree, 
				ct = ct, 
				label = "relative_proportion", 
				value = tibble::tibble(
					sample = factor(), 
					ct = factor(), 
					relative_proportion=double(), 
					absolute_proportion=double()
				)
			)
	}
	
	# Add the first proportions to tree
	tree = 
		add_info_to_tree(
			tree,
			ct = tree$name,
			label = "relative_proportion", 
			value = (
				mix %>% 
					dplyr::distinct(sample) %>% 
					dplyr::mutate(
						ct = factor(tree$name), 
						relative_proportion = 1,
						absolute_proportion = 1
					)
			)
		)
	
	tree
}


#' Select the proportions for one specific tree form the table of proportionsd
#'
#' @param node A node from the tree
#' @param name A char
#' @return A tree
pick_proportion_for_specific_tree = function(node, name){
	
	# If the children has proportion information
	node$relative_proportion = 
		node$relative_proportion %>%
		dplyr::filter(sample == name)
	
	# Recursive selection
	if(length(node$children)>0)
		node$children = 
			lapply(
				node$children, 
				pick_proportion_for_specific_tree, 
				name
			)
	node
}

#' Divide one tree with probability tables in many trees with probability reals
#'
#' @param node A node from the tree
#' @return A list of trees
divide_trees_proportion_across_many_trees = function(node){
	
	my_tree_names = 
		node$relative_proportion %>% 
		dplyr::pull(sample) %>% 
		levels()
	
	trees = 
		lapply(
			my_tree_names, 
			function(m)  pick_proportion_for_specific_tree(node, m)
		)
	
	names(trees) = my_tree_names
	
	trees
}

get_tree_hypoth_test = function(tree_out, tree_in){
	
	for( ct in tree_out$Get('name') ){
		
		if(ct %in% c("dendritic","TME")) next
		ti = node_from_name(tree_in, ct)
		if(length(ti[["stats"]])==0) next
		s = ti[["stats"]] %>% dplyr::mutate_if(is.factor, as.character)
		
		data.tree::FindNode(tree_out, ct)$Set(estimate_extrinsic =      round(s$mcap_human_readable, 2), filterFun = function(x) x$name == ct)
		data.tree::FindNode(tree_out, ct)$Set(std_error_extrinsic =     round( s$ecap,2), filterFun = function(x) x$name == ct)
		data.tree::FindNode(tree_out, ct)$Set(direction_extrinsic =     s$direction, filterFun = function(x) x$name == ct)
		data.tree::FindNode(tree_out, ct)$Set(pvalue_extrinsic =        s$pcap_human_readable, filterFun = function(x) x$name == ct)
		data.tree::FindNode(tree_out, ct)$Set(significance_extrinsic =  s$symbol_pcap, filterFun = function(x) x$name == ct)
		
	}
	
	tree_out
	
}

#' Get information from tree recursively for every node
get_distributions_from_tree = function(node){
	
	
	
	if(length(node$children)>0){

			node[["summary"]] %>%
			dplyr::bind_rows(
				do.call(dplyr::bind_rows, lapply(node$children, get_distributions_from_tree))
			)

	}
	else node[["summary"]]
}

get_gene_distributons_recursive = function(node, ref_tbl){

	if(length(node$children)>0){
		
		doParallel::registerDoParallel(length(node$children))
		
		node$children = foreach(nn = node$children) %do% {
				get_gene_distributons_recursive(nn, ref_tbl)
		}
		
		node$summary =
				do.call(dplyr::bind_rows, lapply(node$children, function(n) n$summary)) %>% 
					dplyr::group_by(gene) %>%
					# Safer alternative to calculation if NA are present
					dplyr::summarise(
						log_mean =  ifelse(
							all(is.na(log_sd)), 
							NA, 
							mean(log_mean, na.rm = T)  
						),
						log_sd = ifelse(
							all(is.na(log_sd)), 
							NA, 
							mean(log_sd, na.rm = T)  # sqrt( sum(log_sd^2, na.rm = T))  
						),
						ct = node$name
					) %>%
					dplyr::ungroup() %>%
					dplyr::mutate(original = F)
	}
	
	else {
		node$summary = 
			ref_tbl %>% 
				dplyr::filter(ct == node$name) %>%
				dplyr::mutate(ct = as.character(ct)) %>%
				dplyr::mutate(original = T)
	}
	
	node
}



