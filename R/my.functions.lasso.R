# node_size = gr[3, ]; max_size=100;
getSizeWeight <- function(node_size, type, max_size) {
	node_size[which(node_size > max_size)] = max_size;
	wt = c();
	if (type == 1) {
		wt = node_size^0.8;
	} else if (type == 2) {
		wt = node_size^0.5;
	} else if (type == 3) {  # seldom work
		wt = node_size^0.5;
		wt[which(node_size > 100)] = 50;
	} else if (type == 4) {
		node_size[which(node_size < 2)] = 2;  # leave nodes
		wt = log(node_size, 2);
	} else if (type == 5) {  # seldom work
		node_size[which(node_size < 2)] = 2;  # leave nodes
		wt = (log(node_size, 2))^0.5;
	} else {
		wt = node_size;
		index_large = which(node_size > max_size);
		wt[index_large] = max_size;	
		a = 10 / max(wt);
		wt = wt * a;
		wt = wt^2;
	}
	return(wt);
}

# arg_response=response; arg_feature=feature; arg_gr=gr; arg_expected=0.01; arg_covar_cnt=0;   ## arg_covar_cnt=20;
checkGroups <- function(arg_response, arg_feature, arg_gr, arg_expected, arg_covar_cnt) {
	cnt_samples = length(arg_response);
	cnt_features = nrow(arg_feature);
	
	cnt_expected = 1; # absolute count
	if (arg_expected < 1)  { # percentage
		cnt_expected = (cnt_features - arg_covar_cnt )* arg_expected;
		if (cnt_expected < 5) {
			cnt_expected = 5;
		}
	} else {
		cnt_expected = arg_expected;
	}
	cnt_expected = round(cnt_expected + arg_covar_cnt, 0);

	par = seq(0.1, 0.99, 0.01); # lambda
	opts=list();
	opts[['init']]=2;        # starting from a zero point
	opts[['tFlag']]=5;       # run .maxIter iterations
	opts[['maxIter']]=100;   # maximum number of iterations
	opts[['rFlag']]=1;       # use ratio
	opts[['nFlag']]=0;		# do not normalize
	opts[['ind']] = arg_gr[c(1,2,4), ];

	p_lb = 1;
	diff = cnt_features;
	beta = c();
	for (p in 1:length(par)) {
		lasso.output = treeLinearLasso(arg_feature, arg_response, par[p], opts);
		w=abs(lasso.output[[1]]);
		beta = cbind(beta, w);  # features on row and iterations on column
		cnt_selected_features = length(which(w!=0));
		d = abs(cnt_selected_features - cnt_expected);
		if (d < diff & cnt_selected_features >= cnt_expected) {
			p_lb = p;
			diff = d;
		} else if (cnt_selected_features < cnt_expected) {
			break;
		}
	}
	beta_p = abs(beta[, p_lb:ncol(beta)]);
	beta_max = apply(beta_p, 1, max);
	selected_features = which(beta_max != 0);
	
#	gr = arg_gr(:, 2:end);
	gr = arg_gr;
	beta_gr = c();
	for (g in 1:ncol(gr)) {
		beta_gr = c(beta_gr, sum(beta_max[gr[1, g]:gr[2, g]]));
	}
	gr = rbind(gr, beta_gr);
	
	selected_group_size = gr[3, which(gr[5, ] > 0)];
	group_size_ub = -100;
	if (length(selected_group_size) > 0) {
		group_size_ub = max(selected_group_size);
	} ## !!!!!!!!!!! Add check larger vs. small groups !!!!!!!
	
	cat(paste('group size upper bound ... ', group_size_ub, ', (', par[p_lb], ', ', length(selected_features), ')\n', sep=' '));

	return(group_size_ub);
}

# arg_phenotype=response; arg_genotype=feature; arg_opts=opts; arg_expected=selected_cnt_group; arg_covar_cnt=0; tree=F; ## arg_covar_cnt=20;
determineLambdaRange <- function(arg_phenotype, arg_genotype, arg_opts, arg_expected, arg_covar_cnt, tree=F) {
	cnt_samples = length(arg_phenotype);
	cnt_features = ncol(arg_genotype);

	cnt_expected = 1; # absolute count
	if (arg_expected < 1)  { # percentage
		cnt_expected = (cnt_features - arg_covar_cnt ) * arg_expected;
		if (cnt_expected < 5) {
			cnt_expected = 5;
		}
	} else {
		cnt_expected = arg_expected;
	}
	cnt_expected = cnt_expected + arg_covar_cnt;
	cnt_expected = round(cnt_expected, 0);
	
	# determine the lower bound of lambda
	par = seq(0.01, 0.99, 0.01); # lambda
	cat(paste('determine the lower bound of lambda ... ', arg_expected, '(', cnt_expected, ')\n', sep=' '));
	lambda_lb = binary.search(par, arg_genotype, arg_phenotype, arg_opts, cnt_expected, tree);
	cat(paste('lambda:', lambda_lb, '\n', sep=' '));
	
	# customize range of lambda
	par = seq(lambda_lb, 0.99, 0.1);	
	if(lambda_lb > 0.5) {
		par = seq(lambda_lb, 0.99, 0.05);	
	}
	
	return(par);
}

binary.search <- function(arg_par, arg_genotype, arg_phenotype, arg_opts, arg_expected, tree=F){  
	cnt_feature = ncol(arg_genotype);
	left = 1;
	right = length(arg_par);
	found = 0;
	while(right - left > 1 & found == 0){
		middle = floor((left + right) / 2)
		lasso.output = c()
		if(tree) {
			lasso.output = treeLinearLasso(arg_genotype, arg_phenotype, arg_par[middle], arg_opts);
		} else {
			lasso.output = linearLasso(arg_genotype, arg_phenotype, arg_par[middle], arg_opts);
		}
		w=lasso.output[[1]];	
		cnt_selected_features = length(which(w != 0));
#		cat(paste('[', left, ',', right, '] ...', middle, '-->', cnt_selected_features, 'vs', arg_expected, '(', cnt_selected_features - arg_expected, ')', arg_par[middle], '\n', sep=' '));
		if (cnt_selected_features < arg_expected) {
			right = middle - 1;
		} else if ((cnt_selected_features - arg_expected) < 0.1 * arg_expected) {
			found = middle;
		} else {
			left = middle;
		}
	}
	if(found == 0) {
		found = left;
	}
	lambda_lb = arg_par[found];

    return(lambda_lb)
}

# arg_phenotype=this_response; arg_genotype=this_feature; arg_opts=opts; arg_par=par_tree; tree=F
stabilitySelect_tree <- function(arg_phenotype, arg_genotype, arg_opts, arg_par, tree=F) {
	cnt_samples = length(arg_phenotype);
	cnt_features = ncol(arg_genotype);

	par = c(0.4, 0.5, 0.6, 0.7, 0.8, 0.9); # lambda
	if (length(arg_par) > 0) {
		par = arg_par;
	}
	opts=arg_opts;	
	gr = opts[['ind']];
	
	beta = c();
	for (p in 1:length(par)) {
		lasso.output = c()
		if(tree) {
			lasso.output = treeLinearLasso(arg_genotype, arg_phenotype, par[p], opts);
		} else {
			lasso.output = linearLasso(arg_genotype, arg_phenotype, par[p], opts);
		}
		w = lasso.output[[1]];
		beta = cbind(beta, w);  # features on row and iterations on column
	}
	
	comb_beta = c();
	for (b in 1:ncol(beta)) {
		this_beta = abs(beta[, b]);
		gr_beta = c();
		for (g in 1:ncol(gr)) {
			gr_beta = c(gr_beta, sum(this_beta[gr[1, g]:gr[2, g]])); 
		}
		comb_beta = cbind(comb_beta, gr_beta);  # sum or mean beta, groups on row and iterations on column
	}
	comb_beta.max = apply(comb_beta, 1, max);
	max_beta = t(rbind(gr, comb_beta.max));
	
	return(max_beta);
}


library(TreeGuidedRegression)
library(tictoc);
library(parallel);
# flag = 1; flag_string = 'exp'; output_folder = 'tree.lasso/';  ## for gtex
# flag = 0; flag_string = 'response'; output_folder = './';  ## for simulation
# output_file_name = 'tree-lasso.func5'; 
# weight_flag = 0; 
# max_level = 8;
# skip = 1;  ## 0-not skip; 1-skip
# maxIter = 100;
run.tree.lasso <- function(name, flag, flag_string, path='', output_file_name='tree-lasso.funcr', weight_flag=0, max_level=8, skip=0, maxIter=100, time.limit=0) {
	set.seed(0);
	if(time.limit > 0) {
		setTimeLimit(cpu=time.limit);
	} else {
		setTimeLimit(cpu=Inf);
	}
	tic();
	one_response_file = name;
	one_feature_file = gsub(flag_string, 'geno.norm', one_response_file);
	one_group_file = gsub(flag_string, 'group', one_response_file);
	one_zip_file = one_feature_file;
	one_feature_file <- paste(path, one_feature_file, sep='');
	one_group_file <- paste(path, one_group_file, sep='');
	one_zip_file <- paste(path, one_zip_file, sep='');
	one_response_file <- paste(path, one_response_file, sep='');
	one_lasso_file = gsub(flag_string, output_file_name, one_response_file);
	one_lasso_file = gsub('.gz', '', one_lasso_file);
	
	if(flag == 1) {
		one_response_file <- paste(path, 'gene.exp.orig.maf05.adj/', name, '.exp.csv.gz', sep='');		
		one_feature_file <- paste(path, 'gene.feature.orig.maf05/', name, '.feature.csv.gz', sep='');
		one_group_file <- paste(path, 'gene.group.maf05/', name, '.group.csv.gz', sep='');
		one_zip_file <- paste(path, 'gene.geno.norm.maf05/', name, '.geno.norm.csv.gz', sep='');
		one_lasso_file <- paste(path, 'tree.lasso.maf05/', name, '.lasso.csv', sep='');
	}
	
	cat(paste(one_zip_file, '...\n', sep=' '));
	if (file.exists(one_zip_file) & file.exists(one_group_file)) {
		cat(paste('start processing ...', one_response_file, '\n', sep=' ')); flush.console();
	
		if (file.exists(one_lasso_file) & skip == 1) {
			cat(paste('existing ', one_lasso_file, '\n', sep=' '));  flush.console();
		} else {
			response = read.table(one_response_file, sep=',', header=F, stringsAsFactors=F); 
			feature = read.table(one_feature_file, sep=',', header=F, stringsAsFactors=F); 
			if(flag == 1) {
				feature = read.table(one_zip_file, sep=',', header=F, stringsAsFactors=F); 
			}
			group_all = read.table(one_group_file, sep=',', header=F, stringsAsFactors=F); 
			response = response[, 1];
			feature = as.matrix(feature);
			cnt_sample = length(response);
			cnt_feature = ncol(feature);		

			n = length(which(group_all[, 6] == 1));
			n = round(n * 0.05, 0) + 1;  # 5% of groups
			s = sum(group_all[which(group_all[, 6] == 1), 3]);
			s = round((cnt_feature - s) * 0.01, 0) + 1; # 1% of singles
			m = round(max(group_all[, 3]) * 0.01, 0);  # 1% largest group members
			selected_cnt_leave = n + s;
			selected_cnt_group = n * 3 + s;
				
			max_size = max(group_all[, 3]);

			# type_size = c();
			# gr.orig = t(group_all[, 1:3]);
			# for (type in 1:4) {
				# gr = rbind(gr.orig, getSizeWeight(gr[3, ], type, max_size));
				# group_size_ub = checkGroups(response, feature, gr, 0.01, 0);
				# type_size = rbind(type_size, c(type, group_size_ub));
			# }
			# m = which.min(type_size[, 2]);
			# type = type_size[m, 1];
			# group_size_ub = type_size[m, 2];
			# cat('type ', type, '\n'); flush.console();
				
			# if (group_size_ub < 100) {
				# selected_cnt_leave = s;
				# selected_cnt_group = s + 5;
			# } else {
				# selected_cnt_group = n * 3 + s + m;
			# }
			selected_cnt_group = 0.01;
			type = 2;

			# different weight schemes
			group = group_all[which(group_all[, 6] <= max_level), ];
			node_size = group[, 3];
			maf = group[, 4];
			maf = (maf*(1-maf))/0.25;
			func = group[, 5];
			func = 1 - func;
			group_active = list();
			g = 1;
			group_active[[g]] = t(cbind(group[, 1:2], getSizeWeight(node_size, type, max_size))); g = g + 1;  
			group_active[[g]] = t(cbind(group[, 1:2], getSizeWeight(node_size, type, max_size) + maf)); g = g + 1;  
			group_active[[g]] = t(cbind(group[, 1:2], getSizeWeight(node_size, type, max_size) + maf + func*5)); g = g + 1;  
					
			if (weight_flag == 0) {
				group_active = group_active[1:3];
			} else if (weight_flag == 1) {  ## not using
				group_active = group_active[1:7];
			} else if (weight_flag == 1.5) {  ## not using
				group_active = group_active[c(1:3, 8:11)];
			} else if (weight_flag == 10) {  ## not using
				group_active = group_active[1:2];
			}

			opts=list();
			opts[['init']]=2;        # starting from a zero point
			opts[['tFlag']]=5;       # run .maxIter iterations
			opts[['maxIter']]=maxIter;   # maximum number of iterations
			opts[['rFlag']]=1;       # use ratio
			opts[['nFlag']]=0;		# do not normalize
				
			beta_all = group[, 1:6];
			iteration = 3;
#			for (a in 1:length(group_active)) {
			for (a in 1) {
				gr = group_active[[a]]; 
				gr_tree = gr;
				opts[['ind']] = gr_tree;
				tic();
					par_tree = determineLambdaRange(response, feature, opts, selected_cnt_group, 0, tree=T);
				toc();
				
				gr_leave = gr[, 1:cnt_feature];
				opts[['ind']] = gr_leave;
				tic();
					par_leave = determineLambdaRange(response, feature, opts, selected_cnt_leave, 0, tree=F);
				toc();
				
				for (t in 1:iteration) {
					random_index = 1:cnt_sample;
					if (t > 1) {
						random_size = round(cnt_sample * 0.99, 0);
						random_index = sample(1:cnt_sample, random_size);
					}
					this_response = response[random_index];
					this_feature = feature[random_index, ];
					# full tree
					opts[['ind']] = gr_tree;
					tic();
						max_beta = stabilitySelect_tree(this_response, this_feature, opts, par_tree, tree=T);
					toc();
					beta_all = cbind(beta_all, max_beta[, 4]);
					# leaves only - regular lasso 
					opts[['ind']] = gr_leave;
					tic();
						max_beta = stabilitySelect_tree(this_response, this_feature, opts, par_leave, tree=F);
					toc();
					beta_all = cbind(beta_all, c(max_beta[, 4], rep(0, ncol(gr_tree)-cnt_feature)));
				}
			}		
			colnames(beta_all) <- '';

			# write output
			write.table(round(beta_all, 5), one_lasso_file, sep=',', row.names=F, col.names=F, quote=F);			
		}
		toc();
		cat(paste('finished lasso ...', one_response_file, '\n', sep=' ')); flush.console();
	}
}
## test
# setwd('~/my.scratch/projects/GTEx/simulation-2018/tree.lasso/1_causal/'); 
# response.files = list.files('./', pattern='*.response.csv.gz');
# run.tree.lasso(response.files[1], flag=0, flag_string='response', output_folder='./');
# dummy = mclapply(response.files[1:5], run.tree.lasso, flag=0, flag_string='response', output_folder='./');
# setwd('~/my.scratch/projects/GTEx/brain.hippo')
# response.files = list.files('gene.exp.maf05.adj/', pattern='*.exp.csv.gz');
# run.tree.lasso(response.files[1], flag=1, flag_string='exp', output_folder='tree.lasso/');

# input.folder.path='./'; output.folder.path='./';
# input.list <- create.processing.list.gtex(input.folder.path, output.folder.path);
# run.tree.lasso(name=input.list[1], flag=1, flag_string='exp', output_file_name='tree.lasso/');


