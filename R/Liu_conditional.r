require(MASS)

parse.lasso <- function(input, snp.names, group, lasso.file.path, top.file.path, weighted='maf') {
	cat(paste('analyzing lasso output', lasso.file.path, '...\n', sep=' ')); flush.console();
	temp <- colnames(input)[2];
	if(substr(temp, 0, 1) == 'X') {
		snp.names <- gsub(':', '.', snp.names);
		snp.names <- paste('X', snp.names, sep='');
		snp.names <- gsub('XX', 'X', snp.names);
	}
	colnames(input)[1] <- 'resp';
	
	lasso.all <- read.table(lasso.file.path, sep=',', header=F, stringsAsFactors=F);
	col.end <- ncol(lasso.all);	
	if(weighted == 'both') {
		col.end <- 6 + 6*3;
	} else if (weighted == 'maf') {
		col.end <- 6 + 6*2;
	} else if (weighted == 'function') {
		col.end <- 6 + 6*2;
	} else {
		col.end <- 6 + 6;
	}

	max.level <- 8;
	qnt <- 0.95;
	lasso.all <- lasso.all[which(lasso.all[, 6] <= max.level), ];

	tops <- c();
	for(l in seq(7, col.end, 2)) {
		lasso <- lasso.all[, c(1:6, l, l+1)];
		colnames(lasso) <- c('start', 'end', 'size', 'maf', 'func', 'level', 'beta_group', 'beta_single');	
		idx <- which(lasso$size > 1)
		lasso[idx, 'beta_group'] <- lasso[idx, 'beta_group']/lasso[idx, 'size'];
		lasso[idx, 'beta_single'] <- lasso[idx, 'beta_single']/lasso[idx, 'size'];
		lasso$beta_group <- lasso$beta_group * 1/(1 + group$weight);
		lasso$beta_single <- lasso$beta_single * 1/(1 + group$weight);
		lasso$selected <- 0;
		## 1: select leaf SNPs based on beta_single 
			threshold <- quantile(lasso$beta_single, qnt, na.rm=T);
			lasso[which(lasso$beta_single > 0 & lasso$beta_single >= threshold), 'selected'] <- 1;  
		## 2: select node SNPs based on beta_single
#			threshold <- quantile(lasso[which(lasso$size > 1), 'beta_single'], qnt, na.rm=T);
#			lasso[which(lasso$size > 1 & lasso$beta_single > 0 & lasso$beta_single > threshold), 'selected'] <- 3; 
			leaf.pos <- lasso[which(lasso$selected == 1), 'start'];
			leaf.group.index <- unique(unlist(lapply(leaf.pos, function(x) { which(lasso$start <= x & lasso$end >= x & lasso$size > 1) })));
			lasso[leaf.group.index, 'selected'] <- 2;
		## 3: select node SNPs based on beta_group
			threshold <- quantile(lasso[which(lasso$size > 1), 'beta_group'], qnt, na.rm=T);
			lasso[which(lasso$size > 1 & lasso$beta_group > 0 & lasso$beta_group > threshold), 'selected'] <- 3; 
		## 4 & 5: select leaf SNPs based on beta_group 
			threshold <- quantile(lasso[which(lasso$size == 1), 'beta_group'], qnt, na.rm=T);
			lasso[which(lasso$size == 1 & lasso$beta_group > 0 & lasso$beta_group > threshold & lasso$selected == 0), 'selected'] <- 4;  
			for(s in which(lasso$selected == 3)) {
				this.selected <- lasso[s, ];
				leaf.inside <- lasso[which(lasso$size == 1 & lasso$start >= this.selected$start & lasso$start <= this.selected$end), ];
				threshold <- quantile(leaf.inside$beta_group, qnt, na.rm=T);
				leaf.inside.selected <- leaf.inside[which(leaf.inside$beta_group > 0 & leaf.inside$beta_group > threshold), ];
				lasso[which(lasso$size == 1 & lasso$start %in% leaf.inside.selected$start & lasso$selected == 0), 'selected'] <- 5;
				lasso[which(lasso$size == 1 & !lasso$start %in% leaf.inside.selected$start & lasso$selected == 4), 'selected'] <- 0;
			}
		top.lasso <- lasso[which(lasso$selected > 0), ]
		tops <- rbind(tops, top.lasso);		
	}
	top.agg <- aggregate(selected ~ start + end + level, data=tops, 'min');
	colnames(top.agg) <- c('start', 'end', 'level', 'selected');
	top.agg$size <- top.agg$end - top.agg$start + 1;
	top.agg <- top.agg[order(top.agg$start), ]
	top <- top.agg;	
	cat(paste('   top table:', nrow(top), 'rows.\n', sep=' ')); flush.console();
	
	write.table(top, top.file.path, sep='\t', row.names=F, quote=F);
	
	return(top);
}
	
select.within.groups <- function(top, input, snp.names, group, within.file.path) {
	# stepwise selection within nodes
	top.within.group <- c(); ## index in "input" order
	check.group.index <- which(top$size > 1); ## index in "top" order
	check.leaf.index <- c(); ## index in "top" order
	for(t in check.group.index) {
		this.group.start <- top[t, 'start'];
		this.group.end <- top[t, 'end'];
		within.index <- which(top$size == 1 & top$start >= this.group.start & top$start <= this.group.end); ## index in "top" order
		check.leaf.index <- c(check.leaf.index, within.index); ## leaves inside a node, in "top" order
		sig.index.1 <- stepwise.leaf(input, index=top[within.index, 'start'], snp.names, weight=group$weight, type='conditional')
		sig.index.2 <- c();
		sig.index.2 <- stepwise.leaf(input, index=this.group.start:this.group.end, snp.names, weight=group$weight, type='conditional');  
		top.within.group <- unique(c(top.within.group, sig.index.1, sig.index.2));  ## index in "input" order
	}
	check.index <- unique(c(check.group.index, check.group.index)); ## nodes + leaves, index in "top" order
	uncheck.index <- which(!1:nrow(top) %in% check.index);  ## some leaves are not part of a group, index in "top" order
	top.within.group <- c(top.within.group, top[uncheck.index, 'start']);  ## index in "input" order
	top.within.group <- unique(top.within.group);
	top.within.group <- top.within.group[order(top.within.group)];
	
	# add linked SNPs to within-group selections
	top.within.linked.perfect <- top.within.group;
	top.within.linked <- top.within.group;
	for(n in top.within.group) {
		a <- input[, n + 1];
		b <- input[, -1];
		cor <- cor(a, b, use='complete.obs')^2; ## !!!!!!!!!!-------------
		linked.perfect <- which(cor > 0.9999);  ##  !!!!!!!!!!-------------
		linked <- which(cor > 0.8 & cor < 0.9999);  ##  !!!!!!!!!!-------------
		top.within.linked.perfect <- c(top.within.linked.perfect, linked.perfect);
		top.within.linked <- c(top.within.linked, linked);
	}
	top.within.linked.perfect <- unique(top.within.linked.perfect);	
	top.within.linked.perfect <- top.within.linked.perfect[order(top.within.linked.perfect)];
	top.within.linked <- unique(top.within.linked);	
	top.within.linked <- top.within.linked[order(top.within.linked)];
	cat(paste('   top.within.group:', length(top.within.group), '-- top.within.linked:', length(top.within.linked), '-- top.within.linked.perfect:', length(top.within.linked.perfect), '\n', sep=' ')); flush.console();

	selected.within <- list(top.within.group=top.within.group, top.within.linked=top.within.linked, top.within.linked.perfect=top.within.linked.perfect);
	save(selected.within, file=within.file.path);	
	return(selected.within)
}
	
select.cross.groups <- function(selected.within, input, snp.names, group, cross.file.path) {
	top.within.group <- selected.within[['top.within.group']];
	top.within.linked <- selected.within[['top.within.linked']];
	top.within.linked.perfect <- selected.within[['top.within.linked.perfect']];

	# conditional test across groups
	top.cross.group <- stepwise.leaf(input, top.within.group, snp.names, weight=group$weight, type='conditional');
#	top.cross.group <- stepwise.leaf(input, top.within.linked, snp.names, weight=group$weight, type='conditional');
	cat(paste('   top.cross.group:', length(top.cross.group), '\n', sep=' ')); flush.console();

	# add linked SNPs to cross-group selections
	top.cross.linked.strong <- c(top.cross.group);
	top.cross.linked.weak <- c(top.cross.group);
	for(n in top.cross.group) {
		a <- input[, n + 1]
		b <- input[, -1]
		cor <- cor(a, b, use='complete.obs')^2;  ##  !!!!!!!!!!-------------
		linked.strong <- which(cor > 0.99);  ##  !!!!!!!!!!-------------
		linked.weak <- which(cor > 0.95);  ##  !!!!!!!!!!-------------
		top.cross.linked.strong <- c(top.cross.linked.strong, linked.strong);
		top.cross.linked.weak <- c(top.cross.linked.weak, linked.weak);
	}	
	top.cross.linked.strong <- unique(top.cross.linked.strong)
	top.cross.linked.weak <- unique(top.cross.linked.weak);

	top.cross.linked <- top.cross.linked.strong;
	if(length(top.cross.linked) > 0) {
		top.cross.linked <- top.cross.linked[order(top.cross.linked)]
	}
	cat(paste('   top.cross.linked:', length(top.cross.linked), '\n', sep=' ')); flush.console();
		
	final.linked.strong <- unique(c(top.cross.linked.strong, top.within.linked.perfect));
	if(length(final.linked.strong) > 0) {
		final.linked.strong <- final.linked.strong[order(final.linked.strong)]
	}
	final.linked.weak <- unique(c(top.cross.linked.weak, top.within.linked.perfect));
	if(length(final.linked.weak) > 0) {
		final.linked.weak <- final.linked.weak[order(final.linked.weak)]
	}
	final.linked <- final.linked.strong;  ## final.linked.strong is a subset of final.linked.weak (cor>0.98 vs. cor>0.995)
	cat(paste('   final.linked:', length(final.linked), '\n', sep=' ')); flush.console();
	
	selected.cross <- list(top.cross.group=top.cross.group, top.cross.linked=top.cross.linked, top.cross.linked.weak=top.cross.linked.weak, final.linked=final.linked, final.linked.weak=final.linked.weak, top.within.group=top.within.group);
	save(selected.cross, file=cross.file.path);	
	
	return(selected.cross)
}

uni.unorder <- function(resp, ftr, ftr.name='') {
	p.values <- c();
	if(is.vector(ftr)) {
		ftr <- as.data.frame(ftr);
		colnames(ftr) <- ftr.name;
	}
	for(f in 1:ncol(ftr)) {
		temp.data <- data.frame(response=resp, feature=ftr[, f]);
		if(length(unique(temp.data$feature)) == 1) {
			p <- 1;
		} else {
			fit <- lm(response ~ feature, data=temp.data, na.action=na.exclude)
			p <- summary(fit)$coefficients[2, 4];
		}
		p.values <- c(p.values, p);
	}
	result <- data.frame(index=1:ncol(ftr), rs=colnames(ftr), p.value=p.values, stringsAsFactors=F);
	return(result);
}

uni.order <- function(resp, ftr, ftr.name='') {
	result <- uni.unorder(resp, ftr, ftr.name);
	result <- result[order(result$p.value), ]
	return(result);
}

## input: 1st column is the response
## index: start positon of SNPs. So need to add 1 to get the indices in the input matrix
library(MASS)
stepwise.leaf <- function(input, index, snp.names, weight=NULL, type='both') {
	if(length(index) > 0) {
		colnames(input)[1] <- 'resp';
		input.sub <- input[, c(1, index+1)];
		sig.snp.index <- c();
		if(type == 'conditional') {
			cond <- my_conditional(input.sub, weight=weight[index]);
			if(length(cond) > 1) { 
				sig.name <- cond$snp;
				sig.snp.index <- which(snp.names %in% sig.name);
			}
		} else if(type == 'comb') {
			sig.snp.index <- test_comb(input=input, index=index, locus=NULL, cn=2);
		} else if(type == 'my_backward') {
			sig.snp.index <- my_backward(input=input, snp.names=snp.names, index=index);
		} else {
			model.full <- lm(resp ~ ., data=input.sub);
			model.step <- stepAIC(model.full, direction=type, trace=F);
			sig.snp <- rownames(summary(model.step)$coeff);
			sig.snp <- sig.snp[-1];
			sig.snp.index <- which(snp.names %in% sig.snp)
		}
		return(sig.snp.index);
	}
}
## testing
#index=top[within.index, 'start']

my_conditional <- function(data.all, weight=NULL, time=F) {
	colnames(data.all)[1] <- 'resp'; ## 1st column is the response.
	len <- ncol(data.all)
	if(is.null(weight)) {
		weight=rep(1, len-1);
	}
	names(weight) <- colnames(data.all)[-1];  

	kept <- c()
	data1 <- data.all
	time.1 <- as.numeric(Sys.time());
	col_checked <- c();
	flag <- T;
	while(flag) {  
		time.2 <- as.numeric(Sys.time());
		if(time & time.2 - time.1 > 300) {
			cat(paste('  WARNING ! time out, skipping', k, 'th snp and afterwards.\n'));
			break();
		}
		
		if(length(kept) > 0) { ## conditioning
			string <- paste("resp ~ ", paste(kept, collapse="+"), sep="")
			fit <- lm(data=data.all, formula=as.formula(string), na.action=na.exclude)
			data1 <- data.all;
			data1[, 1] <- residuals(fit)
		}

		## Find the best single snp on conditioned response
		data1 <- data1[, which(!colnames(data1) %in% col_checked)];
		len1 <- ncol(data1);
		col1 <- colnames(data1);
		target <- 0
		p_target <- 1
		for(i in 2:len1) {  
			string <- paste("resp ~ ", col1[i], sep="")
			fit <- lm(data=data1, formula=as.formula(string), na.action=na.exclude)
			coeff <- summary(fit)$coefficients
			var_names <- rownames(coeff)[-1];
			p <- coeff[2, 4];
			p_adj <- p * weight[var_names]
			if(!is.na(p_adj) & p_adj < p_target) { 
				p_target <- p_adj
				target <- i
			}
		}
		target_name <- col1[target];
		col_checked <- c(col_checked, target_name);

		cutoff=max(0.001, 0.05/(len1-1))  ## cutoff=0.01
		if(p_target > cutoff) {
			flag <- F;
		} else {
			if(length(kept) == 0) {
				kept <- target_name
			} else {
				kept_test <- c(kept, target_name)
				string <- paste("resp ~ ", paste(kept_test, collapse="+"), sep="")
				fit <- lm(data=data.all, formula=as.formula(string), na.action=na.exclude)
				coeff <- summary(fit)$coefficients[-1, ]
				var_names <- row.names(coeff);
				p <- coeff[ ,4];
				p_adj <- p * weight[var_names]
				sig_names <- var_names[which(p_adj <= 0.05)]
				kept <- sig_names
			}
		}
		
		if(length(col_checked) == (len-1)) {
			flag <- F;
		}
	}           

	if(!exists('kept') || length(kept) == 0 || is.na(kept)) {
		return(0)
	} else {
		fit <- lm(data=data.all, formula = as.formula(paste("resp ~ ", paste(kept, collapse="+"), sep="")), na.action=na.exclude)
		coeff <- summary(fit)$coefficients
		var_names <- row.names(coeff);
		p <- coeff[, 4];
		p_adj <- p * weight[var_names]
		effect <- coeff[, 1];
		result <- data.frame(snp=var_names, p=p, p.w=p_adj, eff=effect, stringsAsFactors=F);
		result <- result[-1, ]
		return(result)
	}
}

