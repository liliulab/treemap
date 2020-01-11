library(BayesFactor)

## post-processing --> multi-threading enabled
## write a gene_name.RData output file
# input.info=finished.info=input.list.this[[n]]; input.folder='./'; weighted=''; covar.file=''; skip=F; 
collect.info <- function(input.info, input.folder='', weighted='', covar.file='', flag=0, delim=',', header=F) {
	gene.name <- '';
	group.file.path <- '';
	response.file.path <- '';
	geno.file.path <- '';
	snp.name.file.path <- '';
	truth.file.path <- '';
	lasso.file.path <- '';
	output.file.path <- '';
	log.file.path <- '';
	
	if(flag == 0) {
		gene.name <- input.info[1, 'gene.name'];
		group.file.path <- paste(input.folder, input.info[1, 'group.file'], sep='');
		response.file.path <- paste(input.folder, input.info[1, 'response.file'], sep='');
		geno.file.path <- paste(input.folder, input.info[1, 'geno.file'], sep='');
		snp.name.file.path <- paste(input.folder, input.info[1, 'snp.name.file'], sep='');	
		truth.file.path <- paste(input.folder, input.info[1, 'truth.file'], sep='');
		lasso.file.path <- paste(input.folder, input.info[1, 'lasso.file'], sep='');
		output.file.path <- paste(input.folder, input.info[1, 'output.file'], sep='');
		log.file.path <- paste(input.folder, input.info[1, 'log.file'], sep='');
	} else {
		gene.name <- input.info;
		group.file.path <- paste(input.folder, 'gene.group.maf05/', gene.name, '.group.csv.gz', sep='');
		response.file.path <- paste(input.folder, 'gene.exp.orig.maf05.adj/', gene.name, '.exp.csv.gz', sep='');		
		geno.file.path <- paste(input.folder, 'gene.geno.orig.maf05/', gene.name, '.geno.csv.gz', sep='');
		snp.name.file.path <- paste(input.folder, 'gene.feature.orig.maf05/', gene.name, '.feature.csv.gz', sep='');
		lasso.file.path <- paste(input.folder, 'tree.lasso.maf05/', gene.name, '.lasso.csv', sep='');
		output.file.path <- paste(input.folder, 'output.maf05/', gene.name, '.RData', sep='');
		log.file.path <- paste(input.folder, 'output.maf05/', gene.name, '.log.txt', sep='');
	}

	cat(paste(gene.name, 'collecting\n')); flush.console();
	
	snp.names <- read.table(snp.name.file.path, sep=delim, header=header, stringsAsFactors=F);	
	snp.anno <- snp.names
	if(ncol(snp.anno) == 1) {
		colnames(snp.anno) <- 'rs';
		snp.anno$id <- 1:nrow(snp.anno);
	} else {
		snp.anno$pos <- as.numeric(do.call(rbind, strsplit(snp.anno$rs, ':'))[, 2])
	}
	snp.names <- snp.names[, 1];
	snp.names <- data.frame(index=seq(1:length(snp.names)), name=snp.names, stringsAsFactors=F);
	
	truth <- data.frame(Causal_SNP=character(0), frequency=numeric(0), Allele_effect=numeric(0));		
	truth.index <- 0;
	if(nchar(truth.file.path) > 0 & file.exists(truth.file.path)) {
		truth <- read.table(truth.file.path, sep='\t', header=T, stringsAsFactors=F);
		truth <- truth[, which(colnames(truth) != 'Index')]
		truth <- merge(truth, snp.names, by.x='Causal_SNP', by.y='name', all=F);
		truth <- truth[order(truth$index), ]
		truth.index <- truth$index;
	}

	snp.names=snp.names$name;
	if(flag == 1) {
		snp.names <- paste('V', snp.names, sep='');
	}
	
	group <- c();
	if(file.exists(group.file.path)) {
		group <- read.table(group.file.path, sep=',', header=F, stringsAsFactors=F);	
		colnames(group) <- c('start', 'end', 'weight.size', 'weight.maf', 'weight.func', 'level', 'index', 'levelG');
		if(weighted == 'both') {
			group$weight <- group$weight.maf*(1-group$weight.maf) + (1-group$weight.func)*5 + 1
		} else if(weighted == 'maf') {
			group$weight <- group$weight.maf*(1-group$weight.maf) + 1
		} else if(weighted == 'function') {
			group$weight <- (1-group$weight.func)*5 + 1
		} else {
			group$weight <- 1;   ## no weight
		}
	}
	
	response <- read.table(response.file.path, sep=',', header=F, stringsAsFactors=F);		
	geno <- read.table(geno.file.path, sep=',', header=!header, stringsAsFactors=F);	
	input <- cbind(response, geno);
	colnames(input)[1] <- 'resp';
	
	covar <- c();
	if(nchar(covar.file) > 0) {
		covar <- read.table(covar.file, sep=',', header=F, stringsAsFactors=F);	
		colnames(covar) <- paste(rep('C', ncol(covar)), 1:ncol(covar), sep='')
	}	

	top <- c();
	top.file.path <- gsub('.RData', '.top.txt', output.file.path);
	if(file.exists(top.file.path)) {
		top <- read.table(top.file.path, sep='\t', header=T, stringsAsFactors=F);	
	}
	
	selected.within <- list();
	within.file.path <- gsub('.RData', '.within.RData', output.file.path);
	if(file.exists(within.file.path)) {
		load(within.file.path);	
	}
	
	selected.cross <- list();
	cross.file.path <- gsub('.RData', '.cross.RData', output.file.path);
	if(file.exists(cross.file.path)) {
		load(cross.file.path);	
	}
	
	selected.locus <- c();
	locus.file.path <- gsub('.RData', '.locus.txt', output.file.path);
	if(file.exists(locus.file.path)) {
		if(length(readLines(locus.file.path)) > 1) {
			selected.locus <- read.table(locus.file.path, sep='\t', header=T, stringsAsFactors=F);	
		}
	}

	selected.comb <- c();
	comb.file.path <- gsub('.RData', '.comb.txt', output.file.path);
	if(file.exists(comb.file.path)) {
		if(length(readLines(comb.file.path)) > 1) {
			selected.comb <- read.table(comb.file.path, sep='\t', header=T, stringsAsFactors=F);
		}
	}

	effect.file.path <- gsub('.RData', '.treemap.out', output.file.path);

	collected <- list(gene.name=gene.name, input=input, snp.anno=snp.anno, snp.names=snp.names, group=group, lasso.file.path=lasso.file.path, top.file.path=top.file.path, top=top, within.file.path=within.file.path, selected.within=selected.within, cross.file.path=cross.file.path, selected.cross=selected.cross, locus.file.path=locus.file.path, selected.locus=selected.locus, comb.file.path=comb.file.path, effect.file.path=effect.file.path, log.file.path=log.file.path, selected.comb=selected.comb, truth=truth);
	cat(paste(gene.name, 'collection done\n')); flush.console();
	return(collected);
}

prepare.postprec.GTEx <- function(input.folder.path, lasso.file.name='tree-lasso', output.folder.path='output/') {
	if(!dir.exists(output.folder.path)) {
		dir.create(output.folder.path, recursive=T);
	}
	
	file.info.list <- list();
	response.files = list.files(paste(input.folder.path, '/gene.exp.maf05.adj/', sep=''), pattern='*.exp.csv.gz');
	for(f in 1:length(response.files)) {
		gene.name <- gsub('.exp.csv.gz', '', response.files[f]);
		one.response.file.path <- paste(c(input.folder.path, '/gene.exp.maf05.adj/', response.files[f]), collapse='');
		one.geno.file.path <- gsub('exp', 'geno', one.response.file.path);
		one.geno.file.path <- gsub('.adj', '', one.geno.file.path);
		one.group.file.path <- gsub('geno', 'group', one.geno.file.path);
		one.name.file.path <- gsub('geno', 'feature', one.geno.file.path);
		one.lasso.file.path <- paste(input.folder.path, '/tree.lasso/', gene.name, '.', lasso.file.name, '.csv', sep='');
		one.output.file.path <- paste(c(output.folder.path, '/', gene.name, '.RData'), collapse='');

		file.info <- data.frame(gene.name=gene.name, group.file=one.group.file.path, response.file=one.response.file.path, geno.file=one.geno.file.path, lasso.file=one.lasso.file.path, snp.name.file=one.name.file.path, truth.file=NA, output.file=one.output.file.path, stringsAsFactors=F);	

		file.info.list[[f]] <- file.info;
	}
	return(file.info.list);
}
## testing
#setwd('~/my.scratch/projects/GTEx/brain.hippo'); 
#setwd('C:/Users/lliu80/Desktop/eQTL/data/GTEx/'); 
#source('../../my.functions.post.processing.R');
#source('../../Yang_conditional.r')
#source('../../Liu_conditional.r')
#input.folder.path='brain.hippo/'; output.folder.path='brain.hippo/output/'; lasso.file.name='tree-lasso.none'
#gtex.file.info.list <- prepare.postprec.GTEx(input.folder.path='brain.hippo/', output.folder.path='brain.hippo/output/', lasso.file.name='tree-lasso.none');
#finished.gtex.info <- mclapply(gtex.file.info.list, post.processing.1, input.folder='./', gz=T, weighted=T)

partition <- function(group, candidates, merge.single.leaf=F) {
	blocks <- c();
	candidates <- unique(candidates)
	for(r in candidates) {
		in.group <- group[which(group$weight.size > 1 & group$start <= r & group$end >= r), ];
		if(nrow(in.group) > 0) {
			blocks <- rbind(blocks, in.group[which.max(in.group$weight.size), c('start', 'end')]);  ## the largest group
		}		
	}
	blocks <- unique(blocks);
	
	list.blocks <- list();
	if(!is.null(blocks)) {
		for(b in 1:nrow(blocks)) {
			one.block <- blocks[b, ];
			in.one.block <- candidates[which(candidates %in% one.block$start:one.block$end)];
			in.one.block <- in.one.block[order(in.one.block)]
			exist <- unlist(lapply(list.blocks, function(x) length(x)==length(in.one.block) && sum(abs(x-in.one.block)) == 0))
			if(!TRUE %in% exist) {
				list.blocks[[b]] <- in.one.block;
			}
		}

		blocks.idx <- unique(unlist(apply(blocks, 1, function(x) return(x[1]:x[2]))))
		leaf <- candidates[which(!candidates %in% blocks.idx)];
		if(length(leaf) > 0) {
			if(merge.single.leaf) {   ## all leaves are merged to a single node
				list.blocks[[b+1]] <- leaf;
			} else {            ## each leaf is converted to a node.
				for (e in leaf) {
					b <- b+1;
					list.blocks[[b]] <- e;
				}
			}
		}
	}
	return(list.blocks);
}

rank.within.block <- function(data.all, top, one.block, cor.cutoff=0) {
	geno <- data.all[, -1];
	snp.names <- colnames(geno);
	ranked <- data.frame(idx=top, snp=snp.names[top], cor=1, rank=1, stringsAsFactors=F);
	one.block <- one.block[which(one.block != top)];
	if(length(one.block) > 0) {
		cor.block <- cor(geno[, one.block], geno[, top])^2;
		cor.block <- data.frame(cor=cor.block, idx=one.block, snp=snp.names[one.block])
		cor.block <- cor.block[which(cor.block$cor > cor.cutoff), ];
		if(nrow(cor.block) > 0) {
			cor.block$rank <- rank(-cor.block$cor, ties.method='min') + 1;
			cor.block <- cor.block[order(cor.block$rank), ]
			ranked <- rbind(ranked, cor.block);
		}
	}
	ranked$lead <- snp.names[top]
	ranked$p.value <- uni.unorder(data.all[, 1], data.all[, ranked$idx+1], snp.names[ranked$idx])$p.value;
	return(ranked);
}

is.colinear <- function(geno, idx.existing, idx.new) {
	data.geno <- data.frame(geno.new=geno[, idx.new], geno[, idx.existing]);
	fit <- lm(geno.new ~ ., data=data.geno, na.action=na.exclude);
	vif <- VIF(fit);
	cat(paste('vif', vif, '\n')); flush.console();
	flag <- F;
	if(vif >= 10) {
		flag <- T;
		cat(paste('multicolinearity [', idx.new, '-', idx.existing, ']...\n')); flush.console();
	}
	return(flag);
}

rsq.sequential <- function(response, geno, index, rsq.cutoff=0.01) {	
	data.input <- cbind(response, geno)
	rsq.pct <- c();
	flag <- T;
	while(flag) {	
		input.sub <- data.input[, c(1, index+1)];
		colnames(input.sub)[1] <- 'resp';
		model.input <- lm(resp ~ ., data=input.sub);
		p <- summary(model.input)$coeff[-1, 4];
		index <- index[order(p)];

		rsq.sub <- c();
		for(i in 1:(length(index))) {
			idx <- index[1:i];
			input.sub <- data.input[, c(1, idx+1)];
			colnames(input.sub)[1] <- 'resp';
			model.input <- lm(resp ~ ., data=input.sub);
			rsq.sub <- c(rsq.sub, summary(model.input)$r.squared)
		}
		rsq.diff <- rsq.sub - c(0, rsq.sub[-length(rsq.sub)])
		rsq.remain <- c(1, 1 - rsq.sub)
		rsq.pct <- rsq.diff/rsq.remain[-length(rsq.remain)]		
		rsq.pct <- data.frame(idx=index, rsq.diff=rsq.pct);
		rsq.pct <- rsq.pct[which(rsq.pct$rsq.diff >= rsq.cutoff), ]
		if(nrow(rsq.pct) > 0) {
			if(length(rsq.pct$idx) == length(index)) {
				flag <- F;
			} else {
				flag <- T;
				index=rsq.pct$idx
			}
		} else {
			flag <- F
		}
	}
	return(rsq.pct);
}

my_backward <- function(input, snp.names, index) {
	colnames(input) <- c('resp', snp.names);
	test.idx <- index;
	flag <- T;
	while(length(test.idx) > 0 & flag) {
		data.sub <- input[, c(1, test.idx+1)];
		fit <- lm(resp ~ ., data=data.sub, na.action=na.exclude); 
		pval <- summary(fit)$coefficients[-1, 4]
		keep.names <- names(pval);
		idx <- intersect(which.max(pval), which(pval > 0.05))
		if(length(idx) > 0) {
			keep.names <- names(pval)[-idx];
		}
		keep.idx <- which(snp.names %in% keep.names);
		if(length(keep.idx) == length(test.idx)) {
			flag = F;
		}
		test.idx <- test.idx[which(test.idx %in% keep.idx)];
	}
	return(test.idx);
}

# index shall include locus.
test_comb <- function(data.all, index, cn=2) {
	if(length(index) > 500) {
		index <- sample(index, 500, replace=F);
	}
	idx.single <- c();
	aic.single <- c();
	snp.names <- colnames(data.all)[-1];
	top <- uni.order(data.all[, 1], data.all[, index+1], snp.names[index])[1, ];
	if(top$p < 0.05) {
		idx.single <- which(snp.names == top$rs);
		aic.single <- AIC(lm(resp ~ ., data=data.all[, c(1, idx.single+1)], na.action=na.exclude));
	}
	
	idx.comb <- c();
	aic.comb <- c();	
	if(length(index) > 1) {
		comb <- combn(index, cn)
		result <- apply(comb, 2, function(x) {
			fit <- lm(resp ~ ., data=data.all[, c(1, x+1)], na.action=na.exclude);  ## change to BayesFactor gives the same results.
			s <- summary(fit)
			if(nrow(s$coeff) == 3 && s$coeff[2, 4] < 0.05 && s$coeff[3, 4] < 0.05) {
				return(c(x, s$adj.r.squared))
			}
		})	
		if(length(result) > 0) {
			if(is.list(result)) {
				result <- do.call(rbind, result);
			} else {
				result <- t(result);
			}
			idx.comb <- result[which.max(result[, 3]), ][1:2]
			aic.comb <- AIC(lm(resp ~ ., data=data.all[, c(1, idx.comb+1)], na.action=na.exclude));
		}
	}
	
	best <- c();
	if(length(aic.single) == 0 & length(aic.comb) > 0) {
		best <- idx.comb;
	} else if (length(aic.single) > 0 & length(aic.comb) == 0) {
		best <- idx.single;
	} else if (length(aic.single) > 0 && length(aic.comb) > 0 && aic.single < aic.comb) {
		best <- idx.single;
	} else if (length(aic.single) > 0 && length(aic.comb) > 0 && aic.comb < aic.single) {
		best <- idx.comb;
	}		
	return(best);
}

# pp=100110061
# snp.anno[which(snp.anno$pos > (pp-1000) & snp.anno$pos < (pp+1000)), ]

report.locus <- function(gene.name, selected.cross, input, snp.names, group, rsq.cutoff=0.01, locus.file.path) { 
	reported <- selected.cross;
	data.all <- input;
	
	top.cross.group <- reported$top.cross.group; 
	top.cross.linked <- reported$top.cross.linked; 
	top.cross.linked.weak <- reported$top.cross.linked.weak; 
	final.linked <- reported$final.linked; 
	final.linked.weak <- reported$final.linked.weak; 
	top.within.group <- reported$top.within.group;
	
	locus.block.all <- c();
#	list.reported <- list(first=top.cross.linked, second=top.cross.linked.weak, third=final.linked, fourth=final.linked.weak);
#	list.reported <- list(first=top.cross.linked, second=final.linked);
	list.reported <- list(first=top.cross.linked, second=top.within.group);
	for(k in 1:length(list.reported)) {
#		cat(k, '*****\n'); flush.console();
		level <- names(list.reported)[k];
		reported <- list.reported[[k]];
		locus.block <- c();
		if(length(reported) > 0) {
			reported <- reported[order(reported)]
			list.blocks <- partition(group=group, candidates=reported, merge.single.leaf=F) 	
#			cat('before in.block\n'); flush.console();
			ll <- 0;
			if(length(list.blocks) > 0) {
				for(b in 1:length(list.blocks)) {
					one.block <- list.blocks[[b]];
	#				cat(length(one.block), ','); flush.console();
					comb <- test_comb(data.all=data.all, index=one.block, cn=2);
					for(tp in comb) {
						ranked <- rank.within.block(data.all=data.all, top=tp, one.block=one.block, cor.cutoff=0.8);
						ll <- ll + 1;
						ranked$locus <- ll;
						ranked$level <- level;
						ranked$cross <- 'in';
						locus.block <- rbind(locus.block, ranked);
					}
				}  
			}
			
#			cat('before remains\n'); flush.console();
			remains <- reported[which(!reported %in% locus.block$idx)];
			flag <- !is.null(remains) && length(remains) > 0
			while(flag) {
				top <- uni.order(data.all[, 1], data.all[, remains+1], snp.names[remains])[1, 'rs'];
				top <- which(snp.names == top);
				fit <- lm(resp ~ ., data=data.all[, c(1, top+1)], na.action=na.exclude)
				coeff <- summary(fit)$coefficients
				p <- coeff[2, 4];
				if(p < 0.05 & summary(fit)$adj.r.squared > rsq.cutoff) {
					ranked <- rank.within.block(data.all, top, remains, cor.cutoff=0.8);
					ll <- ll + 1;
					ranked$locus <- ll;
					ranked$level <- level;
					ranked$cross <- 'remains';
					locus.block <- rbind(locus.block, ranked);
					remains <- remains[which(!remains %in% ranked$idx)]
					flag <- !is.null(remains) && length(remains) > 0
				} else {
					flag <- F;
				}
			}
					
			if(!is.null(locus.block)) {
				locus.block <- data.frame(gene.name=gene.name, locus.orig=as.numeric(locus.block$locus), rank=as.numeric(locus.block$rank), idx=as.numeric(locus.block$idx), snp=locus.block$snp, cor=as.numeric(locus.block$cor), level=locus.block$level, p.value=locus.block$p.value, lead=locus.block$lead, cross=locus.block$cross, stringsAsFactors=F);
				locus.block <- unique(locus.block)
				locus.block <- locus.block[order(locus.block$locus.orig, locus.block$level, locus.block$rank, locus.block$idx), ];
				
	#			cat('before top\n'); flush.console();
				locus.block.top <- c();
				for(lc in 1:max(locus.block$locus.orig)) {
					one <- locus.block[which(locus.block$locus.orig == lc & locus.block$snp == locus.block$lead), ][1, ];
					locus.block.top <- rbind(locus.block.top, one);
				}
				if(nrow(locus.block.top) > 1) {
					fit <- lm(resp ~ ., data=data.all[, c(1, locus.block.top$idx+1)], na.action=na.exclude)
					coeff <- summary(fit)$coeff;
					order <- order(coeff[-1, 4])
					locus.block.top <- locus.block.top[order, ];
				}
				locus.block.top$locus <- 1:nrow(locus.block.top);
				locus.block <- merge(locus.block, locus.block.top[, c('snp', 'locus')], by.x='lead', by.y='snp', all=T);
				locus.block <- locus.block[, c('gene.name', 'locus', 'rank', 'idx', 'snp', 'cor', 'lead', 'cross', 'level')]
				locus.block <- locus.block[order(locus.block$locus, locus.block$rank, locus.block$idx), ];
				locus.block$flag <- ifelse(locus.block$locus <= length(top.cross.group), '*', '-');		
				locus.block$flag <- ifelse(locus.block$cross == 'in', paste(locus.block$flag, '*', sep=''), locus.block$flag);		
				
				candidates <- locus.block[which(locus.block$rank == 1 & locus.block$snp == locus.block$lead), 'idx'];
				rsq.diff <- rsq.sequential(data.all[, 1], data.all[, -1], candidates, rsq.cutoff)
				locus.block <- merge(locus.block, rsq.diff, by='idx', all=T);
				locus.block <- locus.block[order(locus.block$locus, locus.block$rank, locus.block$idx), ];
				
				locus.block.all <- rbind(locus.block.all, locus.block);
			}
		} ## end if length(reported)>0
	}
	cat(paste('   locus.block.all:', nrow(locus.block.all), '\n')); flush.console();
	
	write.table(locus.block.all, locus.file.path, sep='\t', row.names=F, col.names=T, quote=F);
	return(locus.block.all);
}
## testing
# all.locus.blocks.1 <- do.call(rbind, mclapply(file.info.list.1, report.locus, input.folder='./'));
# all.locus.blocks.1.L4 <- do.call(rbind, mclapply(file.info.list.1.L4, report.locus, input.folder='./'));
# gtex.locus.blocks <- do.call(rbind, mclapply(gtex.file.info.list, report.locus, input.folder='./'));
# mclapply(gtex.file.info.list, report.locus, input.folder='./', write=T, mc.cores=20);

combine.locus <- function(gene.name, selected.locus, selected.cross, input, snp.names, group, rsq.cutoff=0.01, comb.file.path, snp.anno) {
	reported <- selected.cross;
	report <- selected.locus;
	data.all <- input;
	
	report.1 <- report[which(report$level == 'first' & report$rsq.diff >= rsq.cutoff), ]
	report.2 <- report[which(report$level == 'second' & report$rsq.diff >= rsq.cutoff), ]
	report.3 <- report[which(report$level == 'third' & report$rsq.diff >= rsq.cutoff), ]
	report.4 <- report[which(report$level == 'fourth' & report$rsq.diff >= rsq.cutoff), ]
#	report.list <- list(first=report.1, second=report.2, third=report.3, fourth=report.4);
	report.list <- list(first=report.1, second=report.2);
	report.all <- do.call(rbind, report.list);
	candidates <- unique(report.all$idx)
	where.signal <- '';
	if(length(candidates) > 0) {
		if(length(candidates) <= 10) {
			bf <- regressionBF(resp ~ ., data=data.all[, c(1, candidates+1)], progress=F)
			bf.best <- head(bf, n=1)
			best <- as.character(names(bf.best)$num)
			best <- trimws(unlist(strsplit(best, '\\+')))
			best <- which(snp.names %in% best)
			report.bf <- report.all[which(report.all$idx %in% best), ]
			report.bf$level <- 'bf';
			report.bf$rank <- 1;
			report.bf$locus <- 1;
			report.bf$cross <- 'in';
			report.bf$flag <- '--*';
			report.bf$rsq.diff <- 1;
			report.bf$lead <- report.bf$snp;
			report.bf$cor <- 1;
			report.bf <- unique(report.bf);
			report.list[[length(report.list)+1]] <- report.bf;
			where.signal <- 'bf';
		} else if(length(candidates) > 10){
			best <- stepwise.leaf(data.all, candidates, snp.names, type='backward')
			report.aic <- report.all[which(report.all$idx %in% best), ]
			report.aic$level <- 'aic';
			report.aic$rank <- 1;
			report.aic$locus <- 1;
			report.aic$cross <- 'in';
			report.aic$flag <- '--*';
			report.aic$rsq.diff <- 1;
			report.aic$lead <- report.aic$snp;
			report.aic$cor <- 1;
			report.aic <- unique(report.aic);
			report.list[[length(report.list)+1]] <- report.aic;
			where.signal <- 'aic'
		}
		best <- -Inf;
		best.l <- 0;
		for(l in 1:length(report.list)) {
			rpt <- report.list[[l]];
			if(nrow(rpt) > 0) {
				fit <- lm(resp ~ ., data=data.all[, c(1, rpt$idx+1)]);
				x <- summary(fit)$adj.r.squared;
				# a <- -AIC(fit); 
				# cat(l, ':', x, '-', a, '\n'); 
				if(x > best) {
					best <- x;
					best.l <- l;
				}
			}
		}
		
		if(best.l > 0) {
			locus.comb <- report.list[[best.l]]
			
			locus.block <- c();
			if(nrow(locus.comb) > 0) {
				base <- reported$final.linked;
				list.blocks <- partition(group=group, candidates=base, merge.single.leaf=F) 	
				candidates <- unique(locus.comb$idx)
				ll <- 0;
				for(t in 1:length(candidates)) {
					locus <- candidates[t];
					for(b in 1:length(list.blocks)) {
						one.block <- list.blocks[[b]];
						if(locus %in% one.block) {
							ranked <- rank.within.block(data.all, locus, one.block, cor.cutoff=0.5);
							ranked$level <- paste('comb', best.l, sep='.');
							ranked$cross <- 'in';  ## everything "in", nothing "out"
							ranked <- ranked[which(ranked$p.value < 0.01), ]
							if(nrow(ranked) > 0) {
								ranked$rank <- 1:nrow(ranked);
								ll <- ll + 1;
								ranked$locus <- ll;
								locus.block <- rbind(locus.block, ranked);
							}
						}
					}  
				}
				if(!is.null(locus.block)) {
					locus.block$rsq.diff <- 1;
					locus.block$gene.name <- gene.name;
					locus.block$flag <- '*^';
					locus.block <- merge(locus.block, snp.anno, by.x='idx', by.y='id', all.x=T, all.y=F);
					locus.block <- locus.block[order(locus.block$locus, locus.block$rank), ]		
					write.table(locus.block, comb.file.path, sep='\t', row.names=F, col.names=T, quote=F);
					return(locus.block);
				} else {
					where.signal <- 'locus.block=0';
				}
			} else {
				where.signal <- 'locus.comb=0';
			}
		} else {
			where.signal <- 'best=0';
		}
	} else {
		where.signal <- 'candidate=0';
	}
	cat(paste('   where.signal:', where.signal, ' -- ', gene.name, '\n')); flush.console();
}


## estimate effect size - gtex.file.info.list[[21726]]; 
## input.folder='./'; p.cutoff=1e-4
estimate.effect <- function(gene.name, input.folder='./', p.cutoff=1e-4) {
	output.file <- paste(input.folder, '/', 'output.maf05/', gene.name, '.comb.txt', sep='');
	geno.file <- paste(input.folder, '/', 'gene.geno.orig.maf05/', gene.name, '.geno.csv.gz', sep='');
	exp.file <- paste(input.folder, '/', 'gene.exp.orig.maf05.adj/', gene.name, '.exp.csv.gz', sep='');
	feature.file <- paste(input.folder, '/', 'gene.feature.orig.maf05/', gene.name, '.feature.csv.gz', sep='');

	ln <- length(readLines(output.file))
	effect <- c();
	if(ln > 1) {
		cat(paste('estimating effect ...', gene.name, '\n', sep='')); flush.console();
		output <- read.table(output.file, sep='\t', header=T, stringsAsFactors=F);	
		geno <- read.table(geno.file, sep=',', header=F, stringsAsFactors=F);	
		response <- read.table(exp.file, sep=',', header=F, stringsAsFactors=F);	
		feature <- read.table(feature.file, sep='\t', header=T, stringsAsFactors=F);	
	
		if(min(output$p.value) <= p.cutoff) {
			snps <- unique(output[which(output$rank == 1), 'snp'])
			flag <- T;
			while(flag & length(snps) > 0) {
				data.sub <- data.frame(response[, 1], geno[, snps]);
				colnames(data.sub) <- c('exp', snps);
				fit <- lm(exp ~ ., data=data.sub);
				effect <- as.data.frame(summary(fit)$coeff[, ])
				effect <- effect[-1, ]
				temp <- effect[order(-effect[, 4]), ]
				if(temp[1, 4] < 0.01) {
					flag <- F
				} else {
					snps <- snps[which(snps != rownames(temp)[1])]
				}
			}
			if(nrow(effect) > 0) {
				effect <- data.frame(var=row.names(effect), effect);
				effect <- cbind(gene.name, effect);
			}
		}
	}
	return(effect);
}

## input.folder='./'; gene.name='ENSG00000124102'; snps=c('V4034', 'V8452'))
estimate.effect.snp <- function(gene.name, input.folder='./', snps=NULL) {
	output.file <- paste(input.folder, '/', 'output.maf05/', gene.name, '.comb.txt', sep='');
	geno.file <- paste(input.folder, '/', 'gene.geno.orig.maf05/', gene.name, '.geno.csv.gz', sep='');
	exp.file <- paste(input.folder, '/', 'gene.exp.orig.maf05.adj/', gene.name, '.exp.csv.gz', sep='');
	feature.file <- paste(input.folder, '/', 'gene.feature.orig.maf05/', gene.name, '.feature.csv.gz', sep='');

	output <- read.table(output.file, sep='\t', header=T, stringsAsFactors=F);	
	geno <- read.table(geno.file, sep=',', header=F, stringsAsFactors=F);	
	response <- read.table(exp.file, sep=',', header=F, stringsAsFactors=F);	
	feature <- read.table(feature.file, sep='\t', header=T, stringsAsFactors=F);	
	
	data.sub <- data.frame(response[, 1], geno[, snps]);
	colnames(data.sub) <- c('exp', snps);
	fit <- lm(exp ~ ., data=data.sub);
	effect <- summary(fit)$coeff[, ];
	effect <- data.frame(var=row.names(effect), effect);
	effect <- cbind(gene.name, effect[-1, ]);
	
	return(effect);
}

estimate.effect.comb <- function(gene.name, selected.comb, input, effect.file.path, p.cutoff=1e-4) {
	if(!is.null(selected.comb) && nrow(selected.comb) > 0) {
		cat(paste('estimating effect ...', gene.name, '\n', sep='')); flush.console();
		output <- selected.comb[, which(!colnames(selected.comb) %in% c('idx', 'level', 'cross', 'flag', 'rs', 'p', 'p.w'))];	
		geno <- input[, -1];
		response <- input[, 1];
	
		effect <- c();
		if(min(output$p.value, na.rm=T) <= p.cutoff) {
			snps <- unique(output[which(output$rank == 1), 'snp'])
			flag <- T;
			while(flag & length(snps) > 0) {
				data.sub <- data.frame(response, geno[, snps]);
				colnames(data.sub) <- c('exp', snps);
				fit <- lm(exp ~ ., data=data.sub);
				effect <- as.data.frame(summary(fit)$coeff[, ])
				effect <- effect[-1, ]
				temp <- effect[order(-effect[, 4]), ]
				if(temp[1, 4] < 0.01) {
					flag <- F
				} else {
					snps <- snps[which(snps != rownames(temp)[1])]
				}
			}
			if(nrow(effect) > 0) {
				effect <- data.frame(snp=row.names(effect), rank=1, effect.est=effect[, 1], effect.p=effect[, 4]);
				effect <- cbind(gene.name, effect);
				comb.effect <- merge(output, effect, by=c('gene.name', 'snp', 'rank'), all=T);
				comb.effect <- comb.effect[order(comb.effect$locus, comb.effect$rank), ];
				write.table(format(comb.effect, digits=3), effect.file.path, sep='\t', row.names=F, col.names=T, quote=F);
				return(comb.effect);
			}
		}
	}
}
