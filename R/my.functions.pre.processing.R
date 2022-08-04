
## pre-processing: top down - graph
library(igraph)
library(intervals)
find.blocks.graph <- function(cor.feature, cutoffs, maf=NULL, func=NULL) {
	cor <- abs(cor.feature)
	colnames(cor) <- 1:nrow(cor);
	rownames(cor) <- 1:nrow(cor);
	diag(cor) <- 0;
	ld <- c();
	for(level in 1:length(cutoffs)) {
		this.cor <- ifelse(cor > cutoffs[level], 1, 0)
		net <- as.undirected(graph_from_adjacency_matrix(this.cor))
#		plot(net);

		cluster <- cluster_fast_greedy(net, merges=F, modularity=T, membership=T, weights=NULL);
#		plot(cluster, net);
		cluster.size <- as.numeric(sizes(cluster));
		large.cluster <- cluster[which(cluster.size > 1)]
		cluster.new <- list();
		cnt <- 1;
		if(length(large.cluster) > 0) {
			for(u in 1:length(large.cluster)) {
				one.cluster <- as.numeric(large.cluster[[u]])
				one.cluster <- one.cluster[order(one.cluster)]
				idx <- 1:(length(one.cluster)-1);
				distance <- one.cluster[idx + 1] - one.cluster[idx];
				break.point <- which(distance > 100);
				if(length(break.point) > 0) {
					prev <- 1;
					for(b in break.point) {
						cluster.new[[cnt]] <- one.cluster[prev:b];
						prev <- b + 1;
						cnt <- cnt + 1;
					}
					cluster.new[[cnt]] <- one.cluster[prev:length(one.cluster)];
					cnt <- cnt + 1;
				} else {
					cluster.new[[cnt]] <- one.cluster;
					cnt <- cnt + 1;
				}
			}
		}
		if(length(cluster.new) > 0) {
			cluster.new.size <- unlist(lapply(cluster.new, length));
			large.cluster.new <- cluster.new[which(cluster.new.size > 1)]
			cluster.min <- unlist(lapply(large.cluster.new, min));
			cluster.max <- unlist(lapply(large.cluster.new, max));
			cluster.range <- data.frame(min=cluster.min, max=cluster.max);
			
			m <- cluster.range[order(cluster.range$min), ]
			prev <- nrow(m);
			this <- 0;
			while(prev != this) {
				prev <- nrow(m);
				interval.1 <- Intervals(m, closed=T, type='Z')
				interval.2 <- interval.1
				overlap <- interval_overlap(interval.1, interval.2);
				overlap <- unique(overlap);
				results.range <- c();
				for(l in 1:length(overlap)) {
					ov <- overlap[[l]];
					ov.range <- range(unlist(m[ov, ]));
					results.range <- rbind(results.range, ov.range);	
				}
				results.range <- unique(results.range);
				m <- results.range;
				rownames(m) <- 1:nrow(m);
				this <- nrow(m);
			}
			ld <- rbind(ld, cbind(m, level));
		}
	}
	
	if(!is.null(ld)) {
		ld.maf <- c();
		if(is.null(maf)) {
			ld.maf <- 1.5;
		} else {
			ld.maf <- apply(ld, 1, function(x) { mean(maf[x[1]:x[2]], na.rm=F); });
		}
		
		ld.func <- c();
		if(is.null(func)) {
			ld.func <- 0;
		} else {
			ld.func <- apply(ld, 1, function(x) { mean(func[x[1]:x[2]], na.rm=F); });
		}
		
		ld <- data.frame(start=ld[, 1], end=ld[, 2], size=ld[, 2]-ld[, 1]+1, maf=ld.maf, func=ld.func, level=ld[, 3]);
	} else {
		ld <- data.frame(start=numeric(0), end=numeric(0), size=numeric(0), maf=numeric(0), func=numeric(0), level=numeric(0));
	}
	return(ld);
}
## testing
#	cor.feature=cor.feature.4; offset.start=1; offset.end=nrow(cor.feature.4); cutoff=0.8; distance=20; maf=rep(1, offset.end-offset.start+1);
#	cor.feature=cor.feature.4; offset.start=8; offset.end=15; cutoff=0.99; distance=5; maf=rep(1, offset.end-offset.start+1);
# setwd('C:/Users/lliu80/Desktop/eQTL/data/GTEx/brain.hippo')
# genotype <- read.table(gzfile('ENSG00000000003.geno.csv.gz'), header=T, sep=','); 
# cor.feature <- cor(genotype, use='complete.obs'); cor.feature <- cor.feature^2; maf <- colSums(genotype)/(2*nrow(genotype)); cutoffs <- c(0.8, 0.9, 0.95, 0.99, 0.998);

library(DMwR);  ## for knnImputation
impute.missing <- function(pheno, geno) {
	missing.pheno.idx <- which(is.na(pheno));
	na.geno.idx <- apply(geno, 2, function(x) { length(which(is.na(x))) });
	low.geno.idx <- which(na.geno.idx > 0.1 * nrow(geno))
	if(length(low.geno.idx) > 0) {
		geno <- geno[, -low.geno.idx];
	}
	
	geno.imputed <- tryCatch( {
		t(knnImputation(t(geno),k=1, meth='median'))
	}, error=function(x) {
		tryCatch( {
			knnImputation(geno,k=1, meth='median')
		}, error=function(y) {
			geno;
		})
	})
	
	missing.geno.idx <- which(is.na(rowSums(geno.imputed)));
	missing.index <- unique(c(missing.pheno.idx, missing.geno.idx));	
	
	if(length(missing.index) > 0) {
		pheno <- pheno[-missing.index];
		geno <- geno.imputed[-missing.index, ]
	} else {
		geno <- geno.imputed;
	}
	
	## remove SNPs where all samples have the same genotype 
	uu <- apply(geno, 2, function(x) length(unique(x)));
	uu.index <- which(uu == 1);
	if(length(uu.index) > 0) {
		geno <- geno[, -uu.index];
	}

	return(list(pheno=pheno, geno=geno));
}

normalize <- function(x) {
	y <- (x - mean(x, na.rm=T))/sd(x, na.rm=T);
	return(y);
}

library(parallel)
pre.processing.1 <- function(input.info, input.folder='', gz=TRUE, size='all', skip=F, ld=T) {
	name <- c();
	data.file <- c();
	pheno.file <- c();
	geno.file <- c();
	feature.file <- c();
	output.file.prefix <- c();
	sep <- ',';
	if(length(input.info) > 1) {
		name <- input.info[1, 'name'];
		data.file <- paste(input.folder, input.info[1, 'data.file'], sep='');
		feature.file <- paste(input.folder, input.info[1, 'feature.file'], sep='');
		output.file.prefix <- input.info[1, 'output.file.prefix'];
		sep <- input.info[1, 'sep'];
		
		test.out <- paste(output.file.prefix, 'geno.norm.csv.gz', sep='.');
		if(file.exists(test.out) & skip) {
			return(0);
		}
	} else {
		name <- input.info;
		pheno.file <- paste(input.folder, '/gene.exp.orig.maf05.adj/', name, '.exp.csv.gz', sep='');
		geno.file <- paste(input.folder, '/gene.geno.orig.maf05/', name, '.geno.csv.gz', sep='');
		feature.file <- paste(input.folder, '/gene.feature.orig.maf05/', name, '.feature.csv.gz', sep='');
		
		test.out <- gsub('.orig', '', geno.file);
		if(file.exists(test.out) & skip) {
	#		cat('exist. return(0)\n'); flush.console();
			return(0);
		}
	}		
	cat(paste('pre-processing', name, '... \n', sep=' ')); flush.console();	

	phenotype <- c();
	genotype <- c();
	if(!is.null(data.file) && nchar(data.file) > 0) {
		data <- read.table(data.file, sep=sep, header=T, row.names=1, stringsAsFactors=F);			
		phenotype <- data[, 1];
		genotype <- data[, -1];
	} else if(nchar(pheno.file) > 0 & nchar(geno.file) > 0) {
		phenotype <- read.table(pheno.file, sep=sep, header=F, stringsAsFactors=F);			
		genotype <- read.table(geno.file, sep=sep, header=F, stringsAsFactors=F);	
	}

	feature <- c();
	if(nchar(feature.file) > 2) {
		feature <- read.table(feature.file, sep='\t', header=T, stringsAsFactors=F);	
		if(ncol(feature) == 1) {
			colnames(feature) <- c('rs');
			feature$func <- 0;
		} else if(ncol(feature) == 2) {
			colnames(feature) <- c('rs', 'func');
		} else if(ncol(feature) == 3) {
			colnames(feature) <- c('id.before', 'idx.orig', 'rs');
			feature$func <- 0;
			feature$rs <- paste(feature$rs, feature$id.before, sep=':');
		} else if(ncol(feature) == 4) {
			colnames(feature) <- c('id.before', 'idx.orig', 'rs', 'func');
			feature$rs <- paste(feature$rs, feature$id.before, sep=':');
		} else if(ncol(feature) == 6) {
			colnames(feature) <- c('id.before', 'idx.orig', 'rs', 'dnase', 'tfbs', 'func');
			feature$rs <- paste(feature$rs, feature$id.before, sep=':');
		}
	} else {
		feature <- data.frame(rs=colnames(genotype), func=0);
	}
	if(size == 'small' & nrow(feature) > 10000) {
#		cat(nrow(feature), ' ---> skipping\n');
		return(name);
	} else if(size == 'big' & (nrow(feature) <= 10000 | nrow(feature) > 20000)) {
#		cat(nrow(feature), ' ---> skipping\n');
		return(name);
	} else if(size == 'huge' & nrow(feature) <= 20000) {
#		cat(nrow(feature), ' ---> skipping\n');
		return(name);
	} else {
		cat(paste(nrow(feature), '\n', sep=' '));
	}
	
	cat(paste('   input', nrow(feature), 'features', length(phenotype), 'phenotypes', nrow(genotype), 'x', ncol(genotype), 'genotypes ...\n', sep=' ')); flush.console();
	colnames(genotype) <- feature$rs;
	genotype[genotype == -1] <- NA;
	imputed <- impute.missing(phenotype, genotype);
	phenotype <- imputed[['pheno']];
	genotype <- imputed[['geno']];	
	cat(paste('   after imputation', nrow(genotype), 'x', ncol(genotype), 'genotypes \n', sep=' ')); flush.console();
	if(ncol(genotype) == 0) {
		cat(paste('   0 genotype, ---> nothing to work with.\n', sep=' '));
		return(name);
	}
	
	if(ld) {
		cor.feature <- cor(genotype, use='complete.obs');
		cor.feature <- cor.feature^2;
		maf <- colSums(genotype)/(2*nrow(genotype));
		feature.cnt <- ncol(genotype)
		snp.names <- data.frame(id=1:feature.cnt, rs=colnames(genotype));
		
		func <- merge(snp.names, feature, by='rs', all.x=T, all.y=F);
		func <- func[order(func$id), ];
		func$func <- ifelse(is.na(func$func), 0, func$func);
		func <- as.numeric(func$func);
		
		## LD graph
		cutoffs <- c(0.8, 0.85, 0.9, 0.95, 0.99, 0.998);
		blocks <- find.blocks.graph(cor.feature, cutoffs, maf, func);
		## re-order clusters and remove redundancy
		blocks.multiple <- data.frame();  ## tag, blocks
		for(level in length(cutoffs):1) {
			sub.blocks <- blocks[which(blocks$level == level), ];
			sub.blocks <- sub.blocks[order(sub.blocks$start), ];
			if(nrow(blocks.multiple) > 0) {  ## comment out this block if not removing redundancy
				existing.blocks <- paste(blocks.multiple[, 1], blocks.multiple[, 2], blocks.multiple[, 3])
				new.blocks <- paste(sub.blocks[, 1], sub.blocks[, 2], sub.blocks[, 3])
				remove.index <- which(new.blocks %in% existing.blocks);    
				if(length(remove.index) > 0) {
					sub.blocks <- sub.blocks[-remove.index, ]
				}
			}
			if(nrow(sub.blocks) > 0) {
				blocks.multiple <- rbind(blocks.multiple, sub.blocks);
			}
		}	
		if(nrow(blocks.multiple) > 0) {
			blocks.multiple$level <- max(blocks.multiple$level) - blocks.multiple$level + 1;
		}
		## add leave nodes
		blocks.multiple <- rbind(data.frame(start=1:feature.cnt, end=1:feature.cnt, size=1, maf=maf, func=func, level=0), blocks.multiple);  
		blocks.multiple$index <- 1:nrow(blocks.multiple);
		blocks.multiple$group <- ifelse(blocks.multiple$level == 0, 300, 100 + blocks.multiple$level);  
		## start, end, size, avg.maf, avg.func, (level,  index, label)
	}
	
	genotype.norm <- apply(genotype, 2, normalize);
	genotype.norm <- round(genotype.norm, digits=3)
	feature.name <- colnames(genotype.norm);	
	pheno.out <- c();
	geno.out <- c();
	geno.norm.out <- c();
	feature.out <- c();
	group.out <- c();
	if(!is.null(output.file.prefix) && nchar(output.file.prefix) > 0) {
		pheno.out <- paste(output.file.prefix, 'response.csv.gz', sep='.');
		geno.out <- paste(output.file.prefix, 'geno.csv.gz', sep='.');
		geno.norm.out <- paste(output.file.prefix, 'geno.norm.csv.gz', sep='.');
		feature.out <- paste(output.file.prefix, 'feature.csv.gz', sep='.');
		group.out <- paste(output.file.prefix, 'group.csv.gz', sep='.');
	} else {
		pheno.out <- gsub('.orig', '', pheno.file);
		geno.out <- gsub('.orig', '', geno.file);
		geno.norm.out <- gsub('geno', 'geno.norm', geno.out);
		feature.out <- gsub('.orig', '', feature.file);
		group.out <- gsub('feature', 'group', feature.out);
	}
	write.table(phenotype, gzfile(pheno.out), sep=',', row.names=F, col.names=F, quote=F)
	write.table(genotype, gzfile(geno.out), sep=',', row.names=F, col.names=T, quote=F)
	write.table(genotype.norm, gzfile(geno.norm.out), sep=',', row.names=F, col.names=F, quote=F)
	write.table(feature.name, gzfile(feature.out), sep=',', row.names=F, col.names=F, quote=F);	
	if(ld) {
		write.table(blocks.multiple, gzfile(group.out), sep=',', row.names=F, col.names=F, quote=F);
	}
	gc();
	return(0);
}
## testing
# simulation
# finished <- mclapply(input.list, pre.processing.1, gz=F)
# gtex
#lapply(input.list[1:3], pre.processing.1, gz=T, size='all')
#finished <- mclapply(input.list[61:51653], pre.processing.1, gz=T)
# finished <- mclapply(input.list[1:100], mc.cores=20, pre.processing.1, gz=T, size='all')
# sbatch batch.preprocessing.r.2.bash

## input.folder.path='simulated/1_causal/0_0.1/'; output.folder.path='output/1_causal/0_0.1/'; weighted=''
create.processing.list.simulation.2019 <- function(input.folder.path, pattern='.simulated.txt', output.folder.path, weighted='') {
	result <- c();
	input.files <- list.files(input.folder.path, pattern);
	gene.names <- gsub(pattern, '', input.files);
	for(f in 1:length(input.files)) {
		one.input.file <- input.files[f];
		gene.name <- gene.names[f]
		cat('listing ...', f, '...', gene.name, '\n'); flush.console();
		one.input.file.path <- paste(input.folder.path, one.input.file, sep='/');
		data.file <- one.input.file.path;
		
		feature.file <- '';
		if(weighted == 'func') {
			feature.file = gsub(pattern, '.anno.txt', data.file);
		}
		output.folder <- output.folder.path;
		if(!dir.exists(output.folder)) {
			dir.create(output.folder, recursive=T);
		}
		output.file.prefix <- paste(c(output.folder, '/', gene.name), collapse='');

		response.file = paste(output.file.prefix, '.response.csv.gz', sep='');
		lasso.file = paste(output.file.prefix, '.tree-lasso.csv', sep='');
		name.file = paste(output.file.prefix, '.feature.csv.gz', sep=''); 
		group.file = paste(output.file.prefix, '.group.csv.gz', sep=''); 
		geno.file = paste(output.file.prefix, '.geno.csv.gz', sep=''); 
		truth.file = gsub(pattern, '.causal.txt', one.input.file.path); 
		output.file = paste(output.file.prefix, '.RData', sep=''); 
		log.file = paste(output.file.prefix, '.log.txt', sep=''); 
		if(nchar(weighted) > 0) {
			output.file <- gsub('RData', paste(weighted, 'RData', sep='.'), output.file);
			log.file <- gsub('log.txt', paste(weighted, 'log.txt', sep='.'), log.file);
		}
		result <- rbind(result, c(f, data.file, feature.file, output.file.prefix, sep='\t', gene.name=gene.name, group.file=group.file, response.file=response.file, geno.file=geno.file, lasso.file=lasso.file, snp.name.file=name.file, truth.file=truth.file, output.file=output.file, log.file=log.file, output.path=output.folder));
	}
	result <- as.data.frame(result, stringsAsFactors=F);
	colnames(result) <- c('name', 'data.file', 'feature.file', 'output.file.prefix', 'sep', 'gene.name', 'group.file', 'response.file', 'geno.file', 'lasso.file', 'snp.name.file', 'truth.file', 'output.file', 'log.file', 'output.path');
	result <- split(result, seq(nrow(result))); ## convert to list
	return(result);
}

create.batch.processing.list <- function(input.folder.path, pattern, output.folder.path, weighted='') {
	input.info.list <- create.processing.list.simulation.2019(input.folder.path, pattern, output.folder.path, weighted=weighted)
	return(input.info.list);
}

adj.response <- function(response.file, covar.file) {
	cat('processing', response.file, '...\n'); flush.console();
	geno.file <- gsub('exp', 'geno', response.file);
	if(file.exists(geno.file)) {
		covar <- read.table(covar.file, sep=',', header=F, stringsAsFactors=F);
		colnames(covar) <- paste(rep('C', ncol(covar)), 1:ncol(covar), sep='')
		output.file <- gsub('gene.exp.orig.maf05', 'gene.exp.orig.maf05.adj', response.file);

		response <- read.table(gzfile(response.file), sep=',', header=F, stringsAsFactors=F);	
		geno <- read.table(gzfile(geno.file), sep=',', header=F, stringsAsFactors=F);	
		result <- c();
		for(idx in 1:ncol(geno)) {
			input <- cbind(response, covar, geno[, idx]);
			colnames(input)[1] <- 'resp';
			model <- lm(resp ~ ., data=input);
			coeff <- summary(model)$coeff;
			p <- coeff[nrow(coeff), 4];
			result <- c(result, p);
		}
		result <- data.frame(idx=1:length(result), p=result);
		result <- result[order(result$p), ];
		
		idx.best <- result[1, 'idx'];
		input <- cbind(response, covar, geno[, idx.best]);
		colnames(input)[1] <- 'resp';
		model <- lm(resp ~ ., data=input);
		coeff <- summary(model)$coeff;
		b <- coeff[1:(nrow(coeff)-1), 1];
		resp.adj <- c();
		for(i in 1:nrow(covar)) {
			x <- c(1, unlist(covar[i, ]));
			one.adj <- response[i, 1] - sum(b * x);
			resp.adj <- c(resp.adj, one.adj); 
		}
		write.table(resp.adj, gzfile(output.file), sep=',', quote=F, row.names=F, col.names=F);
	}
}
## testing
# setwd('~/my.scratch/projects/GTEx/brain.hippo')
# source('../my.functions.pre.processing.R')
# response.files <- paste('gene.exp.orig.maf05/', list.files('gene.exp.orig.maf05/'), sep='/');
# covar.file <- 'Brain_Hippocampus.v7.covariates.csv';
# finished <- mclapply(response.files, mc.cores=28, adj.response, covar.file=covar.file)

## input.folder='simulated/1_causal/0_0.1/'; pattern='.simulated.txt'; anno.file.name='../simulation/real_genotype/DHS_variants_in_simulated_cis_regions';
annotate.features = function(input.folder, pattern='.simulated.txt', anno.file.name) {
	anno = read.table(anno.file.name, header=F, sep='\t', stringsAsFactors=F)
	func.var = anno[, 1]
	input.files = list.files(input.folder, pattern=pattern);
	input.files = paste(input.folder, input.files, sep='/');
	output.files = gsub(pattern, '.anno.txt', input.files)
	for(f in 1:length(input.files)) {
		data = read.table(input.files[f], sep='\t', header=T, nrows=3)
		features = data.frame(rs=colnames(data[, -c(1,2)]), func=0);
		features[which(features$rs %in% func.var), 'func'] = 1;
		write.table(features, output.files[f], sep='\t', row.names=F, quote=F);
	}
}

