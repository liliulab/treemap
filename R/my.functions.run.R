library(parallel);

# i=3; n=136; input.info <- input.list.all[[i]][[n]]; input.folder='./'; write=F; rsq.cutoff=0.01; weighted=''; input.info$gene.name; flag=0; time.limit=0; delim=','; header=F; skip=0;
# input.info="ENSG00000163071"; weighted=''; steps=3:8; flag=1; time.limit=0; delim='\t'; header=T; skip=1;
run.1 <- function(input.info, weighted='', steps=1:9, flag=0, time.limit=0, delim=',', header=F, skip=0, p.cutoff=1e-4) {
	## pre-processing
	time.0 <- Sys.time()
	stop.signal <- F;
	if(1 %in% steps) {
		stop.signal <- tryCatch({
			pre.processing.1(input.info);
		}, error = function(e) {
			stop.signal <- T;
			cat(paste(input.info$gene.name, ' preprocessing error\n'));
			print(e);
			return(stop.signal);
		})
	}
	time.pre <- Sys.time()
	
	if(!stop.signal) {
		## outer layer
		if(2 %in% steps) {
			if(flag==0) {
				run.tree.lasso(name=input.info[1, 'response.file'], flag=flag, flag_string='response', output_file_name='tree-lasso', time.limit=time.limit, weight_flag=0, skip=skip);
			} else {
				run.tree.lasso(name=input.info, flag=flag, flag_string='exp', output_file_name='tree.lasso.maf05', time.limit=time.limit, weight_flag=0, skip=skip);
			}
		}
		time.lasso <- Sys.time()
		cat('before collecting\n');
		collected <- collect.info(input.info, weighted=weighted, covar.file='', flag=flag, delim='\t', header=header);
		cat('after collecting\n');
		gene.name <- collected[['gene.name']];
		input <- collected[['input']];
		snp.names <- collected[['snp.names']];
		snp.anno <- collected[['snp.anno']];
		group <- collected[['group']];
		truth <- collected[['truth']];
		lasso.file.path <- collected[['lasso.file.path']]
		top.file.path <- collected[['top.file.path']]
		within.file.path <- collected[['within.file.path']]
		cross.file.path <- collected[['cross.file.path']]
		locus.file.path <- collected[['locus.file.path']]
		comb.file.path <- collected[['comb.file.path']]
		effect.file.path <- collected[['effect.file.path']]
		log.file.path <- collected[['log.file.path']]	
		top <- collected[['top']];
		selected.within <- collected[['selected.within']]
		selected.cross <- collected[['selected.cross']]
		selected.locus <- collected[['selected.locus']]
		selected.comb <- collected[['selected.comb']]
		
		run.time <- c(time.lasso - time.pre);
		if(file.exists(comb.file.path) & skip==1) {
			cat(paste('skipping', comb.file.path, '\n')); flush.console();
			run.time <- 0;
		} else if(file.exists(lasso.file.path)) {		
			if(3 %in% steps) {
				top <- parse.lasso(input=input, snp.names=snp.names, group=group, lasso.file.path=lasso.file.path, top.file.path=top.file.path, weighted=weighted)
			}
			time.top <- Sys.time()
			
			## middle layer
			if(4 %in% steps) {
				selected.within <- select.within.groups(top=top, input=input, snp.names=snp.names, group=group, within.file.path=within.file.path)
			}
			time.within <- Sys.time()	
			
			## inner layer
			if(5 %in% steps) {
				selected.cross <- select.cross.groups(selected.within=selected.within, input=input, snp.names=snp.names, group=group, cross.file.path=cross.file.path)
			}
			time.cross <- Sys.time()	
			if(6 %in% steps) {
				selected.locus <- report.locus(gene.name=gene.name, selected.cross=selected.cross, input=input, snp.names=snp.names, group=group, rsq.cutoff=0.01, locus.file.path);
			}
			time.locus <- Sys.time()
			if(7 %in% steps) {
				selected.comb <- combine.locus(gene.name=gene.name, selected.locus=selected.locus, selected.cross=selected.cross, input=input, snp.names=snp.names, group=group, rsq.cutoff=0.01, comb.file.path=comb.file.path, snp.anno=snp.anno)
			}
			time.comb <- Sys.time()
			if(8 %in% steps) {
				effect <- estimate.effect.comb(gene.name=gene.name, selected.comb=selected.comb, input=input, effect.file.path=effect.file.path, p.cutoff=p.cutoff)
			}
			time.effect <- Sys.time()
			cat(paste(gene.name, 'successfully finished\n')); flush.console();
			
			run.time <- c(time.comb-time.0, time.comb-time.pre, time.pre-time.0, time.lasso-time.pre, time.top-time.lasso, time.within-time.top, time.cross-time.within, time.locus-time.cross, time.comb-time.locus);	
			if(9 %in% steps) {
				write.table(t(run.time), log.file.path, sep='\t', row.names=F, col.names=F, quote=F);
			}
		}
		return(run.time);
	}
}
# for(n in 1:200) { run.1(input.list.all[[i]][[n]], steps=4:7) }

retry <- function(input.info, weighted='', steps=1:8) {
	comb.file.path <- gsub('.RData', '.comb.txt', input.info[1, 'output.file']);
	if(!file.exists(comb.file.path)) {
		run.time <- run.1(input.info, weighted=weighted, steps=steps);	
		return(run.time);
	}
}

retry.steps <- function(input.info, time.baseline, weighted='', all.steps=F) {
	time.baseline <- strptime(time.baseline,"%Y-%m-%d %H:%M:%S");
	output.file <- input.info[1, 'output.file'];
	locus.file <- gsub('RData', 'txt', output.file);
	comb.file <- gsub('RData', 'comb.txt', output.file);
	output.file.time <- file.info(output.file)$mtime;
	locus.file.time <- file.info(locus.file)$mtime;
	comb.file.time <- file.info(comb.file)$mtime;
	
	if(all.steps) {
		if(comb.file.time < time.baseline) {
			run.time <- run.1(input.info, weighted=weighted, steps=1:6);	
			return(run.time);
		}
	}
	
	if(output.file.time < time.baseline) {
		run.time <- run.1(input.info, weighted=weighted, steps=3:5);	
		return(run.time);
	}

	if(locus.file.time < time.baseline) {
		run.time <- run.1(input.info, weighted=weighted, steps=4:5);	
		return(run.time);
	}

	if(comb.file.time < time.baseline) {
		run.time <- run.1(input.info, weighted=weighted, steps=5);	
		return(run.time);
	}
}

# i=5; n=184; input.list.this <- input.list.all[[i]]; input.info=finished.info=input.list.this[[n]]; input.folder='./'; input.info$gene.name; 
conditional.1 <- function(input.info, input.folder='') {
	stop.signal <- F;
	stop.signal <- tryCatch({
		pre.processing.1(input.info, ld=F);
	}, error = function(e) {
		stop.signal <- T;
		cat(paste(input.info$gene.name, ' preprocessing error\n'));
		return(stop.signal);
	})

	gene.name <- input.info[1, 'gene.name'];
	geno.file <- input.info[1, 'geno.file'];
	response.file <- input.info[1, 'response.file'];
	output.file <- input.info[1, 'output.file'];
	cat(paste('conditional analysis for', gene.name, '...\n', sep=' ')); flush.console();

	time.input <- Sys.time()
	fn <- paste(input.folder, geno.file, sep='')
	geno <- read.table(fn, sep=',', header=T, stringsAsFactors=F);	
	snp.names <- data.frame(snp=colnames(geno), idx=1:ncol(geno));
	response <- read.table(paste(input.folder, response.file, sep=''), sep=',', header=F, stringsAsFactors=F);	
	data.all <- data.frame(resp=response[, 1], geno);
	
	time.cond <- Sys.time()
	sig <- my_conditional(data.all);
	
	time.cor <- Sys.time()
	sig <- merge(snp.names, sig, by='snp', all.x=F, all.y=T);
	sig$cor <- 1;
	sig$rank <- 1;
	result <- c();
	for(s in 1:nrow(sig)) {
		this.sig <- sig[s, ]
		this.idx <- sig[s, 'idx'];
		this.cor <- cor(geno, geno[, this.idx])^2;
		this.cor <- data.frame(snp=row.names(this.cor), p=NA, p.w=NA, eff=NA, cor=this.cor);
		this.cor <- this.cor[which(this.cor$cor >= 0.5), ]
		this.cor <- merge(snp.names, this.cor, by='snp', all.x=F, all.y=T);
		this.cor <- this.cor[order(-this.cor$cor), ]
		this.cor <- this.cor[which(this.cor$idx != this.idx), ]
		this.sig.cor <- this.sig;
		if(nrow(this.cor) > 0) {
			this.cor$rank <- 2:(nrow(this.cor)+1)
			this.sig.cor <- rbind(this.sig, this.cor);
		}
		this.sig.cor$locus <- s;
		this.sig.cor$gene.name <- gene.name;
		result <- rbind(result, this.sig.cor);
	}
	if(!is.null(result) && nrow(result) > 0) {
		result <- result[, which(!colnames(result) %in% c('p.w', 'idx'))]
	}
	time.end <- Sys.time()
	run.time <- c(time.cond-time.input, time.cor-time.cond, time.end-time.cor, time.end-time.input);	
	
	write.table(format(result, digits=3), paste(input.folder, gsub('.RData', '.cond.out', output.file), sep=''), sep='\t', row.names=F, quote=F);
	write.table(t(run.time), paste(input.folder, gsub('.RData', '.cond.log', output.file), sep=''), sep='\t', row.names=F, col.names=F, quote=F);

	return(result);
}

#' Perform TreeMap analysis for fine mapping.
#' @param input.folder: The directory where input files are available. Default is the current working directory. Each input file is a tab-delimited text file. The first line has column headers with no specific restrictions. The remaining lines have sample ID in the 1st column, gene expression value in the 2nd column, and genotype data in the other columns. Genotypes values are coded as 0/1/2 for homozygous reference, heterozygous and homozygous alternative, respectively. Each row has data from one individual. A sample input file "SMNT.simulated.txt" is available in the "examples/" folder. 
#' @param pattern: A string that matches part of the input file names. This string shall be present in only input files, and absent from other file names in the same folder. Default is '.txt'.
#' @param output.folder: The directory where output files and intermediate files are to be saved. If the specified directory does not exist, it will be automatically created. Default is the current working directory.
#' @param steps: A vector of integers corresponding to the functions that will be executed. Default is 1:9. The final output is written to a file named xx.treemap.out where xx matches the name of the input file.
#'	1: pre-processing. This step removes variants with missing genotypes in >50% samples, and remove samples with missing genotypes of >50% variants. It performs z-transformation to scale each genotype column. It groups variants into a hierarchical structure based on linkage disequilibrium. This step produces several intermediate files including xx.response.csv.gz, xx.geno.csv.gz, xx.geno.norm.csv.gz and xx.group.csv.gz that correspond to gene expression values, original genotypes, normalized genotype values, LD groups, respectively.
#'	2: outer-layer. This step performs tree-guided group lasso. It produces an intermediate file xx.tree-lasso.csv
#'	3: outer-layer. This step parses lasso output. It produces an intermediate file xx.top.txt.
#'	4: middle-layer. This step selects variants within an LD block. It produces an intermediate file xx.within.RData.
#'	5: inner-layer. This step selects variants across multiple LD blocks and individual variants. It produces an intermediate file xx.cross.RData.
#'	6: inner-layer. This step produces the global optimal combination of multiple loci. It produces an intermediate file xx.locus.txt.
#'	7: inner-layer. This step ranks variants within each selected locus. It produces an intermediate file xx.comb.txt.
#'	8: estimate effect: This step esimates the effect size of each top-ranked variant at a locus via multivariate linear regression. It produces the final output file xx.treemap.out.
#'	9: log run-time: This step writes run time of each of the previous step in a log file xx.log.txt.
#' @param mc.cores: The number of parallel processes that will be forked for batch processing. This is useful for multi-core desktops or server clusters. Default is 1. 
#' @return NULL
#'	The final fine mapping result is written to the file xx.treemap.out. The "locus" column indicates unique causal locus. The "rank" column indicates ranking of variants within each locus. Variants with rank of 1 are lead eVars at a locus. Lead eVars are included in the global optimal solution, and their estimated effect sizes are in the "effect.est" column. Variants with rank >1 are suboptimal solutions that are linked to the lead eVar at each locus, with r2 values shown in the "cor" column.
treemap <- function(input.folder='./', pattern='.txt', output.folder='./', steps=1:9, mc.cores=1) {
	input.list <- create.batch.processing.list(input.folder.path=input.folder, pattern=pattern, output.folder.path=output.folder); 
	treemap.output <- mclapply(input.list, run.1, weighted='', mc.cores=mc.cores, steps=steps)
}

#' Perform simple stepwise conditional analysis for fine mapping.
#' @param input.folder: The directory where input files are available. Default is the current working directory.
#' @param pattern: A string that matches part of the input file names. This string shall be present in only input files, and absent from other file names in the same folder. Default is '.txt'.
#' @param output.folder: The directory where output files and intermediate files are to be saved. If the specified directory does not exist, it will be automatically created. Default is the current working directory.
#' @param mc.cores: The number of parallel processes that will be forked for batch processing. This is useful for multi-core desktops or server clusters. Default is 1. 
#' @return NULL
#'	The final fine mapping result is written to the file xx.treemap.out. The "locus" column indicates unique causal locus. The "rank" column indicates ranking of variants within each locus. Variants with rank of 1 are lead eVars at a locus. P-values and effects of lead eVars are in the "p" column and the "eff" column, respectively. Variants with rank >1 are ordered based on their linkage to the lead eVar at each locus, with r2 values shown in the "cor" column.
conditional <- function(input.folder='./', pattern='.txt', output.folder='./', mc.cores=1) {
	input.list <- create.batch.processing.list(input.folder.path=input.folder, pattern=pattern, output.folder.path=output.folder); 
	conditional.output <- mclapply(input.list, conditional.1, mc.cores=mc.cores)
}
