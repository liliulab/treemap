## TreeMap - A Structured Approach to Fine Mapping of eQTL Variants
TreeMap prioritizes putative causal variants in cis-eQTL accounting for multisite effects and genetic linkage at a locus. 

It takes a 3-layer nested design to remove uninformative variants and reduce redundancies among informative variants progressively.  At the outer layer, tree-guided penalized regression selects groups of variants in LD or individual variants associated with transcriptional changes. At the middle layer, stepwise conditional analysis iterates combinations of variants within each LD block to identify a block-specific optimal solution. At the inner layer, variants selected from the previous layers are aggregated and passed through a Bayesian multivariate analysis to derive a global optimal solution. 

Users can provide an optinal annotation file containing functional scores of variants, which will be used to weigh variants based on functional importance.

TreeMap supports parallel processing in cluster environment.

## Install. 
````
 devtools::install_github('liliulab/treemap')
````

## Usage
```` 
library(treemap)

# perform treemap analysis. 
# A sample input file "SMNT.exp_geno.txt" is available in the "examples" folder. 
treemap(input.folder='examples/', pattern='.exp_geno.txt', output.folder='examples/output/', steps=1:9, mc.cores=1) 
# Compare the output file "SMNT.treemap.out" from the above command with the "examples/SMNT.treemap.out" file in the "examples/" folder. 
# Several intermediates files are produced that can be reused during re-analysis to save time. 
# For re-analysis, users can specify "steps" to execute outer-, middle- or inner-layer functions. 
# Additional usages can be found via using "?treemap" after installing the R treemap package.

# perform treemap analysis with functional weighting. 
# set weighted='func'.
# A sample annotation file "SMNT.anno.txt" is available in the "examples" folder. 
treemap(input.folder='examples/', pattern='.exp_geno.txt', output.folder='examples/output/', steps=1:9, mc.cores=1, weighted='func')
 
# perform simple stepwise conditional analysis. 
conditional(input.folder='examples/', pattern='.simulated.txt', output.folder='examples/output/', mc.cores=1)
# Compare the output file "SMNT.cond.out" from the above command with the "examples/SMNT.cond.out" file in the "examples/" folder.  
# This function is provided in case users prefer simple stepwise conditional analysis for fine mapping. 
# Additional usages can be found via using "?conditional" after installing the R treemap package.

````

The input file is a tab-delimited text file. The first line has column headers with no specific restrictions. The remaining lines have sample ID in the 1st column, gene expression value in the 2nd column, and genotype data in the other columns. Genotypes values are coded as 0/1/2 for homozygous reference, heterozygous and homozygous alternative, respectively. Each row has data from one individual. A sample input file "SMNT.simulated.txt" is available in the "examples/" folder. Below is a brief example.<br>
ind_id	ILphe	v1	v2	v3<br>
HG00096	4.03	0	2	0<br>
HG00097	-1.24	0	2	0<br>
HG00099	0.60	0	1	0<br>

The optional annotation file is a tab-delimited text file. The first line has column headers with no specific restrictions. The remaining lines have variant IDs in the 1st column and functional scores in the 2nd column. Variants with high functional scores have increased chance of being selected as causal. A sample input file "SMNT.simulated.txt" is available in the "examples/" folder. Below is a brief example.<br>

## Reference

Li Liu, Pramod Chandrashekar, Biao Zeng, Maxwell D. Sanderford, Sudhir Kumar, Greg Gibson (2020) TreeMap: A Structured Approach to Fine Mapping of eQTL Variants. Bioinformatics 37(8):1125-1134 https://doi.org/10.1093/bioinformatics/btaa927 

## Contributors
Li Liu developed the algorithm of TreeMap. Pramod Chandrashekar implemented tree-guided group lasso in R based on matlab codes from the SLEP package. Please contact liliu at asu.edu for any questions or suggestions.