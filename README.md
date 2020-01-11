## TreeMap - A Structured Approach to Fine Mapping of eQTL Variants
TreeMap prioritizes putative causal variants in cis-eQTL accounting for multisite effects and genetic linkage at a locus. 

It takes a 3-layer nested design to remove uninformative variants and reduce redundancies among informative variants progressively.  At the outer layer, tree-guided penalized regression selects groups of variants in LD or individual variants associated with transcriptional changes. At the middle layer, stepwise conditional analysis iterates combinations of variants within each LD block to identify a block-specific optimal solution. At the inner layer, variants selected from the previous layers are aggregated and passed through a Bayesian multivariate analysis to derive a global optimal solution. 

TreeMap supports parallel processing in cluster environment.

## Install. 
````
 ## if  you have not installed devtools, please do
 install.packages("devtools");
 ## then install treemap
 devtools::install_github('liliulab/treemap')
````

## Usage
```` 
library(treemap)

# perform treemap analysis. 
# Example: Please download the sample input file "SMNT.simulated.txt" from the "examples" folder and save it in your current working directory

treemap(input.folder='./', pattern='.simulated.txt', output.folder='./output/', steps=1:9, mc.cores=1)
 
# The treemap command creates a new "output/" folder in which a final output file "SMNT.treemap.out" is produced. 
# You can compare this output file with the "examples/output/SMNT.treemap.out" file on this github page.
# The "output/" folder also contains several intermediates files that can be reused during re-analysis to save time. 
# For re-analysis, users can specify "steps" to execute outer-, middle- or inner-layer functions. 
# Additional usages can be found via using "?treemap" after installing the R treemap package.
 
# perform simple stepwise conditional analysis. 

conditional(input.folder='./', pattern='.simulated.txt', output.folder='./output/', mc.cores=1)

# The conditional command creates a new "output/" folder in which a final output file "SMNT.cond.out" is produced. 
# Compare the output file "SMNT.cond.out" from the above command with the "examples/output/SMNT.cond.out" file on this github page.  
# This function is provided in case users prefer simple stepwise conditional analysis for fine mapping. 
# Additional usages can be found via using "?conditional" after installing the R treemap package.
````

The input file is a tab-delimited text file. The first line has column headers with no specific restrictions. The remaining lines have sample ID in the 1st column, gene expression value in the 2nd column, and genotype data in the other columns. Genotypes values are coded as 0/1/2 for homozygous reference, heterozygous and homozygous alternative, respectively. Each row has data from one individual. A sample input file "SMNT.simulated.txt" is available in the "examples/" folder. Below is a brief example.
ind_id	ILphe	v1	v2	v3	
HG00096	4.03	0	2	0	
HG00097	-1.24	0	2	0
HG00099	0.60	0	1	0	

## Reference
Manuscript is currently under peer-review.

Li Liu, Pramod Chandrashekar, Biao Zeng, Maxwell D. Sanderford, Sudhir Kumar, Greg Gibson (2019) TreeMap: A Structured Approach to Fine Mapping of eQTL Variants.  

## Contributors
Li Liu developed the algorithm of TreeMap. Pramod Chandrashekar implemented tree-guided group lasso in R based on matlab codes from the SLEP package. Please contact liliu at asu.edu for any questions or suggestions.