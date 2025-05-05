# Download and load in R “BLTscore.RData”, the linear discriminant analysis (LDA) model from
# https://github.com/inteomics/WholeOrganMappingProtocol/tree/RNAseq
load("BLTscore.RData")
ls()
#[1] "BLTscoreweights"

## prepare for input data;
#norm_expr_mat: gene expression/normalized RNAseq data matrix with p rows(genes) and n #columns (subjects);
#center and standardize gene expression data before multiplying the weights of 28 genes;
lda_rnaseq_scale <- scale(t(norm_expr_mat[names(BLTscoreweights), ]), scale=T)

## calculate BLT score
BLT_score = lda_rnaseq_scale %*% BLTscoreweight
