# Download and read in R “EMT_score_w.csv”, containing the EMT signature and their #corresponding correlation coefficients
# https://github.com/inteomics/WholeOrganMappingProtocol/tree/RNAseq
EMT_w<-read.csv("EMT_score_w.csv") 

## prepare for input data
# norm_expr_mat: gene expression/normalized RNAseq data matrix with p rows(genes) and n #columns (subjects);
EMT_mark_mat_scale <- scale(t(norm_expr_mat[EMT_w$Gene.Symbol , ]), scale=T)

## calculate the EMT score
EMT_score = as.vector(EMT_mark_mat_scale %*% EMT_w$score)
