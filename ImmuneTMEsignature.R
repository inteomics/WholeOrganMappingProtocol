# Download and read in R “marker_Immune_TME.csv”, containing the Immune TME signature
# https://github.com/inteomics/WholeOrganMappingProtocol/tree/RNAseq
marker_Immune_TME<-read.csv("marker_Immune_TME.csv")  

## prepare for input data
# norm_expr_mat: gene expression/normalized RNAseq data matrix with p rows(genes) and n #columns (subjects);
exp.log2_norm.Immune_TME<- norm_expr_mat[marker_Immune_TME$gene,]

# Immune score function
immune.score.fun<-function(Exp.data) {
  Mi<-apply(Exp.data, 2,median)
  Grand_M<-mean( Mi )
  return( Mi-Grand_M)
}

# Calculate the immune score
ImmuScore.TME <- immune.score.fun(exp.log2_norm.Immune_TME)
