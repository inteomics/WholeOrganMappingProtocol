#################################################################################
#Use nonparametric methods to analyze the monotonic changes in genes' DNA methylation
#Perform SAM analysis: V. Tusher, R. Tibshirani, and G. Chu, PNAS, 98:5116–5121, 2001.
#Created by Ziqiao Wang, PhD (email: wzqjanet@gmail.com)
#################################################################################
load("meth_betavalue_genelevel_all.RData")

#map 19
beta_map19=beta_genelevel[,which(info_all$map=="MAP19" | info_all$group=="control")]
info_map19=info_all[which(info_all$map=="MAP19" | info_all$group=="control"),]
table(info_map19$group)
#control      HG      LG      UC 
#8       3      27       6 
#T test for LG-control, HG-control, UC-control
ttest.list.lg=list()
ttest.list.hg=list()
ttest.list.uc=list()
for(i in 1:nrow(beta_map19)){
  ttest.list.lg[[i]] = t.test(as.numeric(beta_map19[i,which(info_map19$group=="LG" | info_map19$group=="control")])~info_map19$group[which(info_map19$group=="LG" | info_map19$group=="control")])
  ttest.list.hg[[i]] = t.test(as.numeric(beta_map19[i,which(info_map19$group=="HG" | info_map19$group=="control")])~info_map19$group[which(info_map19$group=="HG" | info_map19$group=="control")])
  ttest.list.uc[[i]] = t.test(as.numeric(beta_map19[i,which(info_map19$group=="UC" | info_map19$group=="control")])~info_map19$group[which(info_map19$group=="UC" | info_map19$group=="control")])
  
}

gene_NUs = beta_genelevel_rmdup[,which(info_all_rmdup$group=="control")]
gene_maps=beta_genelevel_rmdup[,-which(info_all_rmdup$group=="control")]
info_maps=info_all_rmdup[-which(info_all_rmdup$group=="control"),]

gene_maps_map19=gene_maps[,which(info_maps$map=="MAP19")]
info_maps19=info_maps[which(info_maps$map=="MAP19"),]
identical(as.character(info_maps19$ID),colnames(gene_maps_map19))
#[1] TRUE

#now conduct lg-control, hg-control, uc-control 3 sams
#LG-Control
library(samr)

gene_NUs = beta_genelevel[,which(info_all$group=="control")]
gene_maps=beta_genelevel[,-which(info_all$group=="control")]
info_maps=info_all[-which(info_all$group=="control"),]

gene_maps_map19=gene_maps[,which(info_maps$map=="MAP19")]
info_maps19=info_maps[which(info_maps$map=="MAP19"),]
identical(as.character(info_maps19$ID),colnames(gene_maps_map19))
#[1] TRUE

sam_input=cbind(gene_maps_map19[,which(info_maps19$group=="LG")],gene_NUs)
info_sam=c(rep(1,27),rep(2,8))
set.seed(10)
samfit_lg_map19=SAM(x=sam_input,y=info_sam,resp.type = "Two class unpaired",fdr.output = 1,random.seed = 10,genenames = rownames(sam_input))
save(samfit_lg_map19,file="samfit_lg_map19.RData")
fdr_lg_map19=samfit_lg_map19$siggenes.table
save(fdr_lg_map19,file="fdr_lg_map19.RData")

#HG-control
sam_input=cbind(gene_maps_map19[,which(info_maps19$group=="HG")],gene_NUs)
info_sam=c(rep(1,3),rep(2,8))
set.seed(10)
samfit_hg_map19=SAM(x=sam_input,y=info_sam,resp.type = "Two class unpaired",fdr.output = 1,random.seed = 10,genenames = rownames(sam_input))
save(samfit_hg_map19,file="samfit_hg_map19.RData")
fdr_hg_map19=samfit_hg_map19$siggenes.table
save(fdr_hg_map19,file="fdr_hg_map19.RData")

#UC-control
sam_input=cbind(gene_maps_map19[,which(info_maps19$group=="UC")],gene_NUs)
info_sam=c(rep(1,6),rep(2,8))
set.seed(10)
samfit_uc_map19=SAM(x=sam_input,y=info_sam,resp.type = "Two class unpaired",fdr.output = 1,random.seed = 10,genenames = rownames(sam_input))
save(samfit_uc_map19,file="samfit_uc_map19_all.RData")
fdr_uc_map19=samfit_uc_map19$siggenes.table
save(fdr_uc_map19,file="fdr_uc_map19_all.RData")

#################################################################################
#Summarize the SAM results
#################################################################################

sig_genes=function(fdr_sam_list,threshold_up=0.05,threshold_lo=0.05){
  fdr_sam_up=data.frame(fdr_sam_list$genes.up)
  fdr_sam_up$Gene.ID=as.character(fdr_sam_up$Gene.ID)
  fdr_sam_up$Score.d.=as.numeric(as.character(fdr_sam_up$Score.d.))
  fdr_sam_up$Fold.Change=as.numeric(as.character(fdr_sam_up$Fold.Change))
  fdr_sam_up$q.value...=as.numeric(as.character(fdr_sam_up$q.value...))
  
  
  fdr_sam_lo=data.frame(fdr_sam_list$genes.lo)
  fdr_sam_lo$Gene.ID=as.character(fdr_sam_lo$Gene.ID)
  fdr_sam_lo$Score.d.=as.numeric(as.character(fdr_sam_lo$Score.d.))
  fdr_sam_lo$Fold.Change=as.numeric(as.character(fdr_sam_lo$Fold.Change))
  fdr_sam_lo$q.value...=as.numeric(as.character(fdr_sam_lo$q.value...))
  
  
  id_up=which(fdr_sam_up$q.value...<=threshold_up)
  id_lo=which(fdr_sam_lo$q.value...<=threshold_lo)
  
  sig_genes=list(fdr_sam_lo[id_lo,],fdr_sam_up[id_up,])
  names(sig_genes)=c("sig_genes_lo","sig_genes_up")
  return(sig_genes)
}

load("fdr_hg_map19.RData")
load("fdr_lg_map19.RData")

lg_map19=sig_genes(fdr_lg_map19)
hg_map19=sig_genes(fdr_hg_map19)
uc_map19=sig_genes(fdr_uc_map19)

lg_map19_all=sig_genes(fdr_lg_map19,threshold_up=1,threshold_lo=1)
hg_map19_all=sig_genes(fdr_hg_map19,threshold_up=1,threshold_lo=1)
uc_map19_all=sig_genes(fdr_uc_map19,threshold_up=1,threshold_lo=1)

lg_map19=sig_genes(fdr_lg_map19,threshold_up=0.2,threshold_lo=0.2)
hg_map19=sig_genes(fdr_hg_map19,threshold_up=0.2,threshold_lo=0.2)
uc_map19=sig_genes(fdr_uc_map19,threshold_up=0.2,threshold_lo=0.2)

#venn diagram
HGIN <- c(hg_map19$sig_genes_lo$Gene.ID,hg_map19$sig_genes_up$Gene.ID)
LGIN <- c(lg_map19$sig_genes_lo$Gene.ID,lg_map19$sig_genes_up$Gene.ID)
UC <- c(uc_map19$sig_genes_lo$Gene.ID,uc_map19$sig_genes_up$Gene.ID)
# Prepare a palette of 3 colors with R colorbrewer:
library(RColorBrewer)
myCol <- brewer.pal(3, "Pastel2")
# Chart
venn.diagram(
  x = list(HGIN, LGIN, UC),
  category.names = c("HGIN (320)" , "LGIN (5899)" , "UC (2011)"),
  filename = 'map19_siggenes_venndiagram_SAM.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.4,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
  
)

# Chart
venn.diagram(
  x = list(HGIN, LGIN, UC),
  category.names = c("HGIN (520)" , "LGIN (6755)" , "UC (2474)"),
  filename = 'map19_siggenes_venndiagram_SAM_0_2.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.4,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
  
)

#LG & HG & UC group
lg_up=intersect(lg_map19$sig_genes_up$Gene.ID,hg_map19$sig_genes_up$Gene.ID)
lg_up=intersect(lg_up,uc_map19$sig_genes_up$Gene.ID) #72
lg_lo=intersect(lg_map19$sig_genes_lo$Gene.ID,hg_map19$sig_genes_lo$Gene.ID)
lg_lo=intersect(lg_lo,uc_map19$sig_genes_lo$Gene.ID) #178
n1_up=length(lg_up);n1_lo=length(lg_lo)
# HG & UC group and not in LG
hg_uc_up=intersect(hg_map19$sig_genes_up$Gene.ID,uc_map19$sig_genes_up$Gene.ID)
hg_uc_up=hg_uc_up[which(!(hg_uc_up %in% lg_map19$sig_genes_up$Gene.ID))] #1
hg_uc_lo=intersect(hg_map19$sig_genes_lo$Gene.ID,uc_map19$sig_genes_lo$Gene.ID)
hg_uc_lo=hg_uc_lo[which(!(hg_uc_lo %in% lg_map19$sig_genes_lo$Gene.ID))] #9
n2_up=length(hg_uc_up);n2_lo=length(hg_uc_lo)

#UC group alone
uc_up=uc_map19$sig_genes_up$Gene.ID
uc_up=uc_up[which(!(uc_up %in% lg_map19$sig_genes_up$Gene.ID))]
uc_up=uc_up[which(!(uc_up %in% hg_map19$sig_genes_up$Gene.ID))] #40

uc_lo=uc_map19$sig_genes_lo$Gene.ID
uc_lo=uc_lo[which(!(uc_lo %in% lg_map19$sig_genes_lo$Gene.ID))]
uc_lo=uc_lo[which(!(uc_lo %in% hg_map19$sig_genes_lo$Gene.ID))] #155
n3_up=length(uc_up);n3_lo=length(uc_lo)

gene_map19_threshold0_05=c(lg_lo,hg_uc_lo,uc_lo,uc_up,hg_uc_up,lg_up) #1658
group=c(rep("LG",n1_lo),rep("HG&UC",n2_lo),rep("UC",n3_lo),rep("UC",n3_up),rep("HG&UC",n2_up),rep("LG",n1_up))
direction=c(rep("hyper",n1_lo+n2_lo+n3_lo),rep("hypo",n1_up+n2_up+n3_up)) #notice that the test return a reverse results of group and control group, thus need to reverse the directions
gene_map19_threshold0_05_final=data.frame(cbind(gene_map19_threshold0_05,group,direction))

gene_map19_threshold0_2=c(lg_lo,hg_uc_lo,uc_lo,uc_up,hg_uc_up,lg_up) #1658
group=c(rep("LG",n1_lo),rep("HG&UC",n2_lo),rep("UC",n3_lo),rep("UC",n3_up),rep("HG&UC",n2_up),rep("LG",n1_up))
direction=c(rep("hyper",n1_lo+n2_lo+n3_lo),rep("hypo",n1_up+n2_up+n3_up)) #notice that the test return a reverse results of group and control group, thus need to reverse the directions
gene_map19_threshold0_2_final=data.frame(cbind(gene_map19_threshold0_2,group,direction))


