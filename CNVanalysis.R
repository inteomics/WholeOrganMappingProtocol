### Step1 (bash script)
# The log R ratios and B-allele frequencies were exported from GenomeStudio as a single table in wide format and then split into sample-specific files using a kcolumn script, which is part of the PennCNV software 
PennCNV-1.0.5/kcolumn.pl GenomeStudio_FinalReport_Wide_v1a.txt split 3 -heading 3 -tab -out MDA_v1 --name_by_header



### Step2 (bash script)
# Adjust the order of the columns in each of the files created in Step1 to match OncoSNP requirements
awk '$2!="0" {print $1,$2,$3,$5,$6}' FS="\t" OFS="\t" MDA_v1.6158600021_R01C01 | awk '!seen[$3]++' > MDA_v1.6158600021_R01C01.prc



### Step3 (bash script)
# Detection of copy-number variable regions conducted using the OncoSNP software (version 1.4)
# OncoSNP is executed for tumor-normal tissue sample pairs defined in the batch file, with both stromal contamination and intratumoral heterogeneity models
# The batch file contains 3 tab-delimited columns: SampleID, TumourFile, NormalFile; e.g.:
#   SampleID        TumourFile               NormalFile
#   Map18_F14       6158600029_R05C01.prc    6158600029_R02C01.prc
#
# The log R ratios and B-allele frequency plots are created using OncoSNP for all chromosomes
OncoSNP_v1.4/run_oncosnp.sh /usr/local/MATLAB/MATLAB_Compiler_Runtime/v82 
	--batch-file batch_file 
	--output-dir output_dir 
	--gcdir OncoSNP_v1.4/b37/ 
	--paramsfile OncoSNP_v1.4/configuration/hyperparameters.dat 
	--levelsfile OncoSNP_v1.4/configuration/levels-610.dat 
	--trainingstatesfile OncoSNP_v1.4/configuration/trainingStates.dat 
	--tumourstatesfile OncoSNP_v1.4/configuration/tumourStates.dat 
	--hgtables OncoSNP_v1.4/configuration/hgTables_b37.txt 
	--logfile log_file



### Step4 (R script)
# Process the OncoSNP CNV data using R to extract 
# The script uses *.cnvs and *.qc files created in Step3

library('openxlsx')

# Input data
files <- list.files(path = "output_dir", pattern = ".cnvs$")

# Parameters
cols <- c('Chromosome', 'StartPosition', 'EndPosition', 'CopyNumber', 'LOH', 'Rank', 'Loglik', 
          'nProbes', 'NormalFraction', 'TumourState', 'PloidyNo', 'MajorCopyNumber', 'MinorCopyNumber')
colsqcdat <- c('LogRRatioShift', 'NormalContent', 'Copy Number (Average)', 'Log-likelihood', 
               'OutlierRate', 'LogRRatioStd', 'BAlleleFreqStd', 'PloidyNo')
loglik_thresh <- 0.05

# Initialize result dataframe
CNVdata <- data.frame()

# Process files
for (file in files) {
  sample_name <- gsub('.cnvs$', '', file)
  
  # Read CNV and QC data
  tmpCNV <- read.table(paste0(directory, file), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(tmpCNV) <- cols
  tmpCNV$sample <- sample_name
  
  qcdat <- read.table(paste0(directory, gsub('.cnvs$', '.qc', file)), header = FALSE, sep = "\t")
  colnames(qcdat) <- colsqcdat
  
  # Determine ploidy selection
  relative_loglik_change <- abs(qcdat$`Log-likelihood`[2] - qcdat$`Log-likelihood`[1]) / qcdat$`Log-likelihood`[1]
  ploidySelect <- if (relative_loglik_change < loglik_thresh &&
                      qcdat$`Copy Number (Average)`[2] < qcdat$`Copy Number (Average)`[1]) {
    message(sprintf("%s: %s - Using PloidyNo 2 values. Average CN: %.2f", 
                    which(files == file), file, qcdat$`Copy Number (Average)`[2]))
    2
  } else {
    message(sprintf("%s: %s - Using PloidyNo 1 values. Average CN: %.2f", 
                    which(files == file), file, qcdat$`Copy Number (Average)`[1]))
    1
  }
  
  # Filter and append data
  CNVdata <- rbind(CNVdata, tmpCNV[tmpCNV$PloidyNo == ploidySelect & tmpCNV$Rank == 1, ])
}

# Sort results
CNVdata <- CNVdata[order(CNVdata$sample, CNVdata$Chromosome, CNVdata$StartPosition, CNVdata$Rank), ]

# Save to Excel
wb <- createWorkbook()
addWorksheet(wb, "CNV", gridLines = TRUE)
writeDataTable(wb, sheet = 1, x = CNVdata[, cols], withFilter = FALSE)
saveWorkbook(wb, 'OncoSNP_v1a.xlsx', overwrite = TRUE)



### Step5
# Additional visualizations of copy number-altered regions based on .cnvs and .qc result files provided by OncoSNP 

