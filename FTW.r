#==========================================
# LDGH Project
#==========================================
#
# (/u0254) Copyleft 2017, by LDGH and Contributors.
#
# 
# -----------------
# FTW.r <Fstatiscs Table Workaround>
# -----------------
# GNU GPL 2017, by LDGH and Contributors.
#
# Original Author: Giordano Bruno Soares-Souza
# Contributor(s): 
# Updated by: Giordano Bruno Soares-Souza (04/09/2017)
#
# Command line: R CMD BATCH --no-save --no-restore ./FTW.r <plink_bim_file> <plink_freq_file> <hierfstat_output>
#				Rscript FTW.r <plink_bim_file> <plink_freq_file> <hierfstat_output>	
# Dependencies: R packages: data.table
# Description: Creates a data table encompassing mapping coordinates, F-statistics and cluster frequencies
# INPUTS:
# 1) <plink_bim_file> MAPPING FILE - PLINK BIM FORMAT
# 2) <plink_freq_file> FREQUENCY FILE - OUTPUT OF Plink --frq flag (usually with --within flag) 
# 3) <hierfstat_output> Fstatiscs obtained from HIERFSTAT script (eg. SNPID FST FIS FIT)
# FUTURE DEVELOPMENTS:
# 1) Graphical outputs
################################################################################

if(!require(data.table)) { install.packages("data.table"); require(data.table)}

# GET ARGUMENTS
bimFile = args[1]
freqFile = args[2]
fstFile = args[3]

# OUTPUT FILENAME
outFile = gsub(".out", ".complete.out", fstFile)

# READ DATA FROM FILES

bim = fread(bimFile,sep = "auto", header=FALSE, stringsAsFactors=FALSE,na.strings = "NA")
freq = fread(freqFile,sep = "auto", header=TRUE, stringsAsFactors=FALSE,na.strings = "NA")
fstatiscs = fread(fstFile,sep = "auto", header=TRUE, stringsAsFactors=FALSE,na.strings = "NA")

colnames(bim) = c('CHR','RS','cM','POS','A1','A2')

# Transform Frequency file from Long to Wide by cluster (CLST) column
freqLong = dcast(freq,SNP+A1+A2~CLST,value.var='MAF')

# MERGE FILES
mergeBIMFST = merge(bim[,c(1,2,4,3)],fstatiscs,by.x="V2",by.y="POLYM",sort=F)
mergeBIMFSTFRQ = merge(mergeBIMFST,freqLong,by.x="RS",by.y="SNP",sort=F)

# PRINT OUTPUT
write.table(mergeBIMFSTFRQ,outFile,row.names=F,col.names=T,quote=F)
