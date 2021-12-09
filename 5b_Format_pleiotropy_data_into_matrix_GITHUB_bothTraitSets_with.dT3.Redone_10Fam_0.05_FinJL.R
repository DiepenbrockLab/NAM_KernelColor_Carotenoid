rm(list = ls())
#saved in CBK files as: Create_Pleitropic_QTL_Matrix_for_SI01_data_cbk.r


###Required files:
### (1) JL results in tabular summary - current version organized by QTL left bound, preferably with common.SI tracker 
### (2) ",trait.1,".QTL","that.are.pleiotropic.with", trait.2,".csv", which contains the correlation of original and refit effect estimates


###Output generated:
### (1) Pleiotropic output matrix, which extracts correlations from above files and stores in matrix form (to be used as input file for graph generation)




#setwd("C:\\Users\\ceb19\\Documents\\Gore Lab\\Carotenoid NAM Merged Env")
#pleio.data.dir <- "\\(14)Pleiotropy\\Pleiotropy_Analysis_SI01\\Pleio_Investigation_All_Sig_Test\\"
#setwd("C:/Users/chd45/Documents/Projects/NAM_GWAS/CHD_Tassel3fromSF_modified0.01/")
#home.dir <- getwd()
#tabSummary.path = "/Tabular_Summaries/"

setwd("C:/Users/chdiep/Documents/NAM_KernelColor_Carot/")
home.dir <- getwd()
tabSummary.path = "/Summary_Tables.Figures/"
pleio.dir = "/Summary_Tables.Figures/Pleiotropy/"
pleio.data.dir = "/Summary_Tables.Figures/Pleiotropy/Signif_Pleio_Pairwise_Lists/"
pleio.matrix.dir <- "/Summary_Tables.Figures/Pleiotropy/Matrix/"
trait.set = 'carot'

threshold=0.05

library(gplots)


#Read in the "tabular summary" of JL analysis results
#tab.summary <-read.table(paste(home.dir,"\\(16)Generating Robust Files for Group Review\\Generating overlapping support intervals\\Tab_Sum_Carot_alpha_0.01_20140603_for_Pleio_Script_Comps_Only.txt", sep = ""), head = TRUE) 
 
tab.summary <- read.table(paste(home.dir,tabSummary.path,"Tabular_Summary_of_JL_carot_Results_for_all_traits_SI01_newPVEv2_corrsigns_FinJL.txt",sep=''), sep='\t', head = TRUE)
#tab.summary_carot = read.table(paste(home.dir,tabSummary.path,"Tab_Sum_Carot_alpha.01_GWAS_FamPVE_common_SI_recsuppregions_LODscores_20150612.txt",sep=''), head = TRUE)
#tab.summary = rbind(tab.summary_toco, tab.summary_carot)
#traits <- c("ACAR_RUV", "BCAR_RUV", "BCRY_RUV", "LUT_RUV", "PHYF_RUV", "THLYC_RUV", "TOTCAR_RUV", "ZEA_RUV", "ZEI_RUV", "aT","aT3","dT","dT3","gT","gT3","PC8","totalT","totalT3","totalTocochrs")

#create an initial output matrix. it will be nrow(tab summary)xnrow(tab summary). Populate it with NAs
pleio.output.matrix <- matrix(NA, nrow = nrow(tab.summary), ncol = length(unique(tab.summary[,1])))

number.real.files = 0
for(i in 1:nrow(tab.summary)){ #for loop through i (1 to nrows(tab summary))
  for(j in as.character(unique(tab.summary[,1]))){ #for loop through j in 1:length(unique(tab.summary[,1])))
    correlation <- NA
    print(paste("i = ", i, sep = ""))
    print(paste("j = ", j, sep = ""))
    
    #trait.a will be ith row 1st column
    trait.a <- as.character(tab.summary[i,1])

    #trait.b will be jth row 1st column
    trait.b <- j
    
    #peak.SNP will be ith row 6th column
    peak.SNP <- as.character(tab.summary[i,6])
    
    #if i == j, set the corresponding element equal to 1 and skip to next row     
    if(trait.a == trait.b){
      pleio.output.matrix[i,which(as.character(unique(tab.summary[,1]))== j)] <- 1
      next
    }
    
    #read in "trait.a.QTLthat.are.pleiotropic.withtrait.b.csv" CHD changed to txt in script 7-1, also added strings as factors FALSE because correlations otherwise read as factors (meaningless whole numbers)
    QTL.test.table <- try(read.table(paste(home.dir, pleio.data.dir,trait.a,".QTLthat.are.SIGNIF.pleiotropic.with",trait.b,"_",threshold,"_sign.corrected_FinJL.txt",sep = ""),stringsAsFactors=FALSE))
    
    #Do the following only if the above file exists
    if(class(QTL.test.table)=="try-error"){
      next
    }else{
    number.real.files = number.real.files + 1
    correlation <- as.numeric(QTL.test.table[which(as.character(QTL.test.table[,1]) == peak.SNP),4])
    print(correlation)
    #If such one exists, put the correlation coefficient (column 4 of "trait.a.QTLthat.are.pleiotropic.withtrait.b.csv") 
          #   into the i,jth element of the matrix
    if(length(correlation)>0){pleio.output.matrix[i,which(as.character(unique(tab.summary[,1]))== j)] <- correlation}
    }
  }#end for(j in 1:length(unique(tab.summary[,1])))
}#end for(i in 1:nrows(tab summary))


#put rows and column names on pleio.output.matrix
#output the matrix as a text file
rownames(pleio.output.matrix) <- paste(tab.summary[,1],"_Chr_", tab.summary[,2],"_Pos_",tab.summary[,3],"_Marker_", 
                                       tab.summary[,6], sep = "")

colnames(pleio.output.matrix) <- as.character(unique(tab.summary[,1]))


write.table(pleio.output.matrix, paste(home.dir,pleio.matrix.dir,"Pleiotropic.Output.Matrix.for.",trait.set,".SI01_with.dT3.redone_0.05_sign.corrected_FinJL.txt", sep = ""), 
                          sep = "\t", row.names = TRUE,col.names = TRUE,quote = FALSE)


#generate merged file containing tab summary and pleiotropy matrix
#to be used to generate pleiotropy network graphs
merge.file <- cbind(tab.summary, pleio.output.matrix)

write.table(merge.file, paste(home.dir,pleio.matrix.dir,"Pleiotropic.Output.Matrix.for.",trait.set,".SI01.TAB.SUM.MERGE_with.dT3.redone_0.05_sign.corrected_FinJL.txt", sep = ""), 
                          sep = "\t", row.names = FALSE,col.names = TRUE,quote = FALSE)
