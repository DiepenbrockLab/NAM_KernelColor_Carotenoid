rm(list = ls())
##########################################################################

###Required files:
### (1) Tabular summary file with overlapping support intervals generated in script 3a
### (2) candidate gene list (see sample)
### (3) GWAS Results Summaries generated from script 2a, "Complete Results..." file

extractinate.GWAS.SNPs <- function(){
  #Function 1 should begin here
  
  #Obtain chromosome, start and stop support interval positions, from the "tabular summary"
  row.of.interest <- which((tabular.summary[,1] == trait)  & (tabular.summary[,6] == QTL))
  chr <- tabular.summary[row.of.interest, 2]
  QTL.start.bp <- as.numeric(tabular.summary[row.of.interest, 4])
  QTL.stop.bp <- as.numeric(tabular.summary[row.of.interest, 5])
  
  #To be rigorous, add an object that has all of the names of the NAM founders
  pop.seq <- as.data.frame(as.factor(c("pop01", "pop03", "pop08", 
                                       "pop10", "pop12", "pop13",
                                       "pop20", "pop21", "pop23",  
                                       "pop25")))
  founder.names <- as.data.frame(c("B97", "CML228", "CML52", 
                                  "HP301", "KI11","KI3", "NC350", "NC358", "OH7B","TX303"))
  NAM.pops <- cbind(pop.seq, founder.names)
  colnames(NAM.pops) <- c("Pop.num", "Pop.Founders")
  
  #Obtain the GRZM IDs of genes that are in the support interval
 
   cand.chr <- candidate.gene.list[which(candidate.gene.list[,4] == chr),] #CHd changed from 2 to 4 to fit data set, 11/4/2019
   cand.chr[,5] = as.numeric(cand.chr[,5])
   cand.chr[,6] = as.numeric(cand.chr[,6])

   cand.gene.names <- cand.chr[ which(( (cand.chr[,5] > QTL.start.bp) & (cand.chr[,6] < QTL.stop.bp) ) |  
             ( (cand.chr[,5] < QTL.start.bp) & (cand.chr[,6] > QTL.start.bp)  ) |
             ( (cand.chr[,5] < QTL.stop.bp) & (cand.chr[,6] > QTL.stop.bp) ) ), 3]
   
   genes.identified.in.SI <- NULL
   for(k in 1:length(cand.gene.names)) genes.identified.in.SI <- paste(genes.identified.in.SI,", ", cand.gene.names[k], sep = "")
  #Final product of this phase: genes.identified
 
  #Obtain the GWAS results for this trait, and filter out all SNPs with RMIP < cutoff specified at top of script (5 is default)
  #if (trait %in% multicollinear.traits){
  #GWAS.results <- read.table(paste(home.dir,location.of.GWAS.results,trait,"_no_collin\\Complete_Results_",trait,"_no_collin.txt", sep = ""),head = TRUE)} else {
  if(trait == "dT3"){
    GWAS.results <- read.table(paste(GWAS.results.path,"Complete_Results_",trait,"_Redone.txt", sep = ""),head = TRUE)#}
  }else{
   GWAS.results <- read.table(paste(GWAS.results.path,"Complete_Results_",trait,".txt", sep = ""),head = TRUE)#}
  }
  GWAS.results <- GWAS.results[-which(GWAS.results[,ncol(GWAS.results)] < RMIP_cutoff),]
  
  #If there are no SNPs within the JL support interval, STOP the script. Print out a message saying that JL and GWAS results do not overlap.
  GWAS.results.within.SI <- GWAS.results[which((GWAS.results[,1] == chr) & (GWAS.results[,2] > QTL.start.bp) & 
                                         (GWAS.results[,2] < QTL.stop.bp)),]
  
  number.of.SNPs.in.SI <- nrow(GWAS.results.within.SI)
  
  return(list(genes.identified.in.SI = genes.identified.in.SI, number.of.SNPs.in.SI = number.of.SNPs.in.SI))
}#end extractinate.GWAS.SNPs 

##########################################################################
##########################################################################
##########################################################################
#Set the working directory
#setwd("C:/Users/chdiep/Documents/Model_Fitting_2019_from.Server/")
setwd("C:/Users/mishi/OneDrive/Documents/R Projects/3c/")
home.dir <- getwd()
#geno.path = paste(home.dir, "/Geno.Pheno_Inputs/",sep='') #CHD commented 11/4/19, not used in this script
#trait <- c("aT","aT3","dT","dT3","gT","gT3","PC8","totalT","totalT3","totalTocochrs") 
trait <- c("a_carotene","b_carotene","b_cryp","lutein","phytofluene","zeaxanthin","zeino","total_carot","traitKC")
#trait = "dT3" #CHD added 19 May for re-do
#GWAS.results.path = "C:/Users/chdiep/Documents/Model_Fitting_2019_from.Server/Master_Files_forGWAS/GWAS_Results_fromCBSU_CNVsRemoved/Summaries_by_Trait/" #CHD commented 1/4/2020 to run with upliftable GWAS results
#GWAS.results.path = "C:/Users/chdiep/Documents/Model_Fitting_2019_from.Server/Master_Files_forGWAS/GWAS_Results_fromCBSU_CNVs.Removed_Upliftable/Summaries_by_Trait/"
#GWAS.results.path = "C:/Users/chdiep/Documents/Model_Fitting_2019_from.Server/Master_Files_forGWAS/GWAS_Results_fromCBSU_CNVs.Removed_Upliftable/Summaries_by_Trait/"
GWAS.results.path =  "C:/Users/mishi/Dropbox/NAM_KernelColor_HPLC/Results/GWAS_Results_forMishi/Summaries_by_Trait/"
trait.set <- "carot"
RMIP_cutoff = 5
#tabSummary.path = "C:/Users/chdiep/Documents/Model_Fitting_2019_from.Server/Summary_Tables.Figures/"
tabSummary.path = "C:/Users/mishi/OneDrive/Documents/R Projects/3c/"

#Read in the appropriate files
# summary with common support intervals
tabular.summary <- read.table(paste(tabSummary.path,"Tab_Sum_Final_carot_0904_with_Common_SI_Info_left_bound.txt", sep = ""), sep="\t",head = TRUE)

#Read in the candidate genes
candidate.gene.list <- read.table(paste(tabSummary.path,"Supplem.Data.Set.1_carot_candidate_genes DEAN.EDITS_DONE_forR.txt",sep=''), head = TRUE)

#Run this for loop
the.GRZM.results <- NULL
the.SNP.results <- NULL
for(i in 1:nrow(tabular.summary)){
  trait <- as.character(tabular.summary[i,1])
  if(trait == "Total_Tocopherols"){trait = "totalT"}
  if(trait == "Total_Tocotrienols"){trait = "totalT3"}
  if(trait == "Total_Tocochromanols"){trait = "totalTocochrs"}
  QTL <- as.character(tabular.summary[i,6])
  #For loop through all rows of the tabluar summary
  the.GWAS.SNP.numbers.and.GRZMs.in.SI <- extractinate.GWAS.SNPs()
  the.GRZM.results <- c(the.GRZM.results, the.GWAS.SNP.numbers.and.GRZMs.in.SI$genes.identified.in.SI)
  the.SNP.results <- c(the.SNP.results, the.GWAS.SNP.numbers.and.GRZMs.in.SI$number.of.SNPs.in.SI)
}
######################################################################
#Append the.GRZM.results and the.SNP.results to the tabular.summary
tabular.summary <- cbind(tabular.summary, the.GRZM.results, the.SNP.results)

#Output the updated tabular summary
write.table(tabular.summary, paste(tabSummary.path,"Tab_Sum_",trait.set,"_alpha.01_SI_with_Overlapping_GWAS_SNPs_2019_Carot_Upliftable.txt",sep=''), sep = "\t", row.names = FALSE, quote = FALSE) #CHD added _Upliftable suffix 1/4/2020

