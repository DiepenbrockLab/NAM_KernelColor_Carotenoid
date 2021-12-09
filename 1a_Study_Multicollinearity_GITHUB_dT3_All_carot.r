rm(list = ls())

###########NOTE this has been changed from original on 04/01/15 to only contain correlation analysis for JL markers
### Follow this analysis up by selecting which markers appear to be MC and excluding them from marker vector which will be used in TASSEL rescan function
###Required files:
### (1) imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt (GBS SNPs)
### (2) transformed BLUEs
### (3) JL model output, with NA in pop term row (modified version) and residual row removed from file

setwd("C:/Users/chdiep/Documents/Model_Fitting_2019_from.Server/") #changed CHD 5-11
home.dir <- getwd()
geno.path = home.dir
#trait <- c("aT","aT3","dT","dT3","gT","gT3","PC8","totalT","totalT3","totalTocochrs") #added CHD 5-11
trait <- c("a_carotene","b_carotene","b_cryp","lutein","phytofluene","zeaxanthin","zeino","total_carot")

#trait = "dT3" #CHD added 19 May for re-do
#dT3.new.path = paste(home.dir,"/Methods/dT3_removeExtremeVal_test/new.trans_new.perm_FINAL/",sep='')
#correlation.path = dT3.new.path
BLUE.or.BLUP <- "BLUE"  #Options are "BLUE" and "BLUP"
testing.correlation <- TRUE
transformed = TRUE

if(transformed == FALSE){
  #For JL on untransformed BLUEs, only for backtransformation validation
  JL.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/JL_Analysis_Scripts.Results/UNTRANSF_JLResults/",sep='')
  pheno.path = "/Users/anybody/Documents/Toco_NAM/Geno.Pheno_Inputs/BLUES_Untransf/"
  #pheno.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/Validation_Backtransformation/BLUES_Untransf/",sep='')
}else{
#For master JL analysis
  #JL.path = paste(home.dir,"/validate_CBK.AEL/CHD_Tassel3fromSF_modified0.01/Final_Models_forMultiCollCheck/",sep='')
  JL.path = paste(home.dir,"/Master_Files_forJL/",sep='')
  correlation.path = JL.path
  #pheno.path = paste(home.dir, "/Phenotypes/",sep ='')
  pheno.path = home.dir
}

#library(multtest)

#Read in the gneotypic data
#setwd(paste(home.dir, "\\validate_CBK.AEL", sep = "")) #changed CHD 5-11
genotypes <- read.table(paste(geno.path,"/imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt",sep=''), head = TRUE)

#Eventually, loop around the traits. This loop will begin here

for(i in trait){

#Read in the trait under study
pheno.data = read.table(paste(pheno.path,"/Master_Files_forJL/BLUEs_No_Outliers_Transformed_all_20190904_",i,".txt", sep = ""), head = TRUE)
  
  head(pheno.data)

#Read in the results from JL analysis
if(transformed == FALSE){
  TASSEL.model.results = read.table(paste(JL.path,i,"_model_3fromSF_0.01_UNTRANSF_Final_forR.txt",sep=''),head=TRUE)
}else{
    TASSEL.model.results <- read.table(paste(JL.path,i,"/",i,"_anova_for.R.txt", sep = "") , head = TRUE,stringsAsFactors = FALSE)
  }

if(is.numeric(TASSEL.model.results[,3]) == FALSE){
 print(paste("The chromosome column in ", i," model results needs to be numeric. Please change this before proceeding.", sep = ""))
 break;
}

#### #########################################################################################################################################
#Get the appropriate SNPs and merge the phenotypic and genotypic data together.
  
#Parse out all of the SNPs that were selected in the model
  
## CHD tested addition of chr number 8/27/19, but not needed given that there are no redundant SNP pos across chrs in the 14k map (and we care about correls, no matter the chr.)
#geno.names.to.add.chr = genotypes[,1]
#geno.chrs = genotypes[,3]
#geno.names.with.chr = rep(NA,length(geno.names.to.add.chr))
#for(m in 1:length(geno.names.to.add.chr)){
#  name.with.chr = paste("S",geno.chrs[m],"_",substr(geno.names.to.add.chr[m],start=3,stop=20),sep='')
#  geno.names.with.chr[m] = name.with.chr 
#}
  
geno.reduced <- genotypes[which(genotypes[,1] %in% TASSEL.model.results[,2]),-c(2:5)]
geno.reduced.formatted <-as.data.frame(t(geno.reduced[,-1]))
colnames(geno.reduced.formatted) <- as.character(t(geno.reduced[,1]))

#pheno.data will always have more data becuase IBM is included in the phenotypic data.
colnames(pheno.data)[1] = "Geno_Code" #added by CHD 5-12 so that pheno data has header script is looking for (otherwise is "X.trait.")
geno.and.pheno <- merge(pheno.data, geno.reduced.formatted, by.x = "Geno_Code", by.y = "row.names")

#Add a population column
geno.and.pheno <- cbind(geno.and.pheno[,1], as.factor(substr(as.character(geno.and.pheno[,1]), start = 3, stop = 4)), geno.and.pheno[,c(2:ncol(geno.and.pheno))])
colnames(geno.and.pheno)[2] <- "pop"

##############################################################################################################################################
#Calculate a correlation matrix between correlated SNPs, and print them out to a data set.
if(testing.correlation == TRUE){
  #uncomment two lines below if testing correlation across all significant markers
  TASSEL.model.results.2 <- TASSEL.model.results[3:nrow(TASSEL.model.results),]
  SNP.IDs <- as.vector(unique(TASSEL.model.results.2[,2]))
  Description.of.SNPs = "allSignifSNPs" #line added by CHD 5-12
  
  SNP.set <- geno.and.pheno[,which(colnames(geno.and.pheno) %in% SNP.IDs)]
  correlation.matrix <- cor(SNP.set, use = "pairwise.complete.obs")
  write.table(correlation.matrix, paste(correlation.path,"/Correlation_Matrices/Correlation.matrix_",Description.of.SNPs,"_",i,".txt", sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)
}
}# End for(i in trait)