rm(list = ls())

###Required files:
### (1) imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt (GBS SNPs)
### (2) transformed BLUEs
### (3) JL model output, with NA in pop term row (modified version) and residual row removed from file

#Use the multtest library to obtain FDR-adjusted P-values

setwd("/Users/MaryLaPorte/Documents/Model_Fitting_2019_from.Server")
home.dir <- getwd()
geno.path = home.dir
pheno.path = paste(home.dir, "/Master_Files_forJL/",sep='')
JL.path = paste(home.dir,"/Master_Files_forJL/",sep='')
popByMarker.path = paste(home.dir,"/Master_Files_forJL/popByMarker/",sep='')
residuals.path = paste(home.dir,"/Master_Files_forJL/Residuals.from.R/",sep='')
library(multtest)

trait <- c('total_carot')

#Read in the genotypic data
genotypes <- read.table(paste(geno.path,"/imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt",sep=''), head = TRUE)

pheno.data = read.table(paste(pheno.path,"BLUEs_No_Outliers_KernelColor.Carot_Transformed_all_forR.txt", sep = ""), head = TRUE)
head(pheno.data)

for(i in trait){

#Read in the results from JL analysis
#setwd(paste(home.dir, "\\Data_for_1pct_Support_Intervals\\JL_output\\JL_model_modified_SI01\\", sep = ""))
#new path for SI01 required here (CBK note)
  TASSEL.model.results <- read.table(paste(JL.path,i,"_JLProcessedData_anova_forR.txt", sep = "") , head = TRUE,stringsAsFactors = FALSE)

#############################################################################################################################################
#Get the appropriate SNPs and merge the phenotypic and genotypic data together.

#Parse out all of the SNPs that were selected in the model
geno.reduced <- genotypes[which(genotypes[,1] %in% TASSEL.model.results[2,]),-c(2:5)]
geno.reduced.formatted <-as.data.frame(t(geno.reduced[,-1]))
colnames(geno.reduced.formatted) <- as.character(t(geno.reduced[,1]))

#pheno.data will always have more data becuase IBM is included in the phenotypic data.
colnames(pheno.data)[1] = "Geno_Code"
geno.and.pheno <- merge(pheno.data, geno.reduced.formatted, by.x = "Geno_Code", by.y = "row.names")

#Add a population column
geno.and.pheno <- cbind(geno.and.pheno[,1], as.factor(substr(as.character(geno.and.pheno[,1]), start = 3, stop = 4)), geno.and.pheno[,c(2:ncol(geno.and.pheno))])
colnames(geno.and.pheno)[2] <- "pop"

write.table(geno.and.pheno,paste(popByMarker.path,"geno.and.pheno_",i,".txt",sep=''),quote=F,sep='\t')

##############################################################################################################################################
#Fit the JL model and obtain pop*marker terms
#attach(geno.and.pheno)   #cHD commented out 5/14; attaching causes masking and other problems, avoiding. geno.and.pheno was already specified below where needed.

#Get the model order of the model in R to be the same as that in TASSEL.
seq <- (seq(1:nrow(TASSEL.model.results))-2)[-c(1,2)]                 #BY CBK: obtain number of sig SNPs from JL model, but remove 1st pop row) 
model.order <- cbind(seq, TASSEL.model.results[-c(1,2),c(2,3,4)])     #BY CBK: bind num of sig SNPs with model results (chr, position, SNP_name) 
#Sort by chromosome and bp
model.order <- model.order[order(model.order[,4]),]              #BY CBK: order model based on position (seq now different); CHD changed from 3 to 4 on 9/2/19 to match carot 2019 files, for position
model.order <- model.order[order(model.order[,3]),]              #BY CBK: order model by chromosome (and position); CHD changed from 2 to 3 on 9/2/19 to match carot 2019 files, for chromosome
model.order <- cbind(model.order, seq(1:nrow(model.order)))      #BY CBK: add new sequence numbering to existing order model 
#Sort by the first column of marker order so that the markers will be put into the model in the same order as they were selected.
model.order <- model.order[order(model.order[,1]),]              #BY CBK: obtain original numbering
index <- model.order[,5]+3 #3 is added so because the first SNP is on the fourth column             #by CBK: just the model order by chr and position

base.model <- paste(i, " ~ pop",sep = "")                    #By CBK: paste trait and pop+   - note that tilde means "modeled as"; CHD removed "+" after pop 5-14 because + already appended with each marker
#for(k in index){                                                 #BY CBK: cycle through rows in index, which is the sequence ordered by chr and position
 # base.model <- paste(base.model,"+pop:",colnames(geno.and.pheno)[k],sep = "")   #BY CBK:  writing out model of trait modeled by pop plus pop by SNPs in model
#} #end for(k in 4:ncol(geno.and.pheno))   #By MFL: "hard coding" this in because of the issues with half of the files
base.model <- paste(base.model, "+pop:a_carotene+pop:b_carotene+pop:b_cryp+pop:zeaxanthin+pop:phytofluene+pop:lutein")
#JL.model <- lm(AT ~ pop+pop:S_38836669, data = geno.and.pheno) 
JL.model <- lm(paste(base.model, sep = ""), data = geno.and.pheno)    # linear model of trait modeled by pop, pop:SNPs specifying data in geno.and.pheno

pop.by.marker.effects <- summary(JL.model)$coefficients[-c(1:25),]      # obtain descriptive statistics for linear model and model coefficients for pop and SNP within pop, access individual columns with model coefficients
         #CBK changed 26 back to 25 - THIS CONTRIBUTED TO THE LOSS OF POP1:snp1 EFFECT IN ALL POP.BY.MKR.EFFECT FILES
head(pop.by.marker.effects)                                              #header for this data matrix pop term, estimate, stderr, tvalue, prob, FDR adj pvalue

#Perform the FDR procedure on the P-values from pop.by.marker.effects
 res <- mt.rawp2adjp(pop.by.marker.effects[,4], "BH")
 adjp <- res$adjp[order(res$index), ]
                                       
 pop.by.marker.effects <- cbind(pop.by.marker.effects, adjp[,2])
 colnames(pop.by.marker.effects)[5] <- "FDR_Adjusted_P-values"

#Write out the results, put the in the same directory as the JL results from TASSEL
 pop.by.marker.effects <- cbind(rownames(pop.by.marker.effects), pop.by.marker.effects)
 colnames(pop.by.marker.effects)[1] <- "Term"
 
 #setwd(paste(home.dir, "\\validate_CBK.AEL\CHD_Tassel3fromSF_modified0.01\Allelic_Effect_Estimates.no.MultiColl\\", sep = ""))        #commented out CHD 5-11--moved up to top
 write.table(pop.by.marker.effects, paste(popByMarker.path,"Pop.by.Marker.Effect.Estimates.from.R.",i,".SI01.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE) 

#################################################################################################################################################################################
#For loop through all chromosomes, which will fit a separate GLM model for each chromosome and output all residuals

#setwd(paste(home.dir,"\\Residuals_for_GWAS_SI01" ,sep = "")) #CHD commented out 5/14--moved to top.
#test:---model.order <- model.order[-c(1:10), ]

#strategy: for loop through each chromosome. Use the "which" statement to get rid of SNPs on the chromsome being tested. The model.order object will be utilized for this analysis step.
 # Thn fit a model based off of the remaining SNPs in model.order
geno.and.pheno.mod <- geno.and.pheno[complete.cases(geno.and.pheno[,3]),]   #CBK ADDED to obtain geno.and.pheno set without missing trait BLUEs

#resid <- as.data.frame(cbind(seq(1:nrow(geno.and.pheno)),as.character(geno.and.pheno[,1])))         # has only 2 columns - order and geno_code [3432,2]
 
 #using only genos with trait BLUES
resid <- as.data.frame(cbind(seq(1:nrow(geno.and.pheno.mod)),as.character(geno.and.pheno.mod[,1])))

for(j in 1:10){
  model.order.temp <- model.order
  if(length(which(model.order.temp[,3] == j))>0){ #CHD changed from 2 to 3 on 9/2/19 to match chr column in 2019 carot files
   model.order.temp <- model.order.temp[-which(model.order.temp[,3] == j),] #CHD changed from 2 to 3 on 9/2/19 to match chr column in 2019 carot files
  } #end if
  index <- model.order.temp[,5]+3
  base.model <- paste(i, " ~ pop",sep = "")   #CHD removed "+" after pop 5/14 because one is already appended by terms below
 # for(k in index){
#   base.model <- paste(base.model,"+pop:",colnames(geno.and.pheno.mod)[k],sep = "")
 # } #end for(k in 4:ncol(geno.and.pheno))
  base.model <- paste(base.model, "+pop:zeino+pop:a_carotene+pop:b_carotene+pop:phytofluene+pop:lutein+pop:b_cryp+pop:zeaxanthin")
  #Fit the model specific to the chromosome being tested                                                     
  Residual.model <- lm(paste(base.model, sep = ""), data = geno.and.pheno.mod)    #dim(geno.and.pheno... {3432,16}
  #colnames(resid) <- c("V3","V2","V1")             #added by CBK
  
  resid.chr <-as.data.frame(as.matrix(resid(Residual.model)))     #only has V1, residuals for linear model of trait modeled by marker by pop [3314,1]
  #resid.chr <-cbind(geno.and.pheno.mod, as.data.frame(as.matrix(resid(Residual.model) ) ) )

  resid <- merge(resid, resid.chr, by.x = "V1", by.y = "row.names")
  #resid <- cbind(resid, resid.chr)  #ML 11/5: This line is throwing errors! switched with 126
 
}#end for(j in 1:10)

  resid <- resid[,-1]
  colnames(resid) <- c("Sample", "Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10")

  write.table(resid, paste(residuals.path,"Resid_from.R_",i,".SI.01.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
 #write.table(resid, paste(residuals.path,"Resid.",i,".SI.01_redo_ExtrValRemoved_newTrans.txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE) #for dT3 redo

}# End for(i in trait)