rm(list = ls())
#saved in cbk files as: Assess_Pleiotropy_Identify_Genes_in_Common_QTL_Support_Intervals_Non_Sig_Correl_Included_Too_cbk.r


###Required files:
### (1) pop by marker effect estimates
### (2) transformed BLUEs
### (3) imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt (GBS SNPs)
### (4) JL results in tabular summary - current version organized by QTL left bound, preferably with common.SI tracker 


###Output generated:
### (1) "Results_from",trait.1,"JL.Model.txt", which are sublists of trait QTL
### (2) "Results_from",trait.2,"fitted.to.",trait.1,".Markers.txt", which are effect estimates from regression of trait.2 data onto trait.1 markers
### (3) ",trait.1,".QTL","that.are.pleiotropic.with", trait.2,".csv", which contains the correlation of original (trait.1 data onto trait.1 markers) and refit (trait.2 data onto trait.1 markers) effect estimates


#Important comment to read around line 141 regarding significance threshold - depends on type of analysis being run

###############################################################################
pleiotropy.calculator.within.supp.int <- function(trait.1, trait.2, threshold, validation = FALSE){
  print(trait.1)
  print(trait.2)
  print(paste("Threshold is ", threshold, sep=''))
      #Written by Alex Lipka 12/14/2013      
      for(counter in 1:2){
      #Read in the results for Trait 1.

       Results.trait.1 <- read.table(paste(home.dir,pop.X.mkr.dir,"Pop.by.Marker.Effect.Estimates.from.R.",
                                     trait.1, ".SI01.txt", sep = ""), sep='\t',head = TRUE)
       
      if(validation){
         write.table(Results.trait.1, paste(home.dir,pleio.dir,"Results_from",trait.1,"JL.Model.txt",
                      sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
      }
      
      markers.in.trait.1 <- unique(substr(as.character(Results.trait.1[,1]), 7, nchar(as.character(Results.trait.1[,1]))))
      
      #Read in the phenotypic values for Trait 2
      pheno.data = read.table(paste(home.dir,pheno.dir,"BLUEs_No_Outliers_KernelColor.Carot_Transformed_all_forR_",trait.2,".txt", sep = ""), head = TRUE)
      
      colnames(pheno.data) <- c("Geno_Code", trait.2)
      
      ###################################
      #Fit the markers from Trait 1 to Trait 2
      #Parse out all of the SNPs that were selected in the model
      geno.reduced <- genotypes[which(genotypes[,1] %in% markers.in.trait.1),-c(2:5)]
      geno.reduced.formatted <-as.data.frame(t(geno.reduced[,-1]))
      colnames(geno.reduced.formatted) <- as.character(t(geno.reduced[,1]))
      
      
      #pheno.data will always have more data becuase IBM is included in the phenotypic data.
      
      geno.and.pheno <- merge(pheno.data, geno.reduced.formatted, by.x = "Geno_Code", by.y = "row.names")
      
      #Add a population column
      geno.and.pheno <- cbind(geno.and.pheno[,1], as.factor(substr(as.character(geno.and.pheno[,1]), start = 3, stop = 4)), geno.and.pheno[,c(2:ncol(geno.and.pheno))])
      colnames(geno.and.pheno)[2] <- "pop"
      
      
      base.model <- paste(trait.2, " ~ pop",sep = "") #CHD removed '+' here because is already added below
      for(k in markers.in.trait.1){
        base.model <- paste(base.model,"+pop:",k,sep = "")
      } #end for(k in 4:ncol(geno.and.pheno))
      
      #JL.model <- lm(AT ~ pop+pop:S_38836669, data = geno.and.pheno) 
      JL.model <- lm(paste(base.model, sep = ""), data = geno.and.pheno) 
      
      effects.trait.2.fitted.to.trait.1.markers <- summary(JL.model)$coefficients[-c(1:10),]
      
      if(validation){
         write.table(effects.trait.2.fitted.to.trait.1.markers, paste(home.dir,pleio.dir,"/Pleio_Eff.Ests/Results_from",trait.2,"fitted.to.",trait.1,".Markers.txt", sep = ""), sep = "\t", row.names = TRUE, quote = FALSE)
      }
      
      
      head(effects.trait.2.fitted.to.trait.1.markers)
      
      #Perform the FDR procedure on the P-values from effects.trait.2.fitted.to.trait.1
       res <- mt.rawp2adjp(effects.trait.2.fitted.to.trait.1.markers[,4], "BH")
       adjp <- res$adjp[order(res$index), ]
       
       effects.trait.2.fitted.to.trait.1.markers <- cbind(effects.trait.2.fitted.to.trait.1.markers, adjp[,2])
       colnames(effects.trait.2.fitted.to.trait.1.markers)[5] <- "FDR_Adjusted_P-values"
      
      
      #########################
      #Within each marker, calculate the correlation between the effect estimates
      
      #Object names: effects.trait.2.fitted.to.trait.1.markers  Results.trait.1
      
      #Check to see if the pop*marker effects in each row are equal. If not, break out of the loop
      if(length(which(Results.trait.1[,1] != rownames(effects.trait.2.fitted.to.trait.1.markers))) > 0){
          print(paste("When ", trait.2, " was fitted with the effects of ", trait.1, " the pop.marker terms did not match up", sep = ""))
          break
      }                                         
      

      
      #Obtain separate pop and marker terms
      pop.terms <- substr(as.character(Results.trait.1[,1]), 1, 5)
      marker.terms <- substr(as.character(Results.trait.1[,1]), 7, nchar(as.character(Results.trait.1[,1])))
      
      cor.results <- NULL
      #For loop to cycle through the markers identified
      for(i in unique(marker.terms)){
        #In both data sets parse out only effects from the identified marker
        Results.trait.1.subset <- Results.trait.1[which(marker.terms == i), 1:2]
        Results.trait.1.subset[,1] <- as.character(Results.trait.1.subset[,1])
        effects.trait.2.fitted.to.trait.1.markers.subset <-   effects.trait.2.fitted.to.trait.1.markers[which(marker.terms == i), 1]
        
        data.for.correl <- merge(as.data.frame(Results.trait.1.subset), as.data.frame(effects.trait.2.fitted.to.trait.1.markers.subset),
                                 by.x = "Term", by.y = "row.names")
        
        #Generate scalar from lambda values accounting for any sign changes                         #cbk added 7/9/15
        lambda.trait.1 <- lambda.values[which(lambda.values[,2] == trait.1),3]
        lambda.trait.2 <- lambda.values[which(lambda.values[,2] == trait.2),3] 
        scalar <- if(lambda.trait.1 * lambda.trait.2 < 0){-1}else{1}  
        
        #Calculate Pearson correlation coefficient, the P-value from the F statistic and append the results.
        correl_not.sign.corrected <- cor(data.for.correl[,2], data.for.correl[,3])
        print(paste("Correl before sign correction is ", correl_not.sign.corrected,sep=""))
        correl <- correl_not.sign.corrected*scalar
        print(paste("Correl after sign correction is", correl, sep=''))
        dfr <- nrow(Results.trait.1.subset) - 2
        print(dfr)
        r2 <- correl^2
        Fstat <- r2 * dfr / (1 - r2)
        P.val <- 1 - pf(Fstat, 1, dfr)
        cor.results <- rbind(cor.results, c(i, genotypes[which(as.character(genotypes[,1]) == i  ),3],  
                             genotypes[which(as.character(genotypes[,1]) == i  ),4], correl, P.val))
        }#end  for(i in unique(marker.terms))
        
       #write.table(cor.results, paste(home.dir,"//Pleiotropy_Investigation\\Correl.between.pop.marker.from.",trait.1,".in.",
       #            trait.1, ".vs.", trait.2,".txt", sep = ""), sep = "\t", row.names = FALSE, quote = FALSE)
       
       
       
       #Print out the number of QTL with P-values > 0.05
       p.values.cor.results <- as.numeric(as.vector(cor.results[,5]))
       FDR_procedure = mt.rawp2adjp(p.values.cor.results,"BH")
       adj_p_ofCorrels <- FDR_procedure$adjp[order(FDR_procedure$index),2]
       print(adj_p_ofCorrels)
       
       cor.results = cbind(cor.results,adj_p_ofCorrels)
       
       #####EXTREMELY IMPORTANT
       #### Uncomment next line if all results are desired with specified significance threshold, removing all markers that have P-value < the threshold
       cor.results.signif <- as.matrix(cor.results[which(as.numeric(cor.results[,6]) < threshold),])
       #cor.results.signif <- as.matrix(cor.results[which(p.values.cor.results < threshold),])
       #cor.results.signif <- as.matrix(cor.results)
      
       
       #if(length(cor.results.signif) == 5) cor.results.signif <- t(cor.results.signif)            ### Removed use of this line
       
       
       if( (length(cor.results.signif)>0) ){ 
         if(ncol(cor.results.signif)==1){cor.results.signif = t(cor.results.signif)}
        colnames(cor.results.signif) <- c("Peak_SNP", "Chr", "bp", "Correlation_Coefficient", "Raw_P_value_from_F_Statistic","FDR-Adjusted_PValue_ofCorrel")
         
        #Get only the JL results for only trait.1
        JL.Results.trait.1 <- JL.Results[which(JL.Results[,1] == trait.1),]
      
        #Adding additional QTL information
      
        #Gene.list <- NULL
        Left.supp.int <- NULL
        Right.supp.int <- NULL
        PVE_trait1 <- NULL
        for(i in as.factor(cor.results.signif[,1])){ #used to be as.factor
            #Gene.list <- c(Gene.list, as.character(JL.Results.trait.1[which(JL.Results.trait.1[,6] == i),9]))
            Left.supp.int  <- c(Left.supp.int, as.character(JL.Results.trait.1[which(JL.Results.trait.1[,6] == i),4]))
            Right.supp.int  <- c(Right.supp.int, as.character(JL.Results.trait.1[which(JL.Results.trait.1[,6] == i),5]))
            PVE_trait1  <- c(PVE_trait1, as.character(JL.Results.trait.1[which(JL.Results.trait.1[,6] == i),8]))      
        } # end for(i in as.factor(cor.results.signif[,1]))
      
        #cor.results.signif <- cbind(cor.results.signif, Left.supp.int, Right.supp.int, PVE, Gene.list)
        cor.results.signif <- cbind(cor.results.signif, Left.supp.int, Right.supp.int, PVE_trait1)
      
        write.table(cor.results.signif, paste(home.dir,pleio.dir,"Signif_Pleio_Pairwise_Lists/",trait.1,".QTL",
                                              "that.are.SIGNIF.pleiotropic.with",trait.2,"_0.05_sign.corrected_FinJL.txt", sep = ""), row.names = FALSE, sep = "\t",  quote = FALSE)
        }#end if(nrow(cor.results.signif)>0)
      
      #Swap the traits
      old.trait.1 <- trait.1
      old.trait.2 <- trait.2
      
      trait.1 <- old.trait.2
      trait.2 <- old.trait.1
      
      counter <- counter + 1
      
      } #end for(counter in 1:2)
      
      
}###End pleiotropy calculator





###############################################################################


#source("https://bioconductor.org/biocLite.R")
#biocLite("multtest")
library(multtest)
library(stringr)
setwd("C:/Users/chdiep/Documents/NAM_KernelColor_Carot/")
home.dir <- getwd()
tabSummary.path = "/Summary_Tables.Figures/"
pleio.dir = "/Summary_Tables.Figures/Pleiotropy/"
pop.X.mkr.dir = "/popByMarker/"
pheno.dir = "/Master_Geno.Pheno_Inputs/"
geno.dir = "/Master_Geno.Pheno_Inputs/"
lambda.dir = "/Master_Geno.Pheno_Inputs/"

#Read in the gneotypic data
#setwd(paste(home.dir, "\\(9)JL Analysis\\Permutations\\GBS_SNPs\\", sep = ""))
genotypes <- read.table(paste(home.dir,geno.dir,"imputedMarkers.allchr.0.1cm.final.Panzea.consolidated.B.txt",sep=''), head = TRUE)

#Read in the sheet of JL Results
#JL.Results <- read.table("Tab_Sum_Carot_alpha_0.01_20140603_for_Pleio_Script_Comps_Only.txt", head = TRUE)
JL.Results <- read.table(paste(home.dir,tabSummary.path,"Tabular_Summary_of_JL_carot_Results_for_all_traits_SI01_newPVEv2_corrsigns_FinJL.txt",sep=''), sep='\t',head = TRUE)

lambda.values <-read.table(paste(home.dir,lambda.dir,"BLUEs_lambda_values_used_all_forR.txt", sep = ""), head = TRUE)

#trait.multicollinearity <- c("DT3", "AT3", "TOTAL_TOCOCHROMANOLS", "DT3_BY_AT3") #Not needed anymore because we are using alpha = 0.01 threshold
#BLUE.or.BLUP <- "BLUE"  #Options are "BLUE" and "BLUP"
 
#Estimate Pleiotropy
traits <- c('a_carotene',	'b_carotene',	'b_cryp',	'lutein',	'phytofluene',	'zeaxanthin',	'zeino',	'total_carot','traitKC')

#for(p in 1:(length(traits)-1)){
  #for(q in (p+1):length(traits)) pleiotropy.calculator.within.supp.int(trait.1 = traits[p],trait.2 =  traits[q], threshold = 0.01, validation = FALSE)
#} # end  for(p in 1:(length(traits)-1))

for(p in 1:(length(traits)-1)){
  for(q in (p+1):length(traits)) pleiotropy.calculator.within.supp.int(trait.1 = traits[p],trait.2 =  traits[q], threshold = 0.05, validation = FALSE)
} # end  for(p in 1:(length(traits)-1))

