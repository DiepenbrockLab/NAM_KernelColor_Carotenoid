#Set the working directory
setwd("/Users/MaryLaPorte/Documents/Model_Fitting_2019_from.Server")
home.dir <- getwd()
tabSummary.path = paste(home.dir,"/Summary_Tables.Figures/",sep='')
trait.set = "carot"

###Required files:
### (1) Tabular summary file generated in script 3a; IMPORTANT - must be sorted by chromosome and left bound of support interval -- MANUAL STEP
the.data = read.table(paste(tabSummary.path,"Tabular_Summary_of_JL_carot_Results_for_all_traits_SI01_newPVEv2_corrsigns.Sorted.txt", sep =""), sep = "\t",stringsAsFactors=F,head=T)

Common.SI.tracker <- NULL
counter <- 0
Common.bp.left.vector <- NULL
Common.bp.right.vector <- NULL
#For loop through the chromosomes
for(i in unique(the.data[,2])){
  #Obtain all of the rows from the ith chromosome
  the.data.one.chr <- the.data[which(the.data[,2] == i),]
  
  for(j in 1:nrow(the.data.one.chr)){#For loop through the ith chromosome
    if(j == 1){
      #Set the tracker equal to the the next number, use c() to add it to the Common SI tracker
      counter <- counter + 1
      Common.SI.tracker <- c(Common.SI.tracker, counter)
      #Set the bounds of the common support interval
      common.left <- the.data.one.chr[j,4]
      common.right <- the.data.one.chr[j,5]
      Common.bp.left.vector<- c(Common.bp.left.vector, common.left)
      Common.bp.right.vector<- c(Common.bp.right.vector, common.right)
      #move on to the next iteration of this for loop
      next
    }#end if(j == 1)
    #Record the start and end bp of the row. This will be the "common support interval"
    start.bp.this.row <- the.data.one.chr[j,4]
    stop.bp.this.row <- the.data.one.chr[j,5]
    #Below are the four situations in which the JL support interval overlaps
    Situation.1 <- (start.bp.this.row <= common.left)&(stop.bp.this.row >= common.left) #New row SI flanks left common SI boarder
    Situation.2 <- (start.bp.this.row >= common.left)&(stop.bp.this.row <= common.right) #New row SI completely within common SI
    Situation.3 <- (start.bp.this.row <= common.right)&(stop.bp.this.row >= common.right) #New row SI flanks right common SI boarder
    Situation.4 <- (start.bp.this.row <= common.left)&(stop.bp.this.row >= common.right) #common SI completely within new row SI
    if((Situation.1)|(Situation.2)|(Situation.3)|(Situation.4)){
      #Call it the same interval
      Common.SI.tracker <- c(Common.SI.tracker, counter)
      #Obtain the new bounds of the "common support interval"
      common.left <- min(common.left,start.bp.this.row)
      common.right <- max(common.right, stop.bp.this.row)
      Common.bp.left.vector <- c(Common.bp.left.vector, common.left)
      Common.bp.right.vector <- c(Common.bp.right.vector, common.right)      
    }else{
      #Take the previous Common bp.left.vector and Common.bp.right.vectors, and adjust the corresponding support intervals
      Common.bp.left.vector[which(Common.SI.tracker == counter)] = min(Common.bp.left.vector[which(Common.SI.tracker == counter)])
      Common.bp.right.vector[which(Common.SI.tracker == counter)] = max(Common.bp.right.vector[which(Common.SI.tracker == counter)])
      #Start a new common interval
      counter <- counter+1
      Common.SI.tracker <- c(Common.SI.tracker, counter)
      #Record the start and stop bp position as the new "common support interval"
      common.left <- start.bp.this.row
      common.right <- stop.bp.this.row 
      #Add the common positions to Common.bp.left.vector and Common.bp.right.vector
      Common.bp.left.vector<- c(Common.bp.left.vector, common.left)
      Common.bp.right.vector<- c(Common.bp.right.vector, common.right)
    }#end if statement
    
  }#End for(j in 1:nrow(the.data.one.chr))
}#End for(i in unique(the.data[,2]))

#Append the information obtained in this loop to the input data
the.data.final <- cbind(the.data, Common.SI.tracker, Common.bp.left.vector, Common.bp.right.vector)

#Write this information out into a table.
write.table(the.data.final,paste(tabSummary.path,"Tab_Sum_Final_carot_0904_with_Common_SI_Info_left_bound.txt",sep=''), 
            sep = "\t", col.names = TRUE,row.names = FALSE, quote = FALSE)
