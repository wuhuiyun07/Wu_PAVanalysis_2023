##Updated9/14/23
##dPCR Data back calculation Pipeline
library(tidyverse)
library(dplyr)

#####set work directory to combine the plates####
setwd("/Users/huiyunwu/Desktop/Virus_particle/data/raw/")

# df <-
#   list.files( pattern = "*.csv") %>% 
#   map_df(~read_csv(.))
# df
# 
# write.csv(df,file="/Users/huiyunwu/Desktop/Virus_particle/data/processed_data/S3toS32_raw_combined.csv",row.names = F)

#####set work directory to the original one
setwd("/Users/huiyunwu/Desktop/Virus_particle/")

read_csv(file="data/processed_data/S3toS32_raw_combined.csv")

## File path in, current location of the .csv file
FilePathIn <-"data/processed_data/S3toS32_raw_combined.csv"
x<-FilePathIn
##File Path out, include desired file name with .csv designation
FilePathOut <-"/Users/huiyunwu/Desktop/Virus_particle/data/processed_data/S3to32_back.cal_dup.csv"
y<-FilePathOut
######Load sample volume#####
Vol.<-read_csv(file = "data/processed_data/Vol.Sample.csv",
               col_types = cols(collection.date =col_date(format = "%m/%d/%y"),
                 name = col_character()))
#####################################


##Pull Data into environment
rawdata=read_csv(file=x) %>% 
  rename_all(tolower) ##to lower case
colnames(rawdata)
# [1] "run"         "date"        "instrument"  "plate"      
# [5] "group"       "name"        "well"        "total"      
# [9] "dye"         "target"      "conc. cp/ul" "sd"         
# [13] "cv%"         "95%ci"       "positives" 

#####generate full dataset for back calculation###
back.cal<-left_join(rawdata,Vol.)
dPCRPrepVol.ul<- 10
dPCRNAVol.ul<-5
NAExntVol.ul<- 100
back.cal$finalcon.cp.L<-(back.cal$`conc. cp/ul`*dPCRPrepVol.ul)/dPCRNAVol.ul*NAExntVol.ul/back.cal$vol.L
back.cal$finalcon.cp.L<-format(back.cal$finalcon.cp.L,scientific = TRUE)
# # colnames(back.cal)
# [1] "run"                           "date"                         
# [3] "instrument"                    "plate"                        
# [5] "group"                         "name"                         
# [7] "well"                          "total"                        
# [9] "dye"                           "target"                       
# [11] "conc. cp/ul"                   "sd"                           
# [13] "cv%"                           "95%ci"                        
# [15] "positives"                     "wwtp"                         
# [17] "collection.date"               "pore.size"                    
# [19] "extn.type"                     "vol.L"                        
# [21] "HFUF eluate.mL"                "HFUF eluate for extraction.mL"
# [23] "Extraction Date"               "comments"                     
# [25] "finalcon.cp.L" 

##Create and print output table
#Final<-data.frame(rawdata$SpeciemenNameColumn,rawdata$NameofOutputColumn,Results)
# output=back.cal[,c("run","name","wwtp","collection.date","pore.size","vol.L","conc. cp/uL","positives","finalcon.cp.L")]
output<-back.cal %>% 
  select(run, name, wwtp, collection.date, pore.size, vol.L, `conc. cp/ul`,positives, total, finalcon.cp.L)

##Summary statistics##

#Check for those will less than 3 wells (in positives column) here

qualitycheck<-ifelse(output$positives >=3 & output$total >= 20000,"Pass","Fail")

outputqc=cbind(output,qualitycheck)

##Save data to Specified output
write.csv(outputqc,file= y,row.names=FALSE)

