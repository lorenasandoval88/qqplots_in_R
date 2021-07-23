
# Added api calls to retrieve BigQuery tables
#install.packages("tidyverse")
# install.packages("bigrquery")
# install.packages('devtools')
# devtools::install_github("r-dbi/bigrquery")
library(bigrquery)
library(boxr)
library(tidyverse)
library(lubridate)
library(jsonlite)
setwd("C:/Users/sandovall2/Downloads")
# data from ER- GWAS ------------------------------------------------------------

# bigquery authorization - please select your NIH email in the console below
bq_auth()
# project
project <- "nih-nci-dceg-connect-dev"
# query
sql <- "SELECT * FROM `nih-nci-dceg-connect-dev.ER_data.ER_neg_qqplot`"
tb <- bq_project_query(project, sql)
BigQuery_data = bq_table_download(tb, bigint = c("character"))
data=BigQuery_data[order(BigQuery_data$P_value),]


# create and save expected p-values in a variable
exp_P_value= ppoints(length(data$P_value))

# save the log of expected p-values in dataframe
data$exp_P_value_log = -log10(exp_P_value)

# pull ~10,000 smallest expected p-values 
low_e = head(data, 9998)

# weighted expected p values --Brian's instructions

# get index of a weighted random sample from the full exp p values
random_index = floor((as.numeric(data$ID) * length(data$P_value) ^ 0.5 / 10000) ^ 2) # get ids for weighted exp p values

# subset index between 9999-length(exp p values vector)
random_index_subset = random_index[random_index > 9998 & random_index <= length(data$exp_P_value_log)]# filter out pval
high_e = data[random_index_subset,]

# combine lowest 9,998 p values with random sample of 9,661 higher p values
new_data = rbind(low_e,high_e)

# add log of pvalues column and remove p_value column
new_data$observed_P_value_log = -log10(new_data$P_value)
library(dplyr)
new_data = dplyr::select(new_data, -c(P_value))

# # distribution of -log10 expected pvalues
# hist(new_data$exp_P_value_log)
# 
# # qqplot function
# qq <- function(pvector, ...)
# {
#   # observed p values
#   o <- -log10(new_data$P_value)
# 
#   # expected p values
#   e= new_data$exp_P_value_log
#  
#   plot(e,o,pch=19,cex=1, xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))), xlim=c(0,max(e)), ylim=c(0,max(o)), ...)
#   abline(0,1,col="gray")
#   mtext(paste("lambdaGC =",round(qchisq(1-median(pvector),1)/qchisq(0.5,1),4)),line=-1)
#   mtext(paste("N =", length(pvector)))
# }
# 
# # run qq plot function
# #observed pval 
# o = data$P_value
# qq(o)

# format data to JSON like file----------------

######### make the column names a row and then remove them
newdata = new_data
names <- colnames(newdata)
# data_a[2:nrow(data_a)+1,] <- data_a
# data_a[1,] <- names
colnames(newdata) <- NULL

library('rlist')
# make list of lists for data:
lis = list()
for(i in seq(1,nrow(newdata))){
  #print(data_a[i,])
  lis = list.append(lis, newdata[i,])
  
  #print(lis)
  #break
}

# add key "data" to list of lists
# add key "columns" and list of column names
library(jsonlite)
lst <- list(columns = names,data = lis)
lst = toJSON(lst, pretty=T,digits=400)


# save data as json file
#  #############  oauth for box
box_auth(client_id = "xxx", client_secret = "xxx")

write(lst, "ERneg_qqplot_data.json")
write.csv(lst, "ERneg_qqplot_data.csv")
box_write(lst, "ERneg_qqplot_data.json", dir_id=141602388994)
