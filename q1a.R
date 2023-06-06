# Load libraries
library(robumeta)
library(metafor)
library(dplyr)
library(effsize)
library(MASS)
library(gplots)
library(RColorBrewer)
library(sfsmisc)
library(ggplot2)
library(ggrepel)
library(sfsmisc)
library(compositions)
library(pcaPP)

#setwd("C:\\Users\\anany\\Desktop\\semester8\\hmds\\Ass2")
# Load data
load("Assignment2.RData")

# Filter data based on selected studies and age
df_filtered <- df_ForMetaAnalysis %>%
  filter(study_name %in% SelectedStudies & age >= 60)



compute_meta_corr_group <- function(data,feature_list,metadata_var,grouping_var,grouping_list)
{
  return_out <- as.data.frame(matrix(NA,length(feature_list),10))
  rownames(return_out) <- feature_list
  colnames(return_out) <- c("beta","pval","ci.ub","ci.lb","tau2","QE","QEp","qval","dir","consistency")
  return_out[,1] <- 0
  return_out[,2] <- 1
  return_out[,3] <- 0
  return_out[,4] <- 0
  return_out[,5] <- 0
  return_out[,6] <- 0
  return_out[,7] <- 1
  return_out[,10] <- 0
  
  for(i in 1:length(feature_list))
  {
    species_name <- feature_list[i]
    #print(species_name)
    tryCatch(               
      expr = {                     
        temp_res <- compute_meta_corr(data,species_name,metadata_var,grouping_var,grouping_list)
        print(species_name)
        return_out[i,"beta"] <- temp_res$model$beta
        return_out[i,"pval"] <- temp_res$model$pval
        return_out[i,"ci.ub"] <- temp_res$model$ci.ub
        return_out[i,"ci.lb"] <- temp_res$model$ci.lb
        return_out[i,"tau2"] <- temp_res$model$tau2
        return_out[i,"QE"] <- temp_res$model$QE
        return_out[i,"QEp"] <- temp_res$model$QEp
        return_out[i,"consistency"] <- length(which(sign(temp_res$df_studies[temp_res$df_studies$ri!=0,"ri"])==sign(as.numeric(temp_res$model$beta))))/length(temp_res$df_studies[temp_res$df_studies$ri!=0,"ri"])
        
      },
      error = function(e){    
        print(e)
        print("Error observed. Moving to next")
      },
      finally = {            
        print("finally Executed")
      }
    )
  }
  return_out$qval <- p.adjust(return_out$pval,method="fdr")
  return_out$dir <- ifelse(return_out$qval <= 0.1,3*sign(return_out$beta),ifelse(return_out$pval <= 0.08,2*sign(return_out$beta),sign(return_out$beta)))
  return_list <- list("model" = temp_res$model,"df_studies"=temp_res$df_studies)
  #return(return_list)
  return(return_out)
}

compute_meta_corr <- function(data,var1,var2,grouping_variable,grouping_list)
{
  temp_meta <- data.frame(matrix(NA,length(grouping_list),3))
  colnames(temp_meta) <- c("dataset","ri","ni")
  for(i in 1:length(grouping_list))
  {
    group <- grouping_list[i]
    temp_meta[i,1] <- group
    dat1 <- data[data[,grouping_variable]==group,var1]
    dat2 <- data[data[,grouping_variable]==group,var2]
    temp_meta[i,2] <- cor.fk(dat1,dat2)
    temp_meta[i,3] <- length(dat1)
  }
  temp_meta <- mutate(temp_meta,study_id=grouping_list)
  rownames(temp_meta) <- grouping_list
  print(temp_meta)
  temp_meta <- escalc(measure="ZCOR",ri=ri,ni=ni,data=temp_meta)
  temp_meta <- temp_meta[!is.na(temp_meta$ri),]
  res <- rma(yi, vi, data=temp_meta)
  res$ids <- rownames(temp_meta)
  res$slabs <- rownames(temp_meta)
  return_list <- list("df_studies"=temp_meta,"model"=res)
  return(return_list)
}

compute_detection <- function(data,var1_list,grouping_variable,grouping_list)
{
  detection_matrix <- data.frame(matrix(0,length(var1_list),length(grouping_list)))
  rownames(detection_matrix) <- var1_list
  colnames(detection_matrix) <- grouping_list
  for(i in 1:length(var1_list))
  {
    var1 <- var1_list[i]
    for(j in 1:length(grouping_list))
    {
      group <- grouping_list[j]
      detection_matrix[i,j] <- length(which(data[data[,grouping_variable]==group,var1]>0))/length(data[data[,grouping_variable]==group,var1])
    }
  }
  return(detection_matrix)
}


# Compute meta-analysis results for selected species
# data,feature_list,metadata_var,grouping_var,grouping_list
significant_results <- compute_meta_corr_group(df_filtered, SelectSpecies, "age", "study_name",SelectedStudies )


# Filter for significant results
#significant_results <- significant_results %>%filter(pval <=  0.05)

significant_results <- significant_results %>%
  filter(pval <= 0.05) %>%
  mutate(pval_corrected = p.adjust(pval, method = "fdr"))
 
# Print significant results
print(significant_results)

