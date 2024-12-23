#### This script takes the parameter estimates from the best runs of the best models for FSC runs on
####    Larix datasets and converts the units and makes a table of them
####    this *should* work with any number of populations with divergence times, secondary contact times
####    and migration rates: requires that migration rates are specified such that the parameter names include the numbers of 
####    the populations between which migration is occurring and that these same numbers are included in the 
####    population names in a consistent way: e.g., MIG10 is migration bwtween POP1 and POP0; MIG102 is migration between 
####    an ancestral POP10 (or POP10anc - just so long as the numbers an only those numbers are there) amd POP2
####    (in an alternate tree structure, that could be migration between POP1 and POP02 instead - I beleive I have this correctly handled to parse out the correct pops)


#### This should be effective for any number of populations (of course, assuming they're structured exactly the way I 
####      name things in FSC), including migration among ancestral populations



library(dplyr)
library(stringr)

#### Specify the generation time
gen_time <- 180


# ##### set up directories/files here
# bestlhoods_dir <- "~/Active_Research/Larix_data_and_outs/ipyrad_31k_out/"
# boot_in_dir <- "~/Active_Research/Larix_data_and_outs/ipyrad_31k_boot_inputs" # directory where we will create input files to run bootstraps
# sfs_file <- "~/Active_Research/Larix_data_and_outs/SFS/fastsimcoal2/LarK3_MSFS.obs" # directory of SFS input for original modeling
# 
# # # Files from the bestrun of the best model
# # bestrun_dir <- paste0(bestlhoods_dir, "Rate_noMK3")
# 
# # name the output file for converted parameter estimates
# out_ests <- "Larix_ipyrad_31k_pars_conv.csv"



# ##### set up directories/files here
# bestlhoods_dir <- "~/Active_Research/Larix_data_and_outs/stacks_31k_out/"
# boot_in_dir <- "~/Active_Research/Larix_data_and_outs/stacks_31k_boot_inputs/" # directory where we will create input files to run bootstraps
# sfs_file <- "~/Active_Research/Larix_data_and_outs/SFS_stacks/fastsimcoal2/LarStkK3_MSFS.obs" # directory of SFS input for original modeling
# 
# # Files from the bestrun of the best model
# bestrun_dir <- paste0(bestlhoods_dir, "31k_AncMCoInExpK3")
# 
# # name the output file for converted parameter estimates
# out_ests <- "Larix_Stacks_31k_pars_conv.csv"


##### set up directories/files here
bestlhoods_dir <- "~/Active_Research/Larix_data_and_outs/ipyrad_31k_no_mono/"
# boot_in_dir <- "~/Active_Research/Larix_data_and_outs/stacks_31k_boot_inputs/" # directory where we will create input files to run bootstraps
sfs_file <- "~/Active_Research/Larix_data_and_outs/X_OTHER_DATASETS_calibrations/SFS/fastsimcoal2/LarK3_MSFS.obs" # directory of SFS input for original modeling

# Files from the bestrun of the best model
# bestrun_dir <- paste0(bestlhoods_dir, "31k_AncMCoInExpK3")

# name the output file for converted parameter estimates
out_ests <- "Larix_ipyrad_31k_no_mono_pars_conv.csv"



maxLfile <- list.files(path = bestrun_dir, pattern = "maxL.par$", full.names = TRUE)
pv_to_copy <- list.files(path = bestrun_dir, pattern = ".pv$", full.names = TRUE)
est_to_copy <- list.files(path = bestrun_dir, pattern = ".est$", full.names = TRUE)
tpl_to_copy <- list.files(path = bestrun_dir, pattern = ".tpl$", full.names = TRUE)





###########################################################################################################################################################
### Convert parameter estimates - this works with relatively simple models that don't have migration between ancestral pops, that will complicate things
############################################################################################################################################################

# set working directory
setwd(bestlhoods_dir)




## List all of the .bestlhoods files
files_lnl <- list.files(recursive=TRUE, pattern=".bestlhoods")
## Read in these files
output <- lapply(files_lnl, read.table, header=TRUE)
names(output) <- gsub("\\.bestlhoods", "", files_lnl) # name each element

## read in the csv of the model fits to look at it
mod_fits <- read.csv("FSC_model_Fits.csv" )
mod_fits

## Let's just get estimates for all of the models - this code will not universally convert, but will handle divergence times, population sizes, and present migration rates
##     have now added in functionality to do migration rates among ancestral populations - this works only if the migration rates are structure such that ancestral migration rates 
##     between population 1 and ancestor of pops 2&3 is written as MIG123$ or MIG132$ and the ancestral pop size of 
mods_convert<-names(output)
converted_ests<-vector("list", length(mods_convert)) # empty list to dump coverted output into 
for(j in 1:length(mods_convert)){ # loop over the models doing the conversions we want
  raw_est<-output[mods_convert[j]][[1]] # where the raw estimates are
  converted<-rep(NA, length(raw_est)-2) # Make an empty object that the converted numbers will end up in 
  names(converted)<-gsub("\\.", "", names(raw_est)[1:(length(raw_est)-2)]) # add on the parameter names
  names(raw_est)<-gsub("\\.", "", names(raw_est)) ## Strip out the . in the original names, too
  
  ### Go through and populate the new matrix with values converted from the raw values
  ind_npop<-grep("NPOP", names(converted)) # indices of elements that are population sizes
  converted[ind_npop]<-round(raw_est[ind_npop]/2, -2) # pop sizes get divided by 2 to be diploid and round to hundreds place
  ind_times<-grep("TDIV|TCONT|TSEP|TRESIZE|TMIG", names(converted))  # index of elements that are times (divergences or secondary contact)
  converted[ind_times]<-round(raw_est[ind_times]*gen_time, -2) ## times get multiplied by generation time - round to hundreds place also
  # migration rates are multipled by 2 and by population size of the population that the migration is into (e.g., MIG20 is migration 2->0 and is multiplied by pop size 0 NPOP0)
  # this is slightly more complex, so handle with a loop
  ind_mig<-grep("MIG", names(converted)) # indices of elements that are migration rates
  ind_migA<-grep("MIGA", names(converted))# after handling standard migration rates, handle migration rates among the same pops at an earlier time - I designate these MIGA
  if(length(ind_migA)>0){ # if we have MIGA in here, then these need to be removed from ind_mig and handled separately
    ind_mig <- ind_mig[!ind_mig %in% ind_migA]
  }
  # for handling any migration among ancestral populations, need to get out the number of digits in each population name to identify ancestral populations - do this outside of the loop below to avoid repeating calculations
  pop_names<-names(converted)[ind_npop] # get population names
  nums_pops<-str_split(str_extract(pop_names, "[[:digit:]]+"), pattern="")
  num_digits_pops<-sapply(nums_pops, function(x) sum(!is.na(x))) # get the number of numerical digits that are not NA
  ind_pop_names_anc<-which(num_digits_pops>1) # indices of pop names that are ancestral pops
  dig_anc<-nums_pops[ind_pop_names_anc]  # digits in ancestral pops
  names(dig_anc)<-pop_names[ind_pop_names_anc] # name those
  for(i in ind_mig){
    if(exists("pop1_anc")){ # if we've run this function before or a previous iteration of the loop defined pop1_anc and pop2_anc, clear them out
      rm(pop1_anc)
    }
    if(exists("pop2_anc")){
      rm(pop2_anc) 
    }
    ## Need to get and count the number of digits in the migration rate name, if it's more than 2, then ancestral populations are involved
    ##     this works for 4 popuations where the tree is symmetric right now... gets even more complicated past that
    nums_mig<-str_split(str_extract(names(converted)[i], "[[:digit:]]+"), pattern="")[[1]]
    if(length(nums_mig)==2){
      pop_mult<-gsub("MIG.", "NPOP", names(converted)[i]) # for each index, replace MIG. (i.e.,. MIG and the following character) with NPOP resulting in e.g., changing MIG20 to NPOP0, the population we want to multiply by
      converted[i]<-round(raw_est[i]*2*raw_est[pop_mult], 2) # multiply the migration rate by 2 and the pop size of interest - round to 2 decimal places
    }
    if(length(nums_mig)>2){ 
      ## we don't know where the split is among the two ancestral populations, so we'll try all possible groupings
      for(splt in 2:length(nums_mig)){
        ptpop1<-nums_mig[1:splt-1] #potential ancestral population 1
        ptpop2<-nums_mig[splt:length(nums_mig)] # potetential ancestral population 2
        ## Now we need to test if each of these populations exists
        if(length(ptpop1)>1){ # if there is more than 1 number, check if all of them are in any of the population names
          for(m in 1:length(dig_anc)){
            anc<-unlist(dig_anc[m])
            # if(length(dig_anc)==1){  ## Think this was a failed attempt at handling what I handle with unlist in the previous line now
            #   anc<-anc[[1]]
            # }
            if(sum(ptpop1 %in% anc)==length(ptpop1) && length(ptpop1)==length(anc)){ # check that all numbers in the potential population based on splitting out migration numbers are in the numbers of an ancestral population, and that the lengths of these are the same (e.g., without checking the lengths, we could find 01 in an ancestral populaion 012, but that's not what we're looking for)
              pop1_anc<-names(dig_anc)[m] # if we've found the right ancestral population, assign it to pop1_anc object
            }
          }
        }else{ # if there is not more than 1 number (i.e., there is only 1) then that number is a single present population and necessarily exists
          pop1_anc<-paste0("NPOP", ptpop1)
        }
        if(length(ptpop2)>1){
          for(m in 1:length(dig_anc)){
            anc<-unlist(dig_anc[m])
            # if(length(dig_anc)==1){  ## Think this was a failed attempt at handling what I handle with unlist in the previous line now
            #   anc<-anc[[1]]
            # }
            if(sum(ptpop2 %in% anc)==length(ptpop2) && length(ptpop2)==length(anc)){ # check that all numbers in the potential population based on splitting out migration numbers are in the numbers of an ancestral population, and that the lengths of these are the same (e.g., without checking the lengths, we could find 01 in an ancestral populaion 012, but that's not what we're looking for)
              pop2_anc<-names(dig_anc)[m] # if we've found the right ancestral population, assign it to pop1_anc object
            }
          }
        }else{
          pop2_anc<-paste0("NPOP", ptpop2)
        }
        if(exists("pop1_anc") && exists("pop2_anc")){ # for a given split, if both pops exist, then break and we're done
          break
        }else{ # if not, remove these objects (if either exists) - in any case where we split out the first or last digit, that population definitely does exist, but we only want to keep COMBINATIONS where both populations exist
          if(exists("pop1_anc")){
            rm(pop1_anc)
          }
          if(exists("pop2_anc")){
            rm(pop2_anc) 
          }
        }
      }
      converted[i]<-round(raw_est[i]*2*raw_est[pop2_anc], 2) # multiply the migration rate by 2 and the pop size of interest - round to 2 decimal places
    }
  }
  if(length(ind_migA)>0){
    for(i in ind_migA){
      if(exists("pop1_anc")){ # if we've run this function before or a previous iteration of the loop defined pop1_anc and pop2_anc, clear them out
        rm(pop1_anc)
      }
      if(exists("pop2_anc")){
        rm(pop2_anc) 
      }
      ## Need to get and count the number of digits in the migration rate name, if it's more than 2, then ancestral populations are involved
      ##     this works for 4 popuations where the tree is symmetric right now... gets even more complicated past that
      nums_mig<-str_split(str_extract(names(converted)[i], "[[:digit:]]+"), pattern="")[[1]]
      if(length(nums_mig)==2){
        pop_mult<-gsub("MIGA.", "NPOP", names(converted)[i]) # for each index, replace MIGA. (i.e.,. MIGA and the following character) with NPOP resulting in e.g., changing MIGA20 to NPOP0, the population we want to multiply by
        converted[i]<-round(raw_est[i]*2*raw_est[pop_mult], 2) # multiply the migration rate by 2 and the pop size of interest - round to 2 decimal places
      }
      if(length(nums_mig)>2){ 
        ## we don't know where the split is among the two ancestral populations, so we'll try all possible groupings
        for(splt in 2:length(nums_mig)){
          ptpop1<-nums_mig[1:splt-1] #potential ancestral population 1
          ptpop2<-nums_mig[splt:length(nums_mig)] # potetential ancestral population 2
          ## Now we need to test if each of these populations exists
          if(length(ptpop1)>1){ # if there is more than 1 number, check if all of them are in any of the population names
            for(m in 1:length(dig_anc)){
              anc<-unlist(dig_anc[m])
              if(length(dig_anc)==1){
                anc<-anc[[1]]
              }
              if(sum(ptpop1 %in% anc)==length(ptpop1) && length(ptpop1)==length(anc)){ # check that all numbers in the potential population based on splitting out migration numbers are in the numbers of an ancestral population, and that the lengths of these are the same (e.g., without checking the lengths, we could find 01 in an ancestral populaion 012, but that's not what we're looking for)
                pop1_anc<-names(dig_anc)[m] # if we've found the right ancestral population, assign it to pop1_anc object
              }
            }
          }else{ # if there is not more than 1 number (i.e., there is only 1) then that number is a single present population and necessarily exists
            pop1_anc<-paste0("NPOP", ptpop1)
          }
          if(length(ptpop2)>1){
            for(m in 1:length(dig_anc)){
              anc<-unlist(dig_anc[m])
              if(length(dig_anc)==1){
                anc<-anc[[1]]
              }
              if(sum(ptpop2 %in% anc)==length(ptpop2) && length(ptpop2)==length(anc)){ # check that all numbers in the potential population based on splitting out migration numbers are in the numbers of an ancestral population, and that the lengths of these are the same (e.g., without checking the lengths, we could find 01 in an ancestral populaion 012, but that's not what we're looking for)
                pop2_anc<-names(dig_anc)[m] # if we've found the right ancestral population, assign it to pop1_anc object
              }
            }
          }else{
            pop2_anc<-paste0("NPOP", ptpop2)
          }
          if(exists("pop1_anc") && exists("pop2_anc")){ # for a given split, if both pops exist, then break and we're done
            break
          }else{ # if not, remove these objects (if either exists) - in any case where we split out the first or last digit, that population definitely does exist, but we only want to keep COMBINATIONS where both populations exist
            if(exists("pop1_anc")){
              rm(pop1_anc)
            }
            if(exists("pop2_anc")){
              rm(pop2_anc) 
            }
          }
        }
        converted[i]<-round(raw_est[i]*2*raw_est[pop2_anc], 2) # multiply the migration rate by 2 and the pop size of interest - round to 2 decimal places
      }
    }
  }
  if("MUTRATE" %in% names(converted)){ # if MUTRATE is a parameter (i.e., if we let FSC estimate the mutation rate rather than fixing it)
    converted["MUTRATE"] <- raw_est["MUTRATE"]   ## then just add in the mutation rate - note that this is still mutation rate per generation
  }  
  converted_ests[[j]]<-unlist(converted)
  names(converted_ests)[[j]]<-mods_convert[j]
}

converted_table<-as.data.frame(do.call(bind_rows, converted_ests)) # combine into 1 big table
rownames(converted_table)<-names(converted_ests) # add rownames onto it
AIC_vals<-read.csv("FSC_model_Fits.csv") # read in the AIC values, use this to order the table ot estimates
converted_table<-converted_table[gsub("\\.bestlhoods", "", AIC_vals[,1]),]
vals_and_AIC<-cbind(AIC_vals[,c("AIC", "deltaAIC", "AICw")], converted_table) # combine the AIC values into it
write.csv(vals_and_AIC, paste0(out_ests))




# ########################################################################################################################
# ######   Set up bootstrapping                                                                                   ########
# ########################################################################################################################


# # read in the SFS file
SFS_all<-readLines(sfs_file)
# get just the line with the actual SFS numbers
SFS<-SFS_all[3]
# split each line into numbers
SFS_nums<-strsplit(SFS, split=" ")
# sum these numbers to find the total number of SNPs - note this includes monomorphics
num_SNPs<-round(sum(as.numeric(SFS_nums[[1]])), 0)


 
# read in each maxL_file change number of independent loci to the number of SNPs to make it ready for bootstrapping, write to where the bootstraps are
# Commented out block below loops this for multiple models, here, I have strong support for a single model, so only going to boot that one model
par_to_edit <- readLines(maxLfile) # read in the file
line_to_edit <- grep("Number of independent loci",par_to_edit)+1 # get the line number for the header, then add 1 to get the line we actually want to edit
par_to_edit[line_to_edit] <- gsub("1 0", paste(num_SNPs ,0), par_to_edit[line_to_edit]) # edit in the number of SNPs
freqline_to_edit <- grep("per Block:data type",par_to_edit)+1 # edit freq to DNA so that SFS will be simulated - start with finding line to edit again
par_to_edit[freqline_to_edit] <- gsub("FREQ", "DNA", par_to_edit[freqline_to_edit])
setwd(boot_in_dir)
writeLines(par_to_edit, basename(maxLfile)) # write this out to the new file

################################################################################################
## Then go use FSC to generate the bootstrap reps. Run the following commands in terminal:
################################################################################################
#$ cd /Users/harrington/Active_Research/Larix_data_and_outs/stacks_31k_boot_inputs
#$ fsc26 -i 31k_AncMCoInExpK3_maxL.par -n100 -j -m -s0 -x –I -q -u
################################################################################################




### This is the old loop - untested recently:
# setwd(mods)  # set the working directory to where the models are
# # loop over the maxL_files to make the necessary edit
# for(i in 1:length(maxL_files)){
#   par_to_edit<-readLines(maxL_files[[i]]) # read in the file
#   line_to_edit<-grep("Number of independent loci",par_to_edit)+1 # get the line number for the header, then add 1 to get the line we actually want to edit
#   par_to_edit[line_to_edit]<-gsub("1 0", paste(num_SNPs ,0), par_to_edit[line_to_edit]) # edit in the number of SNPs
#   freqline_to_edit<-grep("per Block:data type",par_to_edit)+1 # edit freq to DNA so that SFS will be simulated - start with finding line to edit again
#   par_to_edit[freqline_to_edit]<-gsub("FREQ", "DNA", par_to_edit[freqline_to_edit])
#   writeLines(par_to_edit, paste0(boot_in_dir, "/", basename(maxL_files[[i]]))) # write this out to the new file
# }
# 



# Copy necessary files into the directory made by the FSC command above
reps_dir <- list.files(path = boot_in_dir, pattern = "_maxL$", include.dirs = TRUE, full.names = TRUE)

files_copy <- c(pv_to_copy, est_to_copy, tpl_to_copy)

for(i in files_copy){ # copy each file
  new_name <- gsub("\\.", "_maxL.", basename(i)) # add "_maxL" into each file to match the SFS files generated by FSC
  file.copy(i, paste0(reps_dir, "/", new_name)) # make the copy with the new name
}




  
