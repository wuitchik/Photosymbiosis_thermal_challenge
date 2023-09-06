# Much of this script was stolen from Dixon (https://github.com/grovesdixon/Acropora_gene_expression_meta)


# load libraries
library(tidyverse)

# CHOOSE SET OF INPUT FILES TO LOOP RUN -----------------------------------
data_path <- "~/Documents/GitHub/MPCC_2018/Host_analyses/GO_DEGs/Astrangia"
results_files <- dir("~/Documents/GitHub/MPCC_2018/Host_analyses/DESeq2/Astrangia/", pattern = "*results.csv") # get file names
names <- sub("\\.csv.*", "", results_files)

# Load all result files from folder
for(i in names){
  filepath = file.path("~/Documents/GitHub/MPCC_2018/Host_analyses/DESeq2/Astrangia/",paste(i,".csv", sep="")) # sets up .csv files
  assign(i, read.delim(filepath, sep = ",") %>% # loads up files
    mutate(mutated_p = -log(pvalue), # modifies to a signed p-value
           mutated_p_updown = ifelse(log2FoldChange < 0,
                                     mutated_p*-1,
                                     mutated_p*1)) %>%
    select(X, mutated_p_updown) %>%
    rename(GeneID = X) %>%
    na.omit() 
    )
}

# Write these new modified files
for( i in 1:length(names)) {
  write.table(get(names[i]), 
             paste0(names[i],
             "_modified_pvalues.csv"),
             sep = ",",
             row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE) # this is critical for GO MWU to read it properly
}

# Now we can loop our GO analyses

# SET BASIC VARS ---------------------------------------------------------
goDatabase = "go.obo"
goAnnotations = "astrangia_iso2go.tab"
divisions = c('CC',
              'MF',
              'BP')

source("gomwu.functions.R")

####################  OVERALL SETS #################### 
inputFiles = paste0(names, "_modified_pvalues.csv")

# LOOP THROUGH SELECTED INPUT FILES AND GO DIVISIONS ---------------------------------------
for (goDivision in divisions){
  print('--------------')
  print('--------------')
  print("------WORKING ON A NEW GO CATEGORY-----")
  print(goDivision)
  print('--------------')
  print('--------------')
  for (input in inputFiles){
    print('--------------')

    print("input being worked on")
    print(input)
    
    ## Modify your GO MWU parameters here
    gomwuStats(input, goDatabase, goAnnotations, goDivision,
               perlPath = "perl", # replace with full path to perl executable if it is not in your system's PATH already
               largest = 0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
               smallest = 10,   # a GO category should contain at least this many genes to be considered
               clusterCutHeight = 0.01, # threshold for merging similar (gene-sharing) terms. See README for details.
               # Alternative=ALTERNATIVE # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
               # Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
               #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
    )
  }
}



## This script kinda throws files all over the place, and so I just lazily moved everything and renamed files using the terminal because I'm like that.

