# GO analysis (Fishers)

# load libraries
library(tidyverse)


# CHOOSE SET OF INPUT FILES TO loop through from Venn folder -----------------------------------
data_path <- "~/Documents/GitHub/MPCC_2018/Host_analyses/GO_venn"
results_files <- dir("~/Documents/GitHub/MPCC_2018/Host_analyses/Venn_diagrams/", pattern = "*cold.csv") # get file names
names <- sub("\\.csv.*", "", results_files)
iso2go = read.delim("astrangia_iso2go.tab", sep = "\t")


# Load all result files from folder and mutate to have a 1 for presence, and 0 for absence to prepare for GO analysis
for(i in names){
  filepath = file.path("~/Documents/GitHub/MPCC_2018/Host_analyses/Venn_diagrams/",paste(i,".csv", sep="")) # sets up .csv files
  assign(i, read.delim(filepath, sep = ",") %>% # loads up files
           select(Protein_ID) %>%
           mutate(Fishers = 1) %>%
           full_join(iso2go) %>%
           mutate(Fishers = if_else(is.na(Fishers), 0, 1)) %>%
           select(Protein_ID, Fishers))
}

# Write these new modified files
for( i in 1:length(names)) {
  write.table(get(names[i]), 
              paste0(names[i],
                     "_fishers.csv"),
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
inputFiles = paste0(names, "_fishers.csv")

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
               smallest = 5,   # a GO category should contain at least this many genes to be considered
               clusterCutHeight = 0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
               # Alternative=ALTERNATIVE # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
               # Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
               #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
    )
  }
}



### Plotting results, adjust for go category and analysis

input = "white_only_cold_fishers.csv"
goDivision="BP" 

pdf("./figures/BP_white_only_cold.pdf")
gomwuPlot(input,goAnnotations,goDivision,
                  absValue = 0.001,  # genes with the measure value exceeding this will be counted as "good genes". This setting is for signed log-pvalues. Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  #	absValue=1, # un-remark this if you are using log2-fold changes
                  level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.01, # FDR cutoff to print in regular (not italic) font.
                  level3=0.001, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5,)
dev.off()



# text representation of results, with actual adjusted p-values
results[[1]]

# ------- extracting representative GOs

# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-4 # adjusted pvalue cutoff for representative GO
hcut=0.7 # height at which cut the GO terms tree to get "independent groups". 


# plotting the GO tree with the cut level (un-remark the next two lines to plot)
plot(results[[2]],cex=0.6)
abline(h=hcut,col="red")

ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
  rn=names(ct)[ct==ci]
  obs=grep("obsolete",rn)
  if(length(obs)>0) { rn=rn[-obs] }
  if (length(rn)==0) {next}
  rr=results[[1]][rn,]
  bestrr=rr[which(rr$pval==min(rr$pval)),]
  best=1
  if(nrow(bestrr)>1) {
    nns=sub(" .+","",row.names(bestrr))
    fr=c()
    for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
    best=which(fr==max(fr))
  }
  if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}

mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
bestGOs
