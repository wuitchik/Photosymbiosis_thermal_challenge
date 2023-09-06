## GO MWU analyses
library(tidyverse)

# set up 
goDatabase = "go.obo"
goAnnotations = "astrangia_iso2go.tab"
source("gomwu.functions.R")

#set variables to run for each module
ll = load('moduleInputFiles.Rdata')
ll
inputFiles
divisions=c("BP", "CC", "MF")
#c('input1.csv', 'input2.csv')


# run go_mwu for each -----------------------------------------------------
for (goDivision in divisions){
  print('==============')
  print(goDivision)
  #set BP smallest to 50 so you don't get too many
  if (goDivision=='BP'){
    SMALLEST=50
  } else {
    SMALLEST=10
  }
  for (input in inputFiles){
    print('--------------')
    print(paste(input, '...', sep='.'))
    # Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
    gomwuStats(input, goDatabase, goAnnotations, goDivision,
               perlPath="perl", # replace with full path to perl executable if it is not in your system's PATH already
               largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
               smallest=SMALLEST,   # a GO category should contain at least this many genes to be considered
               clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
               Alternative="g" # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
               # Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
               #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
    )
  }
}


##### Left off here
save.image(file = "feb_3.RData")

xxxnxNXxxx# record how many significant and save
divRec = c()
inRec = c()
sigRec = c()
for (goDivision in divisions){
  for (input in inputFiles){
    resName = paste(paste('MWU', goDivision, sep = "_"), input, sep = "_")
    go.res=read.table(resName, header = T)
    totSig = sum(go.res$p.adj < 0.1)
    divRec = append(divRec, goDivision)
    inRec = append(inRec, input)
    sigRec = append(sigRec, totSig)
    sig=go.res[go.res$p.adj < 0.1,]
    sigOut=paste( c('./resultsFiles/', goDivision, input, '_sigGos.tsv'), collapse='')
    if (nrow(sig)>0){
      sig %>% 
        write_tsv(file=sigOut)
    }
  }
}
res = tibble('goDivision'=divRec,
             'input'=inRec,
             'nSig'=sigRec)
res %>% 
  write_tsv(path='./resultsFiles/gomwu_results_summary.tsv')


# PLOT A SINGLE DATASET ---------------------------------------------------
input='brown_moduleInput.csv'; goDivision = "MF"

pdf("brown_MF.pdf")
gomwuPlot(input,goAnnotations,goDivision,
          absValue=0.00001,  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
          level1=.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
          level2=.05, # FDR cutoff to print in regular (not italic) font.
          level3=.01, # FDR cutoff to print in large bold font.
          txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
          treeHeight=0.5, # height of the hierarchical clustering tree
          #colors=c(COL,COL,COL,COL) # these are default colors, un-remar and change if needed
)
dev.off()

# plot results for each module with any significant enrichment -------------------
sig = res %>% 
  filter(nSig>=2)

for (i in 1:nrow(sig)){
  row=sig[i,]
  goDivision=row['goDivision']
  input = row['input']
  figFileName = paste('./resultsFiles/', sep='', paste(paste(goDivision, input, sep='_'), 'tree.pdf', sep='_'))
  pdf(figFileName)
  gomwuPlot(input,goAnnotations,goDivision,
            absValue=0.00001,  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
            level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
            level2=0.01, # FDR cutoff to print in regular (not italic) font.
            level3=0.001, # FDR cutoff to print in large bold font.
            txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
            treeHeight=0.5, # height of the hierarchical clustering tree
            #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
  )
  dev.off()
}

