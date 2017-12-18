# NIPA a robust set of tools for analyis of gene lists. 

##############################################################
# uncommmment to install any required packages. 
##############################################################
#biocLite("biomaRt")
#biocLite("gage")
#biocLite("pathview")
#biocLite("gageData")
#biocLite("ggplot2")
#biocLite("stringr")
#biocLite("dplyr")
#biocLite("RamiGO")
##############################################################

source("https://www.bioconductor.org/biocLite.R")
library(biomaRt)
library(pathview)
library(gage)
library(gageData)
library(ggplot2)
library(stringr)
library(dplyr)
library(RamiGO)


###############################################################################
## Input Variables -- USER TO CHANGE [START]
## Check all or may fail.
###############################################################################
goi.column = 1                # if results are from analysis and are a column of a larger table give input column else will assume is column 1 or a single column assumes tab delimited
goi.header = "no"             # "yes" or "no" if header on file 

species = "mouse"             #currently one of "mouse", "human", "rat", "pig", "zebrafish, cow, fly, sheep", 

# colour pathways by expression fold change?
keggFC = "no"                 # yes or no. will colour enriched KEGG pathways by FC data [specify column below]
keggFC.col = 14               # if keggFC = yes specify column of input table with FC values  assumes tab delimited


id.type = "ENSG"              # one of
                              # "ENSG" (ensembl gene),
                              # "ENST" (ensembl trasncript),
                              # "ENSP" (ensembl peptide),
                              # "Entrez"
                              # "Uniprot" (UniProt/SwissProt Accession)
                              # "Unigene"
                              # "Refseq_mrna" (RefSeq mRNA [e.g. NM_001195597])
                              # "Refseq_peptide" (RefSeq Protein ID [e.g. NP_001005353])
                              # "hgnc" (HGNC ID [e.g. LIS1])

# set variables for hypergeometric cutoff enrichment qval less than this and with greater or equal to minimum number of genes in pathway or GO term will be drawn
kegg.qval.cutoff = 0.05
GO.cutoff = 0.05              # qvalue cutoff
min.genes.cutoff = 2

# change below to determine which test to conduct.
doGO = "yes"                  # yes or no.       Run hypergeometric test to find enriched GO terms in BP, MF and CC category
doKEGG = "yes"                 # yes or no.     Run hypergeometric test to find and plot enriched KEGG pathways and visualise using PathView



###############################################################################
## Input Variables -- USER TO CHANGE [END]
###############################################################################














##############################################################################
# Dont alter below this line
##############################################################################
goi.list <- file.choose()
this.dir <- dirname(goi.list)
setwd(this.dir)

outfile.prefix <- goi.list # prefix attached to output files. 
###############################################################################
## set variables based on species given 
###############################################################################
if (species == "sheep")
  {
  ensembl.spp <- "oaries_gene_ensembl"
  species.kegg.code = "oas"
  }
 
if (species == "fly")
  {
  ensembl.spp <- "dmelanogaster_gene_ensembl"
  species.kegg.code = "dme"
  }

if (species == "mouse")
  {
  ensembl.spp <- "mmusculus_gene_ensembl"
  species.kegg.code = "mmu"
  }

if (species == "human")
  {
  ensembl.spp <- "hsapiens_gene_ensembl"
  species.kegg.code = "hsa"
  }

if (species == "rat")
  {
  ensembl.spp <- "rnorvegicus_gene_ensembl"
  species.kegg.code = "rno"
  }

if (species == "pig")
  {
  ensembl.spp <- "sscrofa_gene_ensembl"
  species.kegg.code = "ssc"
  }

if (species == "zebrafish")
  {
  ensembl.spp <- "drerio_gene_ensembl"
  species.kegg.code = "dre"
  }

if (species == "cow")
  {
  ensembl.spp <- "btaurus_gene_ensembl"
  species.kegg.code = "bta"
  }

##############################################################################
# Build kegg sets 
##############################################################################
kegg.gsets.spp <- kegg.gsets(species = species.kegg.code, id.type = "kegg") 
kegg.sets.test <- kegg.gsets.spp$kg.sets
kegg.sets.spp = kegg.gsets.spp$sigmet.idx


##############################################################################
# Get Data
##############################################################################

if (goi.header == "yes") {my.data.in <- read.table(goi.list,sep='\t',header = TRUE, quote = "")}
if (goi.header == "no") {my.data.in <- read.table(goi.list,sep='\t',header = FALSE, quote = "")}
myInterestingGenes <- as.vector(unlist(my.data.in[goi.column]))
myInterestingGenes <- unique(myInterestingGenes)




##############################################################################
# ADD entrez ids using ensembl for KEGG and match to gene input list
# ADD GO using ensembl
##############################################################################

if (id.type =="ENSG") {id.lookup = 'ensembl_gene_id'}
if (id.type =="ENSP") {id.lookup = 'ensembl_peptide_id'}
if (id.type =="ENST") {id.lookup = 'ensembl_transcript_id'}
if (id.type == "Entrez") {id.lookup = 'entrezgene'}
if (id.type == "Refseq_mrna") {id.lookup = 'refseq_mrna'}
if (id.type == "Refseq_peptide") {id.lookup = 'refseq_peptide'}
if (id.type == "Unigene") {id.lookup ='unigene'}
if (id.type == "Uniprot") {id.lookup = 'uniprotswissprot'}
if (id.type == "hgnc") {id.lookup = 'hgnc_symbol'}


ensembl = useEnsembl(biomart="ensembl", dataset=ensembl.spp)

# Get Gene Id info from ensembl 
all.genes <- getBM(attributes=c(id.lookup, 'entrezgene', 'external_gene_name'), mart = ensembl)
colnames(all.genes) <- c("ID","Entrez","Name")
all.genes.entrez <- na.omit(all.genes)
all.genes.entrez <- all.genes.entrez[all.genes.entrez$ID!="",]
goi.entrez <-unique(as.character(all.genes.entrez[all.genes.entrez$ID %in% myInterestingGenes,2]))

universe.size <- as.numeric(length(unique(all.genes.entrez$ID)))

# include GO terms from Ensembl for later analysis.
all.genes.GO <- getBM(attributes=c(id.lookup,'go_id','name_1006','namespace_1003'), mart = ensembl)
colnames(all.genes.GO) <- c("ID","GO_ID","GO_Name","GO_component")
all.GO.lookup <- unique(all.genes.GO[c("GO_ID","GO_Name")])


# if keggFC = yes create foldchanges named list of log fold change values
if (keggFC == "yes")
{
  entrez.FC.match <- merge(all.genes.entrez,my.data.in,by.x="ID",by.y=names(my.data.in[goi.column]))
  foldchanges = as.numeric(as.character(unlist(entrez.FC.match[keggFC.col+3])))
  names(foldchanges) = as.character(entrez.FC.match$Entrez)
}



##############################################################################
# start report and set up variables to catch failing sections. 
##############################################################################
run.report = paste(outfile.prefix,"NIPA.report.txt",sep=".")

cat(c("-----------------------------------------------------------","The NIPA run has initiated: Any warnings will appear below.","-----------------------------------------------------------"), file=run.report, append=FALSE, sep = "\n")

if (length(goi.entrez)==0 )
{
  cat(c("The run has terminated","Conversion of gene/peptide list to entrez failed", "Are the IDs properly formatted or possibly too few IDs"),
      file=run.report, append=TRUE, sep='\n')
  stop("Run terminated, see NIPA.report.txt")
}

# set flags to capture failed sections. 
fail.GO.MF = 0
fail.GO.BP = 0
fail.GO.CC = 0
fail.KEGG = 0
stats.KEGG.fail = 0

##############################################################################
#
# part 1 GO analysis
# 
##############################################################################
##############################################################################
##########################################################

if (doGO == "yes")
{
  
  #################################################################################
  # filter GO by type. 
  #################################################################################
  BP.genes.GO <- all.genes.GO[all.genes.GO$GO_component == "biological_process", ]
  MF.genes.GO <- all.genes.GO[all.genes.GO$GO_component == "molecular_function", ]
  CC.genes.GO <- all.genes.GO[all.genes.GO$GO_component == "cellular_component", ]
  
  ########################################################################################################### 
  #  Biological Process test GO enrichment by hypergeometric test
  ########################################################################################################### 
  # generate table of counts per GO term 
  BP.genes.GO.GOI <- BP.genes.GO[BP.genes.GO$ID %in% myInterestingGenes, ]
  
  BP.genes.GO.table <- as.data.frame(BP.genes.GO %>% dplyr::group_by(GO_ID) %>% 
                                      dplyr::summarise(gene_ids = paste(ID, collapse=" ")) %>%
                                      dplyr::mutate(gene_count = str_count(gene_ids, " ")+1))
  BP.genes.GO.table.GOI <- as.data.frame(BP.genes.GO.GOI %>% dplyr::group_by(GO_ID) %>% 
                                       dplyr::summarise(gene_ids = paste(ID, collapse=" ")) %>%
                                       dplyr::mutate(gene_count = str_count(gene_ids, " ")+1))
  
  BP.genes.GO.merge <- merge(BP.genes.GO.table, BP.genes.GO.table.GOI, by="GO_ID",all.y=TRUE)
  BP.genes.GO.merge <- merge(BP.genes.GO.merge,all.GO.lookup, by="GO_ID", all.x=TRUE)
  colnames(BP.genes.GO.merge) <- c("GO_ID","ALL.gene_ids","ALL.gene_count","GOI.gene_ids","GOI.gene_count","GO_Name")
  BP.genes.GO.merge <- BP.genes.GO.merge[-c(2)] # remove all gene names as not needed
  
  
  
  # do for each pathway in list and generate table of pathways passing cut off after FDR qvalue calculation
  working.GO.BP <- unique( BP.genes.GO.merge$GO_ID)
  GO.BP.hypergeometric.results <- data.frame("GO"= character(0),"p.val"= numeric(0),"FDR q.val"= numeric(0))
  
  
  for (i in 1:nrow(BP.genes.GO.merge))
  {
    current.GO.BP = BP.genes.GO.merge[i,]
    
    sample_success = as.numeric(current.GO.BP$GOI.gene_count) #  goi count in GO term 
    population_success = as.numeric(current.GO.BP$ALL.gene_count) # all gene count in GO term
    population_not_success = (universe.size-population_success)
    sample_size = as.numeric(length(unique(myInterestingGenes)))
    
    
    pval <- phyper(sample_success,population_success,population_not_success,sample_size, lower.tail=FALSE,log.p=FALSE)
    qval <- p.adjust(pval, method = "fdr", n = nrow(BP.genes.GO.merge))
    working.results <- cbind(current.GO.BP,pval,qval)
    GO.BP.hypergeometric.results <- rbind(GO.BP.hypergeometric.results,working.results)
  }
  BP.table.out = paste(outfile.prefix,"GO.BP.table",sep=".")
  GO.BP.hypergeometric.results <- GO.BP.hypergeometric.results[order(GO.BP.hypergeometric.results$pval),] # order by Pvalue
  
  #check for significant results
  GO.BP.hypergeometric.results.sig <- GO.BP.hypergeometric.results[GO.BP.hypergeometric.results$qval <= GO.cutoff & GO.BP.hypergeometric.results$GOI.gene_count >= min.genes.cutoff, ]
  
  if (nrow(GO.BP.hypergeometric.results.sig)==0 )
  {
    fail.GO.BP = 1
    cat(c("GO Biological Process search identified no enriched terms passing cutoffs: Probably too few IDs"),
        file=run.report, append=TRUE, sep='\n')
  }
  if (fail.GO.BP !=1)
  {
    write.table(GO.BP.hypergeometric.results.sig, file=BP.table.out, row.names = FALSE, col.names=TRUE,sep = '\t', quote=FALSE)
    top.result.BP <- head(GO.BP.hypergeometric.results.sig,10)
    top.result.BP <- top.result.BP[order(top.result.BP$pval),]
    top.result.BP$GO_Name <- as.factor(top.result.BP$GO_Name)
    top.result.BP$GO_Name <- factor(top.result.BP$GO_Name, levels = top.result.BP$GO_Name)
    
    if (nrow(top.result.BP) > 0)
    {
      current.min.pval <- min(top.result.BP$pval[top.result.BP$pval > 0])
      if (current.min.pval == Inf) {current.min.pval <- 1e-10}
      top.result.BP$pval[top.result.BP$pval == 0 ] <- current.min.pval  # catches any where p value = 0 and makes it equal to smallest p value.
      
      max.y.plot = 1.2*(max(-log10(top.result.BP$pval)))
      sig.BP.plot <-
        ggplot(data = top.result.BP,
               aes(x = as.factor(GO_Name), y = -log10(top.result.BP$pval),
                   size = GOI.gene_count))+
        geom_point() +
        scale_size_continuous(range = c(4,18), name="Gene count")+
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30))+
        geom_hline(yintercept=1.30103,lty=2, color="grey") + # equivalent of p = 0.05
        geom_hline(yintercept=2,lty=4, color="grey") + # equivalent of p = 0.01
        geom_hline(yintercept=3,lty=3, color="grey") + # equivalent of p = 0.001  
        coord_flip()+
        geom_point(stat = "identity",colour="royalblue4") +
        theme_bw() +
        theme(axis.text.x = element_text(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white")) +
        ylim(-0.5,max.y.plot)+
        xlab("") +
        ylab("Enrichment (-log10 pvalue)")
      
      
      BP.plot.out = paste(outfile.prefix,"GO.BP.Significant.enrichment.plot.pdf",sep=".")
      pdf(BP.plot.out)
      print(sig.BP.plot)
      dev.off()
      
      #### plot SVG of Directed Acyclic graph of 15 most signiifcnat GO terms.
      GO.BP.top.DAG <- paste(outfile.prefix,"GO.BP.top.DAG",sep='.')
      svgRes <- getAmigoTree(top.result.BP$GO_ID, color="red", pvalues =top.result.BP$pval, filename=GO.BP.top.DAG, picType="svg", saveResult=TRUE)
    }
  }
  
  ########################################################################################################### 
  #  Molecular Function test GO enrichment by hypergeometric test
  ########################################################################################################### 
  # generate table of counts per GO term 
  MF.genes.GO.GOI <- MF.genes.GO[MF.genes.GO$ID %in% myInterestingGenes, ]
  
  MF.genes.GO.table <- as.data.frame(MF.genes.GO %>% dplyr::group_by(GO_ID) %>% 
                                       dplyr::summarise(gene_ids = paste(ID, collapse=" ")) %>%
                                       dplyr::mutate(gene_count = str_count(gene_ids, " ")+1))
  MF.genes.GO.table.GOI <- as.data.frame(MF.genes.GO.GOI %>% dplyr::group_by(GO_ID) %>% 
                                           dplyr::summarise(gene_ids = paste(ID, collapse=" ")) %>%
                                           dplyr::mutate(gene_count = str_count(gene_ids, " ")+1))
  
  MF.genes.GO.merge <- merge(MF.genes.GO.table, MF.genes.GO.table.GOI, by="GO_ID",all.y=TRUE)
  MF.genes.GO.merge <- merge(MF.genes.GO.merge,all.GO.lookup, by="GO_ID", all.x=TRUE)
  colnames(MF.genes.GO.merge) <- c("GO_ID","ALL.gene_ids","ALL.gene_count","GOI.gene_ids","GOI.gene_count","GO_Name")
  MF.genes.GO.merge <- MF.genes.GO.merge[-c(2)] # remove all gene names as not needed
  
  # do for each pathway in list and generate table of pathways passing cut off after FDR qvalue calculation
  working.GO.MF <- unique( MF.genes.GO.merge$GO_ID)
  GO.MF.hypergeometric.results <- data.frame("GO"= character(0),"p.val"= numeric(0),"FDR q.val"= numeric(0))
  
  
  for (i in 1:nrow(MF.genes.GO.merge))
  {
    current.GO.MF = MF.genes.GO.merge[i,]
    sample_success = as.numeric(current.GO.MF$GOI.gene_count) #  goi count in GO term 
    population_success = as.numeric(current.GO.MF$ALL.gene_count) # all gene count in GO term
    population_not_success = (universe.size-population_success)
    sample_size = as.numeric(length(unique(myInterestingGenes)))
    
    
    pval <- phyper(sample_success,population_success,population_not_success,sample_size, lower.tail=FALSE,log.p=FALSE)
    qval <- p.adjust(pval, method = "fdr", n = nrow(MF.genes.GO.merge))
    working.results <- cbind(current.GO.MF,pval,qval)
    GO.MF.hypergeometric.results <- rbind(GO.MF.hypergeometric.results,working.results)
  }
  MF.table.out = paste(outfile.prefix,"GO.MF.table",sep=".")
  GO.MF.hypergeometric.results <- GO.MF.hypergeometric.results[order(GO.MF.hypergeometric.results$pval),] # order by Pvalue
  
  #check for significant results
  GO.MF.hypergeometric.results.sig <- GO.MF.hypergeometric.results[GO.MF.hypergeometric.results$qval <= GO.cutoff & GO.MF.hypergeometric.results$GOI.gene_count >= min.genes.cutoff, ]
  
  if (nrow(GO.MF.hypergeometric.results.sig)==0 )
  {
    fail.GO.MF = 1
    cat(c("GO Molecular Function search identified no enriched terms passing cutoffs: Probably too few IDs"),
        file=run.report, append=TRUE, sep='\n')
  }
  if (fail.GO.MF !=1)
  {
    write.table(GO.MF.hypergeometric.results.sig, file=MF.table.out, row.names = FALSE, col.names=TRUE,sep = '\t', quote=FALSE)
    top.result.MF <- head(GO.MF.hypergeometric.results.sig,10)
    top.result.MF <- top.result.MF[order(top.result.MF$pval),]
    top.result.MF$GO_Name <- as.factor(top.result.MF$GO_Name)
    top.result.MF$GO_Name <- factor(top.result.MF$GO_Name, levels = top.result.MF$GO_Name)
    
    if (nrow(top.result.MF) > 0)
    {
      current.min.pval <- min(top.result.MF$pval[top.result.MF$pval > 0])
      if (current.min.pval == Inf) {current.min.pval <- 1e-10}
      top.result.MF$pval[top.result.MF$pval == 0 ] <- current.min.pval  # catches any where p value = 0 and makes it equal to smallest p value.
      
      max.y.plot = 1.2*(max(-log10(top.result.MF$pval)))
      sig.MF.plot <-
        ggplot(data = top.result.MF,
               aes(x = as.factor(GO_Name), y = -log10(top.result.MF$pval),
                   size = GOI.gene_count))+
        geom_point() +
        scale_size_continuous(range = c(4,18), name="Gene count")+
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30))+
        geom_hline(yintercept=1.30103,lty=2, color="grey") + # equivalent of p = 0.05
        geom_hline(yintercept=2,lty=4, color="grey") + # equivalent of p = 0.01
        geom_hline(yintercept=3,lty=3, color="grey") + # equivalent of p = 0.001  
        coord_flip()+
        geom_point(stat = "identity",colour="royalblue4") +
        theme_bw() +
        theme(axis.text.x = element_text(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white")) +
        ylim(-0.5,max.y.plot)+
        xlab("") +
        ylab("Enrichment (-log10 pvalue)")
      
      
      MF.plot.out = paste(outfile.prefix,"GO.MF.Significant.enrichment.plot.pdf",sep=".")
      pdf(MF.plot.out)
      print(sig.MF.plot)
      dev.off()
      
      #### plot SVG of Directed Acyclic graph of 15 most signiifcnat GO terms.
      GO.MF.top.DAG <- paste(outfile.prefix,"GO.MF.top.DAG",sep='.')
      svgRes <- getAmigoTree(top.result.MF$GO_ID, color="red", pvalues =top.result.MF$pval, filename=GO.MF.top.DAG, picType="svg", saveResult=TRUE)
    }
  }
  
  ########################################################################################################### 
  #  Cellular Component test GO enrichment by hypergeometric test
  ########################################################################################################### 
  # generate table of counts per GO term 
  CC.genes.GO.GOI <- CC.genes.GO[CC.genes.GO$ID %in% myInterestingGenes, ]
  
  CC.genes.GO.table <- as.data.frame(CC.genes.GO %>% dplyr::group_by(GO_ID) %>% 
                                       dplyr::summarise(gene_ids = paste(ID, collapse=" ")) %>%
                                       dplyr::mutate(gene_count = str_count(gene_ids, " ")+1))
  CC.genes.GO.table.GOI <- as.data.frame(CC.genes.GO.GOI %>% dplyr::group_by(GO_ID) %>% 
                                           dplyr::summarise(gene_ids = paste(ID, collapse=" ")) %>%
                                           dplyr::mutate(gene_count = str_count(gene_ids, " ")+1))
  
  CC.genes.GO.merge <- merge(CC.genes.GO.table, CC.genes.GO.table.GOI, by="GO_ID",all.y=TRUE)
  CC.genes.GO.merge <- merge(CC.genes.GO.merge,all.GO.lookup, by="GO_ID", all.x=TRUE)
  colnames(CC.genes.GO.merge) <- c("GO_ID","ALL.gene_ids","ALL.gene_count","GOI.gene_ids","GOI.gene_count","GO_Name")
  CC.genes.GO.merge <- CC.genes.GO.merge[-c(2)] # remove all gene names as not needed
  
  # do for each pathway in list and generate table of pathways passing cut off after FDR qvalue calculation
  working.GO.CC <- unique( CC.genes.GO.merge$GO_ID)
  GO.CC.hypergeometric.results <- data.frame("GO"= character(0),"p.val"= numeric(0),"FDR q.val"= numeric(0))
  
  
  for (i in 1:nrow(CC.genes.GO.merge))
  {
    current.GO.CC = CC.genes.GO.merge[i,]
    sample_success = as.numeric(current.GO.CC$GOI.gene_count) #  goi count in GO term 
    population_success = as.numeric(current.GO.CC$ALL.gene_count) # all gene count in GO term
    population_not_success = (universe.size-population_success)
    sample_size = as.numeric(length(unique(myInterestingGenes)))
    
    pval <- phyper(sample_success,population_success,population_not_success,sample_size, lower.tail=FALSE,log.p=FALSE)
    qval <- p.adjust(pval, method = "fdr", n = nrow(CC.genes.GO.merge))
    working.results <- cbind(current.GO.CC,pval,qval)
    GO.CC.hypergeometric.results <- rbind(GO.CC.hypergeometric.results,working.results)
  }
  CC.table.out = paste(outfile.prefix,"GO.CC.table",sep=".")
  GO.CC.hypergeometric.results <- GO.CC.hypergeometric.results[order(GO.CC.hypergeometric.results$pval),] # order by Pvalue
  
  #check for significant results
  GO.CC.hypergeometric.results.sig <- GO.CC.hypergeometric.results[GO.CC.hypergeometric.results$qval <= GO.cutoff & GO.CC.hypergeometric.results$GOI.gene_count >= min.genes.cutoff, ]
  
  if (nrow(GO.CC.hypergeometric.results.sig)==0 )
  {
    fail.GO.CC = 1
    cat(c("GO Cellular Compartment search identified no enriched terms passing cutoffs: Probably too few IDs"),
        file=run.report, append=TRUE, sep='\n')
  }
  if (fail.GO.CC !=1)
  {
    write.table(GO.CC.hypergeometric.results.sig, file=CC.table.out, row.names = FALSE, col.names=TRUE,sep = '\t', quote=FALSE)
    top.result.CC <- head(GO.CC.hypergeometric.results.sig,10)
    top.result.CC <- top.result.CC[order(top.result.CC$pval),]
    top.result.CC$GO_Name <- as.factor(top.result.CC$GO_Name)
    top.result.CC$GO_Name <- factor(top.result.CC$GO_Name, levels = top.result.CC$GO_Name)
    
    if (nrow(top.result.CC) > 0)
    {
      current.min.pval <- min(top.result.CC$pval[top.result.CC$pval > 0])
      if (current.min.pval == Inf) {current.min.pval <- 1e-10}
      top.result.CC$pval[top.result.CC$pval == 0 ] <- current.min.pval  # catches any where p value = 0 and makes it equal to smallest p value.
      
      max.y.plot = 1.2*(max(-log10(top.result.CC$pval)))
      sig.CC.plot <-
        ggplot(data = top.result.CC,
               aes(x = as.factor(GO_Name), y = -log10(top.result.CC$pval),
                   size = GOI.gene_count))+
        geom_point() +
        scale_size_continuous(range = c(4,18), name="Gene count")+
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30))+
        geom_hline(yintercept=1.30103,lty=2, color="grey") + # equivalent of p = 0.05
        geom_hline(yintercept=2,lty=4, color="grey") + # equivalent of p = 0.01
        geom_hline(yintercept=3,lty=3, color="grey") + # equivalent of p = 0.001  
        coord_flip()+
        geom_point(stat = "identity",colour="royalblue4") +
        theme_bw() +
        theme(axis.text.x = element_text(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white")) +
        ylim(-0.5,max.y.plot)+
        xlab("") +
        ylab("Enrichment (-log10 pvalue)")
      
      
      CC.plot.out = paste(outfile.prefix,"GO.CC.Significant.enrichment.plot.pdf",sep=".")
      pdf(CC.plot.out)
      print(sig.CC.plot)
      dev.off()
      
      #### plot SVG of Directed Acyclic graph of 15 most signiifcnat GO terms.
      GO.CC.top.DAG <- paste(outfile.prefix,"GO.CC.top.DAG",sep='.')
      svgRes <- getAmigoTree(top.result.CC$GO_ID, color="red", pvalues =top.result.CC$pval, filename=GO.CC.top.DAG, picType="svg", saveResult=TRUE)
    }
  }
  
  
  
  
}
####
##############################################################################
##############################################################################
#
# part 2 KEGG analysis
# 
##############################################################################
##############################################################################

if (doKEGG == "yes")
{
  pathview.goi.entrez <- rep.int(1, length(goi.entrez))
  
  names(pathview.goi.entrez) = goi.entrez 
  
  keggres = gage(pathview.goi.entrez, gsets=kegg.sets.test, same.dir=TRUE) # determine kegg membership of all genes. 
  
  keggres.pathways <- as.data.frame(keggres)
  keggres.pathways.out <- keggres.pathways[keggres.pathways$greater.set.size > 0, ]
  
  if (nrow(keggres.pathways.out) ==0)
  {
    fail.KEGG = 1
    cat(c("KEGG analysis identified no enriched pathways: Probably too few IDs"),
        file=run.report, append=TRUE, sep='\n')
  }
  
  if (fail.KEGG ==0)
  {
    keggres.pathways.out$KEGGpathways <- rownames(keggres.pathways.out)
    matching.kegg.sets.spp <- kegg.sets.test[c(keggres.pathways.out$KEGGpathways)] # named list of matched pathways
    matching.kegg.sets.spp.total.size <- lengths(matching.kegg.sets.spp, use.names = TRUE) # named list of the number of total number of genes in matched pathway.
    
   
    
    matching.kegg.sets.spp.df <- as.data.frame(unlist(matching.kegg.sets.spp, use.names = TRUE))
    matching.kegg.sets.spp.df$kegg.id <-  gsub("\\d+$", "", rownames(matching.kegg.sets.spp.df)) 
    row.names(matching.kegg.sets.spp.df) <- NULL
    colnames(matching.kegg.sets.spp.df) <- c("entrez.id","kegg.id")
    
    # make subset of matching.kegg.sets.spp.df with just genes of interest in it
    goi.matching.kegg.sets.spp.df <- matching.kegg.sets.spp.df[matching.kegg.sets.spp.df$entrez.id %in% goi.entrez, ]
    
    #################################################################################################################
    # Stats details
    #################################################################################################################
    # for each pathway with > 0 goi in it, conduct a hypergeometric test using phyper
    # phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
    # x, q vector of quantiles representing the number of white balls drawn
    # without replacement from an urn which contains both black and white
    # balls.
    # m the number of white balls in the urn.
    # n the number of black balls in the urn.
    # k the number of balls drawn from the urn.
    # if 
    # pop size : 5260 # total number of entrez gene in all pathways
    # sample size : 131 # total goi
    # Number of items in the pop that are classified as successes : 1998 # entrez in a particular pathway
    # Number of items in the sample that are classified as successes : 62 # goi in a particular pathway
    # 
    # phyper(62,1998,5260-1998,131)
    # e.g pathway 100 genes 10 are in goi list of size 400 universe = 20,000
    # phyper(1,100,20000-100,400, lower.tail=FALSE) = 0.597 = probability of finding this many or greater goi in pathway 
    # phyper(80,100,20000-100,400, lower.tail=FALSE) = 4.603708e-122 = probability of finding this many or greater goi in pathway 
    #################################################################################################################
    #
    #################################################################################################################
    
    
    # run phyper for each pathway in list and generate table of pathways passing cut off after FDR qvalue calculation
    working.pathways <- unique(matching.kegg.sets.spp.df$kegg.id)
    pathways.hypergeometric.results <- data.frame("Pathway"= character(0),"p.val"= numeric(0),"FDR q.val"= numeric(0),"ID"= character(0), "entrez.ids"= numeric(0), "external.ids"= character(0), "goi.count"= numeric(0), "all.count"= numeric(0))
    
    for (i in 1:length(working.pathways)){
      current.pathway = working.pathways[i]
      
      sample_success = as.numeric(nrow(goi.matching.kegg.sets.spp.df[goi.matching.kegg.sets.spp.df$kegg.id == current.pathway, ])) 
      population_success = as.numeric(nrow(matching.kegg.sets.spp.df[matching.kegg.sets.spp.df$kegg.id == current.pathway, ]))
      population_not_success = (universe.size-population_success)
      sample_size = as.numeric(length(unique(myInterestingGenes)))
      
      
      pval <- phyper(sample_success,population_success,population_not_success,sample_size, lower.tail=FALSE,log.p=FALSE)
      qval <- p.adjust(pval, method = "fdr", n = length(working.pathways))
      
      current.goi <- goi.matching.kegg.sets.spp.df[goi.matching.kegg.sets.spp.df$kegg.id == current.pathway, ]
      current.goi <- current.goi[1]
      current.goi.entrez.ids <- as.numeric(as.character(current.goi$entrez.id))
      
      current.goi.ens <- all.genes.entrez[all.genes.entrez$Entrez %in% current.goi.entrez.ids ,]
      
      current.goi.ens.ids <- unique(current.goi.ens$ID)
      current.goi.ext.ids <- unique(current.goi.ens$Name)
      
      
      current.goi.ens.ids <- paste(current.goi.ens.ids, collapse=", ")
      current.goi.entrez.ids <- paste(current.goi.entrez.ids, collapse=", ")
      current.goi.ext.ids <- paste(current.goi.ext.ids, collapse=", ")
      
      current.out <- as.data.frame(cbind(current.pathway,pval,qval,current.goi.ens.ids,current.goi.entrez.ids,current.goi.ext.ids,sample_success,population_success))
      
      pathways.hypergeometric.results <- rbind(pathways.hypergeometric.results, current.out)
      }  

    
    colnames(pathways.hypergeometric.results) <- c("Pathway","p.val","FDR q.val","GOI.ids","Entrez.ids","External.ids","goi.count","All.genes.in.pathway.count")

    
    # make FDR q.val and goi count numeric and sort 
    pathways.hypergeometric.results$`FDR q.val` <- as.numeric(as.character(pathways.hypergeometric.results$`FDR q.val`))
    pathways.hypergeometric.results$goi.count <- as.numeric(as.character(pathways.hypergeometric.results$goi.count))
    pathways.hypergeometric.results <-  pathways.hypergeometric.results[with(pathways.hypergeometric.results, order(pathways.hypergeometric.results$`FDR q.val`)), ]
    
    
    kegg.table.out = paste(outfile.prefix,"kegg.pathway.enrichment.table",sep=".")
    write.table( pathways.hypergeometric.results,file=kegg.table.out, row.names = FALSE, col.names = TRUE, quote = FALSE, sep ='\t')
    
    pathways.hypergeometric.results.sig <- pathways.hypergeometric.results[pathways.hypergeometric.results$`FDR q.val` < kegg.qval.cutoff & pathways.hypergeometric.results$goi.count >= min.genes.cutoff, ]
    
    ##############################################################################################  
    # draw plot of enriched pathways
    ############################################################################################## 
   
    if (nrow(pathways.hypergeometric.results.sig)>0)
    {
      
      kegg.sig.table.out = paste(outfile.prefix,"kegg.pathway.significant.enrichment.table",sep=".")
      write.table( pathways.hypergeometric.results.sig ,file=kegg.sig.table.out, row.names = FALSE, col.names = TRUE, quote = FALSE, sep ='\t')
      
      pathways.hypergeometric.results.sig$p.val <- as.numeric(as.character(pathways.hypergeometric.results.sig$p.val))
     
      ## replace FDR qval of 0 with v small number to avoid infinite values. 
      pathways.hypergeometric.results.sig <- within(pathways.hypergeometric.results.sig, `FDR q.val`[`FDR q.val` == 0] <- 1e-10)
      
      pathways.hypergeometric.results.sig <-  pathways.hypergeometric.results.sig[with(pathways.hypergeometric.results.sig, order(pathways.hypergeometric.results.sig$`FDR q.val`)), ]
      
      pathways.hypergeometric.results.sig$goi.count <- as.numeric(as.character(pathways.hypergeometric.results.sig$goi.count))
      pathways.hypergeometric.results.sig <-  pathways.hypergeometric.results.sig[with(pathways.hypergeometric.results.sig, order(pathways.hypergeometric.results.sig$`FDR q.val`)), ]
      
      
      top.pathways.hypergeometric.results.sig <- head(pathways.hypergeometric.results.sig,10)
      top.pathways.hypergeometric.results.sig$Pathway <- factor(top.pathways.hypergeometric.results.sig$Pathway, levels = top.pathways.hypergeometric.results.sig$Pathway)
      
      
      max.y.plot = 1.2*(max(-log10(top.pathways.hypergeometric.results.sig$`FDR q.val`)))
      sig.kegg.plot <-
        ggplot(data = top.pathways.hypergeometric.results.sig,
               aes(x = as.factor(Pathway), y = -log10(top.pathways.hypergeometric.results.sig$`FDR q.val`),
                   size = goi.count))+
        geom_point() +
        scale_size_continuous(range = c(4,18), "Gene count")+
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30))+
        geom_hline(yintercept=1.30103,lty=2, color="grey") + # equivalent of p = 0.05
        geom_hline(yintercept=2,lty=4, color="grey") + # equivalent of p = 0.01
        geom_hline(yintercept=3,lty=3, color="grey") + # equivalent of p = 0.001  
        coord_flip()+
        geom_point(stat = "identity",colour="royalblue4") +
        theme_bw() +
        theme(axis.text.x = element_text(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_rect(fill = "white")) +
        ylim(-0.5,max.y.plot)+
        xlab("") +
        ylab("Enrichment (-log10 pvalue)")
      
      
      kegg.pdf.out = paste(outfile.prefix,"KEGG.Significant.enrichment.plot.pdf",sep=".")
      pdf(kegg.pdf.out)
      print(sig.kegg.plot)
      dev.off()
      stats.KEGG.fail = 1
    }
    
    

    
     
      
##############################################################################################  
# draw Pathview plots of top enriched KEGG pathways
##############################################################################################    
detach("package:dplyr") # to overcome occasional issues of pathview clashing with dplyr

    
    if (stats.KEGG.fail == 1)
    {
    top.pathways.hypergeometric.results.sig$Pathway <- as.character(top.pathways.hypergeometric.results.sig$Pathway)
    
    for (i in 1:nrow(top.pathways.hypergeometric.results.sig))
      {
      current.sig.pathway = top.pathways.hypergeometric.results.sig$Pathway[i]
      pid <- substr(current.sig.pathway, start=1, stop=8) # get kegg ids 
      num.pid <- substr(pid, start=4, stop=8) # get kegg ids

        if(num.pid != "01100") # avoid drawing entire metabolic pathway plot.
          {
            if (keggFC == "yes")
            {
            pathview(gene.data=foldchanges, pathway.id=pid, species=species.kegg.code)
            tmp.xml <- paste(pid,".xml",sep='')
            xml.tmp <- paste(this.dir,tmp.xml,sep='/')
            tmp.png <- paste(pid,".png",sep='')
            png.tmp <- paste(this.dir,tmp.png,sep='/')
            old.pathview = paste(pid,"pathview.png",sep=".")
            new.pathview = paste(outfile.prefix,old.pathview,sep=".")
            file.remove(xml.tmp)
            file.remove(png.tmp)
            file.rename(old.pathview, new.pathview)
            }

            if (keggFC == "no")
            {
            pathview(gene.data=pathview.goi.entrez, pathway.id=pid, species=species.kegg.code)
            tmp.xml <- paste(pid,".xml",sep='')
            xml.tmp <- paste(this.dir,tmp.xml,sep='/')
            tmp.png <- paste(pid,".png",sep='')
            png.tmp <- paste(this.dir,tmp.png,sep='/')
            old.pathview = paste(pid,"pathview.png",sep=".")
            new.pathview = paste(outfile.prefix,old.pathview,sep=".")
            file.remove(xml.tmp)
            file.remove(png.tmp)
            file.rename(old.pathview, new.pathview)
            
            }
        }
        }
    }
  


    if (stats.KEGG.fail == 0)
    {
      cat(c("KEGG analysis: no pathways pass statistical cutoffs"),
          file=run.report, append=TRUE, sep='\n')
    }
    
  }
}

