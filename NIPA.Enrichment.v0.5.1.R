# NIPA a robust set of tools for analyis of gene lists. 

##############################################################
# uncommmment to install any required packages. 
##############################################################
#biocLite("biomaRt")
#biocLite("GOstats")
#biocLite("ReactomePA")
#biocLite("gage")
#biocLite("pathview")
#biocLite("gageData")
#biocLite("ggplot2")
#biocLite("stringr")
#biocLite("dplyr")
#biocLite("RamiGO")
##############################################################

source("http://www.bioconductor.org/biocLite.R")
library(GOstats)
library(biomaRt)
library(pathview)
library(gage)
library(gageData)
library(ReactomePA)
library(ggplot2)
library(stringr)
library(dplyr)
library(RamiGO)

###############################################################################
## Input Variables -- USER TO CHANGE [START]
## Check all or may fail.
###############################################################################
goi.column = 1 # if results are from analysis and are a column of a larger table give input column else will assume is column 1 or a single column assumes tab delimited
goi.header = "no" # "yes" or "no" if header on file 

species = "human"   #currently one of "mouse", "human", "rat", "pig", "zebrafish, cow, fly", 
outfile.prefix <- "ADAC.test" # prefix attached to output files. 

# colour pathways by expression fold change?
keggFC = "no" # yes or no. will colour enriched KEGG pathways by FC data [specify column below]
keggFC.col = 14 # if keggFC = yes specify column of input table with FC values  assumes tab delimited


id.type = "ENSG"      # one of
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
kegg.qval.cutoff = 0.1
GO.cutoff = 0.05
min.genes.cutoff = 2

# change below to determine which test to conduct.
doGO = "yes" # yes or no.       Run GoStats hypergeometric test to find enriched GO terms in BP, MF and CC category
doReactome = "no" # yes or no. Run ReactomePA to find enriched pathways in Reactomedb -- BIT SLOWER
doKEGG = "yes" # yes or no.     Run hypergeometric test to find and plot enriched KEGG pathways and visualise using PathView


###############################################################################
## Input Variables -- USER TO CHANGE [END]
###############################################################################














##############################################################################
# Dont alter below this line
##############################################################################
goi.list <- file.choose()
this.dir <- dirname(goi.list)
setwd(this.dir)

###############################################################################
## set variables based on species given 
###############################################################################
if (species == "fly")
{
  biocLite("org.Dm.eg.db") # for drosophila
  library(org.Dm.eg.db)    
  ensembl.spp <- "dmelanogaster_gene_ensembl"
  species.ens.code = "Dm"
  species.kegg.code = "dme"
  kegg.data.code = "dm"
  reactome.spp = "fly" #one of "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly".
  kegg.gsets.spp <- kegg.gsets(species = "dme", id.type = "kegg")  
}

if (species == "mouse")
{
  biocLite("org.Mm.eg.db") # for Mouse
  library(org.Mm.eg.db)    
  ensembl.spp <- "mmusculus_gene_ensembl"
  species.ens.code = "Mm"
  species.kegg.code = "mmu"
  kegg.data.code = "mm"
  reactome.spp = "mouse" #one of "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly".
  kegg.gsets.spp <- kegg.gsets(species = "mmu", id.type = "kegg")  
}

if (species == "human")
{
  biocLite("org.Hs.eg.db") # for Human
  library(org.Hs.eg.db)    
  ensembl.spp <- "hsapiens_gene_ensembl"
  species.ens.code = "Hs"
  species.kegg.code = "hsa"
  kegg.data.code = "hsa"
  reactome.spp = "human" #one of "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly".
  kegg.gsets.spp <- kegg.gsets(species = "hsa", id.type = "kegg")  
}

if (species == "rat")
{
  biocLite("org.Rn.eg.db") # for Rat
  library(org.Rn.eg.db)    
  ensembl.spp <- "rnorvegicus_gene_ensembl"
  species.ens.code = "Rn"
  species.kegg.code = "rno"
  kegg.data.code = "rno"
  reactome.spp = "rat" #one of "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly".
  kegg.gsets.spp <- kegg.gsets(species = "rno", id.type = "kegg")  
}

if (species == "pig")
{
  biocLite("org.Ss.eg.db") # for Pig
  library(org.Ss.eg.db)    
  ensembl.spp <- "sscrofa_gene_ensembl"
  species.ens.code = "Ss"
  species.kegg.code = "ssc"
  kegg.data.code = "ssc"
  #reactome.spp = "rat" #one of "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly".
  doReactome = "no"
  kegg.gsets.spp <- kegg.gsets(species = "ssc", id.type = "kegg")  
}

if (species == "zebrafish")
{
  biocLite("org.Dr.eg.db") # for Zebrafish
  library(org.Dr.eg.db)    
  ensembl.spp <- "drerio_gene_ensembl"
  species.ens.code = "Dr"
  species.kegg.code = "dre"
  kegg.data.code = "dre"
  #reactome.spp = "rat" #one of "human", "rat", "mouse", "celegans", "yeast", "zebrafish", "fly".
  doReactome = "no"
  kegg.gsets.spp <- kegg.gsets(species = "dre", id.type = "kegg")  
}

if (species == "cow")
{
  biocLite("org.Bt.eg.db") # for Cow
  library(org.Bt.eg.db)    
  ensembl.spp <- "btaurus_gene_ensembl"
  species.ens.code = "Bt"
  species.kegg.code = "bta"
  kegg.data.code = "bta"
  doReactome = "no"
  kegg.gsets.spp <- kegg.gsets(species = "bta", id.type = "kegg")  
}


##############################################################################
# Build kegg sets 
##############################################################################
kegg.sets.test <- kegg.gsets.spp$kg.sets
kegg.sets.spp = kegg.gsets.spp$sigmet.idx

##############################################################################
# Get Data
##############################################################################

if (goi.header == "yes") {my.data.in <- read.table(goi.list,sep='\t',header = TRUE, quote = "")}
if (goi.header == "no") {my.data.in <- read.table(goi.list,sep='\t',header = FALSE, quote = "")}
myInterestingGenes <- as.vector(unlist(my.data.in[goi.column]))
myInterestingGenes <- unique(myInterestingGenes)

species.db <- paste("org",species.ens.code,"eg.db",sep=".")
ensembl = useEnsembl(biomart="ensembl", dataset=ensembl.spp)



##############################################################################
# Convert IDs to Entrez IDs and match to gene input list
##############################################################################

if (id.type =="ENSG")
{
  all.genes <- getBM(attributes=c('ensembl_gene_id', 'entrezgene', 'external_gene_name'), mart = ensembl)
  colnames(all.genes) <- c("ID","Entrez","Name")
  all.genes.entrez <- na.omit(all.genes)
  all.genes.entrez <- all.genes.entrez[all.genes.entrez$ID!="",]
  goi.entrez <-unique(as.character(all.genes.entrez[all.genes.entrez$ID %in% myInterestingGenes,2]))
  }

if (id.type =="ENSP")
{
  all.genes <- getBM(attributes=c('ensembl_peptide_id', 'entrezgene', 'external_gene_name'), mart = ensembl)
  colnames(all.genes) <- c("ID","Entrez","Name")
  all.genes.entrez <- na.omit(all.genes)
  all.genes.entrez <- all.genes.entrez[all.genes.entrez$ID!="",]
  goi.entrez <-unique(as.character(all.genes.entrez[all.genes.entrez$ID %in% myInterestingGenes,2]))
}

if (id.type =="ENST")
{
  all.genes <- getBM(attributes=c('ensembl_transcript_id', 'entrezgene', 'external_gene_name'), mart = ensembl)
  colnames(all.genes) <- c("ID","Entrez","Name")
  all.genes.entrez <- na.omit(all.genes)
  all.genes.entrez <- all.genes.entrez[all.genes.entrez$ID!="",]
  goi.entrez <-unique(as.character(all.genes.entrez[all.genes.entrez$ID %in% myInterestingGenes,2]))
}

if (id.type == "Entrez")
{
  all.genes <- getBM(attributes=c('entrezgene', 'entrezgene', 'external_gene_name'), mart = ensembl)
  colnames(all.genes) <- c("ID","Entrez","Name")
  all.genes.entrez <- na.omit(all.genes)
  all.genes.entrez <- all.genes.entrez[all.genes.entrez$ID!="",]
  goi.entrez <-unique(as.character(all.genes.entrez[all.genes.entrez$ID %in% myInterestingGenes,2]))
}

if (id.type == "Refseq_mrna")
{
  all.genes <- getBM(attributes=c('refseq_mrna', 'entrezgene', 'external_gene_name'), mart = ensembl)
  colnames(all.genes) <- c("ID","Entrez","Name")
  all.genes.entrez <- na.omit(all.genes)
  all.genes.entrez <- all.genes.entrez[all.genes.entrez$ID!="",]
  goi.entrez <-unique(as.character(all.genes.entrez[all.genes.entrez$ID %in% myInterestingGenes,2]))
}

if (id.type == "Refseq_peptide")
{
  all.genes <- getBM(attributes=c('refseq_peptide', 'entrezgene', 'external_gene_name'), mart = ensembl)
  colnames(all.genes) <- c("ID","Entrez","Name")
  all.genes.entrez <- na.omit(all.genes)
  all.genes.entrez <- all.genes.entrez[all.genes.entrez$ID!="",]
  goi.entrez <-unique(as.character(all.genes.entrez[all.genes.entrez$ID %in% myInterestingGenes,2]))
}

if (id.type == "Unigene")
{
  all.genes <- getBM(attributes=c('unigene', 'entrezgene', 'external_gene_name'), mart = ensembl)
  colnames(all.genes) <- c("ID","Entrez","Name")
  all.genes.entrez <- na.omit(all.genes)
  all.genes.entrez <- all.genes.entrez[all.genes.entrez$ID!="",]
  goi.entrez <-unique(as.character(all.genes.entrez[all.genes.entrez$ID %in% myInterestingGenes,2]))
}

if (id.type == "Uniprot")
{
  all.genes <- getBM(attributes=c('uniprotswissprot', 'entrezgene', 'external_gene_name'), mart = ensembl)
  colnames(all.genes) <- c("ID","Entrez","Name")
  all.genes.entrez <- na.omit(all.genes)
  all.genes.entrez <- all.genes.entrez[all.genes.entrez$ID!="",]
  goi.entrez <-unique(as.character(all.genes.entrez[all.genes.entrez$ID %in% myInterestingGenes,2]))
}

if (id.type == "hgnc")
{
  all.genes <- getBM(attributes=c('hgnc_symbol', 'entrezgene', 'external_gene_name'), mart = ensembl)
  colnames(all.genes) <- c("ID","Entrez","Name")
  all.genes.entrez <- na.omit(all.genes)
  all.genes.entrez <- all.genes.entrez[all.genes.entrez$ID!="",]
  goi.entrez <-unique(as.character(all.genes.entrez[all.genes.entrez$ID %in% myInterestingGenes,2]))
}

# if keggFC = yes create foldchanges named list of log fold change values
if (keggFC == "yes")
{
  entrez.FC.match <- merge(all.genes.entrez,my.data.in,by.x="ID",by.y=names(my.data.in[goi.column]))
  foldchanges = unlist(entrez.FC.match[keggFC.col+2])
  names(foldchanges) = entrez.FC.match$Entrez
}






##########################################################
# Set gene "universse" of all genes
universe <- unique(as.character(all.genes.entrez$Entrez))
##########################################################


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
fail.reactome = 0
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
  # Biological Process
  params.BP <- new('GOHyperGParams',
                   geneIds=goi.entrez,
                   universeGeneIds=universe,
                   ontology='BP',
                   pvalueCutoff=GO.cutoff,
                   conditional=F,
                   testDirection='over',
                   annotation=species.db
  )
  hgOver.BP <- hyperGTest(params.BP)
  result.BP <- summary(hgOver.BP)
  
  
  if (nrow(result.BP)==0 )
  {
    fail.GO.BP = 1
    cat(c("GO Biological process search identified no enriched terms","Probably too few IDs"),
        file=run.report, append=TRUE, sep='\n')
  }
  if (fail.GO.BP !=1)
  {
    result.BP <- result.BP[result.BP$Count >= min.genes.cutoff,] # filter those with < cut off count
    result.BP <- result.BP[order(-result.BP$OddsRatio),] # order by odds ratio
    
    top.result.BP <- head(result.BP,15)
    top.result.BP <- top.result.BP[order(top.result.BP$Pvalue),]
    top.result.BP$Term <- as.factor(top.result.BP$Term)
    top.result.BP$Term <- factor(top.result.BP$Term, levels = top.result.BP$Term)
    
    if (nrow(top.result.BP) > 0)
    {
      top.result.BP$Pvalue[top.result.BP$Pvalue == 0 ] <- 1e-10 # catches any where p value = 0
      max.y.plot = 1.2*(max(-log10(top.result.BP$Pvalue)))
      sig.BP.plot <-
        ggplot(data = top.result.BP,
               aes(x = as.factor(Term), y = -log10(top.result.BP$Pvalue),
                   size = Count))+
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
      
      
      BP.plot.out = paste(outfile.prefix,"GO.BP.Significant.enrichment.plot.pdf",sep=".")
      pdf(BP.plot.out)
      print(sig.BP.plot)
      dev.off()
      
      #### plot SVG of Directed Acyclic graph of 15 most signiifcnat GO terms.
      GO.BP.top.DAG <- paste(outfile.prefix,"GO.BP.top.DAG",sep='.')
      svgRes <- getAmigoTree(top.result.BP$GOBPID, color="red", pvalues =top.result.BP$Pvalue, filename=GO.BP.top.DAG, picType="svg", saveResult=TRUE)
      }
    
    # Add gene names to results table
    allgos.BP <- geneIdUniverse(hgOver.BP)
    output.BP.match <- NULL
    for (i in 1:nrow(result.BP))
    {
      
      go.holding = result.BP$GOBPID[i]
      all.entrez.in.GO <- as.vector(unlist(allgos.BP[go.holding]))
      goi.entrez.in.GO <- intersect(all.entrez.in.GO,goi.entrez)
      input.in.GO.IDs <- all.genes[all.genes$Entrez %in% goi.entrez.in.GO, 1]
      input.in.GO.IDs <- unique(input.in.GO.IDs[input.in.GO.IDs != ""])
      input.in.GO.IDs <- paste(input.in.GO.IDs, collapse = " ")
      input.in.GO.external <- unique(all.genes[all.genes$Entrez %in% goi.entrez.in.GO, 3])
      input.in.GO.external <- paste(input.in.GO.external, collapse = " ")
      temp <- cbind(go.holding,input.in.GO.IDs,input.in.GO.external)
      output.BP.match <- rbind(output.BP.match,temp)
    }
    result.BP <- merge(result.BP, output.BP.match, by.x ="GOBPID", by.y="go.holding", all.x=TRUE)
    BP.table.out = paste(outfile.prefix,"GO.BP.table",sep=".")
    result.BP <- result.BP[order(result.BP$Pvalue),] # order by Pvalue
    write.table(result.BP, file=BP.table.out, row.names = FALSE, col.names=TRUE,sep = '\t', quote=FALSE)
    
  }
  
  
  
  # Molecular Function
  params.MF <- new('GOHyperGParams',
                   geneIds=goi.entrez,
                   universeGeneIds=universe,
                   ontology='MF',
                   pvalueCutoff=GO.cutoff,
                   conditional=F,
                   testDirection='over',
                   annotation=species.db
  )
  hgOver.MF <- hyperGTest(params.MF)
  result.MF <- summary(hgOver.MF)
  
  if (nrow(result.MF)==0 )
  {
    fail.GO.MF = 1
    cat(c("GO Molecular Function search identified no enriched terms","Probably too few IDs"),
        file=run.report, append=TRUE, sep='\n')
  }
  if (fail.GO.MF !=1)
  {
    result.MF <- result.MF[result.MF$Count >= min.genes.cutoff,] # filter those with < cut off count
    result.MF <- result.MF[order(-result.MF$OddsRatio),] # order by odds ratio
    
    top.result.MF <- head(result.MF,15)
    top.result.MF <- top.result.MF[order(top.result.MF$Pvalue),]
    top.result.MF$Term <- as.factor(top.result.MF$Term)
    top.result.MF$Term <- factor(top.result.MF$Term, levels = top.result.MF$Term)
    
    if (nrow(top.result.MF) > 0)
    {
      top.result.MF$Pvalue[top.result.MF$Pvalue == 0 ] <- 1e-10 # catches any where p value = 0 
      max.y.plot = 1.2*(max(-log10(top.result.MF$Pvalue)))
      sig.MF.plot <-
        ggplot(data = top.result.MF,
               aes(x = as.factor(Term), y = -log10(top.result.MF$Pvalue),
                   size = Count))+
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
      
      
      MF.plot.out = paste(outfile.prefix,"GO.MF.Significant.enrichment.plot.pdf",sep=".")
      pdf(MF.plot.out)
      print(sig.MF.plot)
      dev.off()
      
      #### plot SVG of Directed Acyclic graph of 15 most signiifcnat GO terms. 
      GO.MF.top.DAG <- paste(outfile.prefix,"GO.MF.top.DAG",sep='.')
      svgRes <- getAmigoTree(top.result.MF$GOMFID, color="red", pvalues =top.result.MF$Pvalue, filename=GO.MF.top.DAG, picType="svg", saveResult=TRUE)
      
    }
    
    # Add gene names to results table
    allgos.MF <- geneIdUniverse(hgOver.MF)
    output.MF.match <- NULL
    for (i in 1:nrow(result.MF))
    {
      go.holding = result.MF$GOMFID[i]
      all.entrez.in.GO <- as.vector(unlist(allgos.MF[go.holding]))
      goi.entrez.in.GO <- intersect(all.entrez.in.GO,goi.entrez)
      input.in.GO.IDs <- all.genes[all.genes$Entrez %in% goi.entrez.in.GO, 1]
      input.in.GO.IDs <- unique(input.in.GO.IDs[input.in.GO.IDs != ""])
      input.in.GO.IDs <- paste(input.in.GO.IDs, collapse = " ")
      input.in.GO.external <- unique(all.genes[all.genes$Entrez %in% goi.entrez.in.GO, 3])
      input.in.GO.external <- paste(input.in.GO.external, collapse = " ")
      temp <- cbind(go.holding,input.in.GO.IDs,input.in.GO.external)
      output.MF.match <- rbind(output.MF.match,temp)
    }
    result.MF <- merge(result.MF, output.MF.match, by.x ="GOMFID", by.y="go.holding", all.x=TRUE)
    MF.table.out = paste(outfile.prefix,"GO.MF.table",sep=".")
    result.MF <- result.MF[order(result.MF$Pvalue),] # order by Pvalue
    write.table(result.MF, file=MF.table.out, row.names = FALSE, col.names=TRUE,sep = '\t', quote=FALSE)
  }
  
  # Cellular Compartment
  params.CC <- new('GOHyperGParams',
                   geneIds=goi.entrez,
                   universeGeneIds=universe,
                   ontology='CC',
                   pvalueCutoff=GO.cutoff,
                   conditional=F,
                   testDirection='over',
                   annotation=species.db
  )
  hgOver.CC <- hyperGTest(params.CC)
  result.CC <- summary(hgOver.CC)
  
  if (nrow(result.CC)==0 )
  {
    fail.GO.CC = 1
    cat(c("GO Cellular location search identified no enriched terms","Probably too few IDs"),
        file=run.report, append=TRUE, sep='\n')
  }
  
  if (fail.GO.CC !=1)
  {
    result.CC <- result.CC[result.CC$Count >= min.genes.cutoff,] # filter those with < cut off count
    result.CC <- result.CC[order(-result.CC$OddsRatio),] # order by odds ratio
    top.result.CC <- head(result.CC,15)
    top.result.CC <- top.result.CC[order(top.result.CC$Pvalue),]
    top.result.CC$Term <- as.factor(top.result.CC$Term)
    top.result.CC$Term <- factor(top.result.CC$Term, levels = top.result.CC$Term)
    
    if (nrow(top.result.CC) > 0)
    {
      top.result.CC$Pvalue[top.result.CC$Pvalue == 0 ] <- 1e-10 # catches any where p value = 0
      max.y.plot = 1.2*(max(-log10(top.result.CC$Pvalue)))
      sig.CC.plot <-
        ggplot(data = top.result.CC,
               aes(x = as.factor(Term), y = -log10(top.result.CC$Pvalue),
                   size = Count))+
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
      
      CC.plot.out = paste(outfile.prefix,"GO.CC.Significant.enrichment.plot.pdf",sep=".")
      pdf(CC.plot.out)
      print(sig.CC.plot)
      dev.off()
      
      #### plot SVG of Directed Acyclic graph of 15 most signiifcnat GO terms. 
      GO.CC.top.DAG <- paste(outfile.prefix,"GO.CC.top.DAG",sep='.')
      svgRes <- getAmigoTree(top.result.CC$GOCCID, color="red", pvalues =top.result.CC$Pvalue, filename=GO.CC.top.DAG, picType="svg", saveResult=TRUE)
      
    }
    
    
    # Add gene names to results table
    allgos.CC <- geneIdUniverse(hgOver.CC)
    output.CC.match <- NULL
    for (i in 1:nrow(result.CC))
    {
      go.holding = result.CC$GOCCID[i]
      all.entrez.in.GO <- as.vector(unlist(allgos.CC[go.holding]))
      goi.entrez.in.GO <- intersect(all.entrez.in.GO,goi.entrez)
      input.in.GO.IDs <- all.genes[all.genes$Entrez %in% goi.entrez.in.GO, 1]
      input.in.GO.IDs <- unique(input.in.GO.IDs[input.in.GO.IDs != ""])
      input.in.GO.IDs <- paste(input.in.GO.IDs, collapse = " ")
      input.in.GO.external <- unique(all.genes[all.genes$Entrez %in% goi.entrez.in.GO, 3])
      input.in.GO.external <- paste(input.in.GO.external, collapse = " ")
      temp <- cbind(go.holding,input.in.GO.IDs,input.in.GO.external)
      output.CC.match <- rbind(output.CC.match,temp)
    }
    result.CC <- merge(result.CC, output.CC.match, by.x ="GOCCID", by.y="go.holding", all.x=TRUE)
    CC.table.out = paste(outfile.prefix,"GO.CC.table",sep=".")
    result.CC <- result.CC[order(result.CC$Pvalue),] # order by Pvalue
    write.table(result.CC, file=CC.table.out, row.names = FALSE, col.names=TRUE,sep = '\t', quote=FALSE)
  }
}


##############################################################################
##############################################################################
#
# part 2 Pathway analysis
# 
##############################################################################
##############################################################################

if (doReactome == "yes")
{
  reactome.out <- enrichPathway(gene=goi.entrez,
                                #pvalueCutoff=0.05,
                                readable=T,
                                organism = reactome.spp,
                                pAdjustMethod = "BH",
                                qvalueCutoff = 0.01,
                                universe = universe
  )
  
  reactome.writeout <- (as.data.frame(reactome.out))
  
  if (nrow(reactome.writeout) == 0)
  {
    fail.reactome = 1
    cat(c("Reactome analysis identified no enriched pathways","Probably too few IDs"),
        file=run.report, append=TRUE, sep='\n')
  }
  
  if (fail.reactome==0)
  {
    reactome.table.out = paste(outfile.prefix,"reactome.pathway.enrichment.table",sep=".")
    write.table(reactome.writeout, file=reactome.table.out, row.names = FALSE, col.names = TRUE, quote=FALSE, sep='\t')
    
    reactome.dot <- dotplot <- dotplot(
      reactome.out,
      showCategory=15,
      font.size = 12
    )
    reactome.plot.out = paste(outfile.prefix,"reactome.pathway.enrichment.dotplot.tiff",sep=".")
    tiff(filename=reactome.plot.out,
         width = 320,
         height = 240,
         units = "mm",
         res=800
    )
    print(reactome.dot)
    dev.off()
    
    reactome.map.out = paste(outfile.prefix,"reactome.pathway.enrichment.enrichmap.tiff",sep=".")
    tiff(filename=reactome.map.out,
         width = 320,
         height = 240,
         units = "mm",
         res=800,
         type = "Xlib",
         pointsize = 12
    )
    enrichMap(reactome.out,
              layout=igraph::layout.kamada.kawai,
              vertex.label.cex = 0.8
    )
    dev.off()
  }
}

####
##############################################################################
##############################################################################
#
# part 3 KEGG analysis
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
    cat(c("KEGG analysis identified no enriched pathways","Probably too few IDs"),
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
    
    
    universe.size = as.numeric(length(universe))
    total.goi.size = as.numeric(length(goi.entrez))
    

# run phyper for each pathway in list and generate table of pathways passing cut off after FDR qvalue calculation
    working.pathways <- unique(matching.kegg.sets.spp.df$kegg.id)
    
    pathways.hypergeometric.results <- data.frame("Pathway"= character(0),"p.val"= numeric(0),"FDR q.val"= numeric(0),"ID"= character(0), "entrez.ids"= numeric(0), "external.ids"= character(0))
    pathways.hypergeometric.results.sig <- data.frame("Pathway"= character(0),"p.val"= numeric(0),"FDR q.val"= numeric(0), "goi.count"= numeric(0))
    
    
    for (i in 1:length(working.pathways)){
      current.pathway = working.pathways[i]
      goi.in.pathway <- as.numeric(nrow(goi.matching.kegg.sets.spp.df[goi.matching.kegg.sets.spp.df$kegg.id == current.pathway, ]))
      total.genes.in.pathway <- as.numeric(nrow(matching.kegg.sets.spp.df[matching.kegg.sets.spp.df$kegg.id == current.pathway, ]))
      
      pval <- phyper(goi.in.pathway,total.genes.in.pathway,(universe.size-total.genes.in.pathway),total.goi.size, lower.tail=FALSE)
      qval <- p.adjust(pval, method = "fdr", n = nrow(keggres.pathways.out))
      
      current.goi <- goi.matching.kegg.sets.spp.df[goi.matching.kegg.sets.spp.df$kegg.id == current.pathway, ]
      current.goi <- current.goi[1]
      current.goi.entrez.ids <- as.numeric(as.character(current.goi$entrez.id))
      
      current.goi.ens <- all.genes.entrez[all.genes.entrez$Entrez %in% current.goi.entrez.ids ,]
      
      current.goi.ens.ids <- unique(current.goi.ens$ID)
      current.goi.ext.ids <- unique(current.goi.ens$Name)
      
      
      current.goi.ens.ids <- paste(current.goi.ens.ids, collapse=", ")
      current.goi.entrez.ids <- paste(current.goi.entrez.ids, collapse=", ")
      current.goi.ext.ids <- paste(current.goi.ext.ids, collapse=", ")
      
      current.out <- as.data.frame(cbind(current.pathway,pval,qval,current.goi.ens.ids,current.goi.entrez.ids,current.goi.ext.ids))
      current.sig.out <- as.data.frame(cbind(current.pathway,pval,qval,goi.in.pathway))
      
      pathways.hypergeometric.results <- rbind(pathways.hypergeometric.results, current.out)
      pathways.hypergeometric.results.sig <- rbind(pathways.hypergeometric.results.sig, current.sig.out)
    }  

    
    colnames(pathways.hypergeometric.results) <- c("Pathway","p.val","FDR q.val","Ensembl.ids","Entrez.ids","External.ids")
    colnames(pathways.hypergeometric.results.sig) <- c("Pathway","p.val","FDR q.val","goi.count")
    
    # make FDR q.val numeric and sort 
    pathways.hypergeometric.results$`FDR q.val` <- as.numeric(as.character(pathways.hypergeometric.results$`FDR q.val`))
    pathways.hypergeometric.results <-  pathways.hypergeometric.results[with(pathways.hypergeometric.results, order(pathways.hypergeometric.results$`FDR q.val`)), ]
    
    kegg.table.out = paste(outfile.prefix,"kegg.pathway.enrichment.table",sep=".")
    write.table(pathways.hypergeometric.results,file=kegg.table.out, row.names = FALSE, col.names = TRUE, quote = FALSE, sep ='\t')
    
    pathways.hypergeometric.results.sig$`FDR q.val` <- as.numeric(as.character(pathways.hypergeometric.results.sig$`FDR q.val`))
    pathways.hypergeometric.results.sig$goi.count <- as.numeric(as.character(pathways.hypergeometric.results.sig$goi.count))
    pathways.hypergeometric.results.sig <- pathways.hypergeometric.results.sig[pathways.hypergeometric.results.sig$`FDR q.val` < kegg.qval.cutoff & pathways.hypergeometric.results.sig$goi.count >= min.genes.cutoff, ]
    
    ##############################################################################################  
    # draw plot of enriched pathways
    ############################################################################################## 
   
    if (nrow(pathways.hypergeometric.results.sig)>0)
    {
      
      
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
      
  


    if (stats.KEGG.fail == 0)
    {
      cat(c("KEGG analysis no terms pass statistical cutoff"),
          file=run.report, append=TRUE, sep='\n')
    }
    
  }
}

