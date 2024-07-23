
#----------------------------------------------------------------------------------------------------------------------------
######this all for 04_modelInterpretation script#####


####FROM HERE#####
topFeats <- readRDS("top5Feats.rds")

function_forFeats <- function(listofFeatures){
  weights <- c(listofFeatures[,1])
  names(weights) <- rownames(listofFeatures)
  #return(weights)
  weights <- sort(weights, decreasing = T) 
  gseGO <- gseGO(geneList = weights,
                 ont = "BP",
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENSEMBL",
                 eps=0,
                 nPermSimple = 100000
  )
  # return(gseGO)
  # return(weights)
  feat <- colnames(listofFeatures[1])
  res <- gseGO@result
  res <- res %>% arrange(desc(enrichmentScore))
  p <- ggplot(res[1:15,], aes(x = enrichmentScore, y = forcats::fct_reorder(Description, enrichmentScore)))+
    geom_point(size=3.5)+
    theme_bw(base_size = 18)+
    scale_y_discrete(labels = label_wrap_gen(40))+
    xlab("Enrichment score")+
    ylab("Description")+
    facet_wrap(~ paste0(colnames(listofFeatures[1])))
  return(paste0(feat) = p)
}

tempRes <- lapply(topFeats, function_forFeats)
temp2 <- function_forFeats(topFeats[[5]])
gseGOplots <- lapply(topFeats, function_forFeats)
# now I wanna plot the encrihment score
firstREsrEs <- firstREsrEs %>% arrange(desc(enrichmentScore))
  
ggplot(firstREsrEs[1:15,], aes(x = enrichmentScore, 
                               y = forcats::fct_reorder(Description, enrichmentScore)))+
  geom_point(size=3.5)+
  theme_bw(base_size = 18)+
  scale_y_discrete(labels = label_wrap_gen(40))+
  xlab("Enrichment score")+
  ylab("Description")





##external validation part 3
humanPhenotype <- read.delim("genes_for_HP_0000924", col.names = c("geneID", "geneSymbol"))
externalDatabase <- read.delim("https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt", header = F, 
                               sep = "\t")
externalDatabase$V6 <- NULL
colnames(externalDatabase)<- c("humanMarkerSymbol", "humanEntrezGeneID", "mouseMarkerSymbol", "mgiMarkerID", "mammalianPhenotypeID")
head(externalDatabase)
#search for skeletal phenotype and 
patterns <- c("MP:0005390", "MP:0005371")

exDb <- externalDatabase %>% filter(grepl(paste(patterns, collapse = '|'), mammalianPhenotypeID))

#so now we want to know... how many of the model's predictions are known to be associated in the
# mammalian and human phenotype ontology
#but first might need to map the external genes to human ensembl IDs since the identifiers are different 
#to what I am working with now (ensembl IDs)
