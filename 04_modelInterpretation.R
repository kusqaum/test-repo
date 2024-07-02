nmfMatrices <- readRDS("nmf_wMatrices.rds")

nmfMatrices2 <- readRDS("nmfDataframes.rds")

makeSmallFctn <- function(matrices){
  matrices <- matrices[,(ncol(matrices)-2):ncol(matrices)]
  return(matrices)
}
makeSmallFctn2 <- function(matrices){
  matrices <- matrices[,(ncol(matrices)-4):ncol(matrices)]
  return(matrices)
}

res <- lapply(nmfMatrices2, FUN = makeSmallFctn)
res2 <- lapply(nmfMatrices2, FUN = makeSmallFctn2)

allNmfMatrix <- list()
allNmfMatrix[[1]] <- res[[1]]
allNmfMatrix[[2]] <- res2[[2]]


library(tidymodels)
library(themis)
library(skimr)
library(tidyverse)
library(randomForest)
library(vip)
library(future)

#function for preprocessing
split_processingData_fctn <- function(data, proportion){
  set.seed(123)
  dataSplit <- initial_split(data, prop = proportion)
  trainData <- training(dataSplit)
  testData <- testing(dataSplit)
  #create recipe
  rec <-recipe(significant~., data = trainData) %>%
    step_downsample(significant, under_ratio = 1, seed = 345)
  #return(list(rec, trainData, testData))
  return(list("train" = as.data.frame(trainData), "test"=as.data.frame(testData), "recipe"=rec))
}
test <- map(.x = allNmfMatrix, .f = split_processingData_fctn, proportion=0.8)
nmfProcessed <- lapply(X=allNmfMatrix, FUN = split_processingData_fctn, proportion =0.8)
temp_rec <- test[[1]]$recipe
temp_rec %>% 
  prep() %>% 
  juice() %>%
  ggplot(aes(factor(significant))) +
  geom_bar(aes(y = (after_stat(count))/sum(after_stat(count))), colour="black",fill="lightgrey") +
  scale_y_continuous(labels = percent) +
  xlab("") + ylab("% of genes") +
  theme_minimal(base_size = 20)

# counts <- training3 %>% group_by(significant) %>%
#   summarise(count=n())
# ggplot(data = counts, aes(x=significant, y=count))+
#   geom_bar(stat='identity')

#
randForest_fctn <- function(preprocessResult, algorithm){
  if(algorithm == "NMF"){
    model_RF <- rand_forest(trees = tune(), mtry = tune(), min_n = tune(), 
                            mode = "classification") %>% set_engine("randomForest", importance = TRUE)
    set.seed(234)
    
    folds <- vfold_cv(data = preprocessResult$train, v=3)
    limit <- (ncol(preprocessResult$train))-1
    tuningGrid <- grid_regular(
      trees(range = c(1,500)),
      mtry(range = c(1,limit)),
      min_n(range = c(1,limit)),
      levels = limit
    )
    #build workflow
    wkflow <- workflow() %>%
      add_recipe(preprocessResult$recipe) %>%
      add_model(model_RF)
    #tune model
    #plan(multisession, workers = 16)
    res <- tune_grid(
      wkflow,
      resamples = folds,
      grid = tuningGrid,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(roc_auc),
    )
    paramsPlot <- autoplot(res)
    paramsPlot2 <- res %>%
      collect_metrics() %>%
      ggplot(aes(x=mtry, y=mean, col=as.factor(min_n)))+
      geom_point()+
      geom_line()
    
    final_model <- res %>% select_best(metric = "roc_auc")
    
    final_fit <- finalize_workflow(wkflow, final_model) %>%
      fit(data = preprocessResult$train) 
    
    importancePlot <- final_fit %>% extract_fit_parsnip()%>%
      vip(geom='point',aes = list(colour="black", fill='lightblue', alpha=0.7))+
      theme_classic()
    
    importanceDf <- final_fit %>% extract_fit_parsnip() %>%
      vi() %>% as.data.frame()
    
    aug <- augment(final_fit, preprocessResult$test)
    aug_m <- aug %>% mutate(dim = ncol(preprocessResult$train), Algorithm = algorithm)
    
    roc_auc <- roc_auc(aug, significant, .pred_TRUE)
    two_classCurve <- roc_curve(aug, truth = significant,
                                .pred_TRUE)
    rocCurve <- autoplot(two_classCurve)
    
    #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
    result <- list("tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
    saveRDS(result, sprintf("%sMLRes_%s.rds",algorithm, ncol(preprocessResult$train)))
    
    return(list("tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m))
    
  }
  else if(algorithm == "PCA"){
    model_RF <- rand_forest(trees = tune(), mtry = tune(), min_n = tune(), 
                            mode = "classification") %>% set_engine("randomForest", importance = TRUE)
    set.seed(234)
    
    folds <- vfold_cv(data = preprocessResult$train, v=3)
    limit <- (ncol(preprocessResult$train))-1
    tuningGrid <- grid_regular(
      trees(range = c(1,500)),
      mtry(range = c(1,limit)),
      min_n(range = c(1,limit)),
      levels = limit
    )
    #build workflow
    wkflow <- workflow() %>%
      add_recipe(preprocessResult$recipe) %>%
      add_model(model_RF)
    #tune model
    #plan(multisession, workers = 16)
    res <- tune_grid(
      wkflow,
      resamples = folds,
      grid = tuningGrid,
      control = control_grid(save_pred = TRUE),
      metrics = metric_set(roc_auc),
    )
    paramsPlot <- autoplot(res)
    paramsPlot2 <- res %>%
      collect_metrics() %>%
      ggplot(aes(x=mtry, y=mean, col=as.factor(min_n)))+
      geom_point()+
      geom_line()
    
    final_model <- res %>% select_best(metric = "roc_auc")
    
    final_fit <- finalize_workflow(wkflow, final_model) %>%
      fit(data = preprocessResult$train) 
    
    importancePlot <- final_fit %>% extract_fit_parsnip()%>%
      vip(geom='point',aes = list(colour="black", fill='lightblue', alpha=0.7))+
      theme_classic()
    
    importanceDf <- final_fit %>% extract_fit_parsnip() %>%
      vi() %>% as.data.frame()
    
    aug <- augment(final_fit, preprocessResult$test)
    aug_m <- aug %>% mutate(dim = ncol(preprocessResult$train), Algorithm = algorithm)
    
    roc_auc <- roc_auc(aug, significant, .pred_TRUE)
    two_classCurve <- roc_curve(aug, truth = significant,
                                .pred_TRUE)
    rocCurve <- autoplot(two_classCurve)
    
    #dims <- c("5pca","10pca","50pca","100pca","150pca","200pca")
    result <- list("tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                   "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m)
    saveRDS(result, sprintf("%sMLRes_%s.rds",algorithm, ncol(preprocessResult$train)))
    
    return(list("tuningPlots" =paramsPlot, "importancePlot"=importancePlot, "importanceDf"=importanceDf,
                "finalFit" = final_fit, "AUC"= roc_auc, "roc_curve" = rocCurve, "aug"= aug_m))
  }
  
}
availableCores()
plan(multisession, workers=availableCores())
# nmfDataMlResList <- lapply(nmfProcessed, FUN=randForest_fctn)
tempRes <- map(.x =test , .f = randForest_fctn)
ttessterr <- map(.x= test, .f = randForest_fctn_temp, algorithm = "NMF")

#----------------------------------------------------------------------------------------------------------------------------
######this all for 04_modelInterpretation script#####

library(tidyverse)
library(clusterProfiler)
# need to read in results from last script
# nmfRes <- readRDS("processed/tempMLResults.rds")
# read in also the list of different dimensions:
# nmfDataframes <- readRDS("processed/nmfDataframes.rds")

#### get all the auc scores ####
auc <-lapply(nmfDataMlResList, function(x){
  x$AUC$.estimate
})

df_auc <- lapply(auc, function(x){
  df<- data.frame(auc_scores = x)
})
View(df_auc)
# here just binding all the scores for each model
aucDf <- do.call("rbind", df_auc)
head(aucDf)
dimension <- data.frame(dimensions = c(2,4), Algorithm = "NMF")
resDf <- cbind(aucDf, dimension)
head(resDf)
pcaAucDf <- data.frame(auc_scores = c(0.5, 0.509), dimensions = c(2,4), Algorithm = 'PCA')
bothAUC <- rbind(resDf, pcaAucDf)

#then plotting
ggplot(bothAUC, aes(x=as.factor(dimensions), y=auc_scores, col=Algorithm))+
  geom_point(size=4.5) +
  theme(panel.border = element_rect(colour = "black", linewidth = 0.35, fill=NA),
        panel.background = element_blank(),
        #strip.text = element_text(),
        #legend.position = "none",
        #legend.text = element_text("dimensionality"),
        #legend.title = element_text("dimension"),
        aspect.ratio = 1)+
  scale_colour_manual(values=c("tan2", "tomato"))+
  xlab(expression(italic(" k") * " dimensions"))+
  ylab("AUC") +
  guides(colour = guide_legend(title = "Algorithm"))+
  theme_bw(base_size = 16)+
  # facet_wrap(~Algorithm)

#finding which dimension the model that gave the highest auc score was trained on
posBestFit <- which.max(resDf$auc)
bestDimension <- resDf[which.max(resDf$auc),]$dimensions
for (d in 1:length(nmfDataMlResList)){
  #print(nmfRes[[d]])
  if(ncol(allNmfMatrix[[d]])-1 == bestDimension){
    # print(allNmfMatrix[[d]])
    bestDf <- as.data.frame(allNmfMatrix[[d]])
  }
}
#find the final fit for that model
bestFit <- nmfDataMlResList[[posBestFit]]$finalFit
#lets look at the variable importance
bestModVarImport <- nmfDataMlResList[[posBestFit]]$importanceDf
bestModVarImport <- bestModVarImport %>% 
  mutate(sign = case_when(Importance<0 ~"negative", TRUE~"positive"))

#subset for the ones that are positive sign
impFeats <- bestModVarImport %>% filter(sign=="positive")
#then extract the variables that contribute to model's predictions
modelFeats <- impFeats$Variable
#find those features 
modelFeatsDf <- bestDf %>% select(all_of(modelFeats))# these are gonna be input for GSEA 

unlabelledGenes$prediction <- predict(bestFit, unlabelledGenes)
ggplot(bestModVarImport, aes(x=Variable, y=Importance, col = sign))+
  geom_point(size=4)+
  theme_bw(base_size = 18)+
  theme(legend.position="none")+
  xlab("")

vip(bestModVarImport, geom = 'point', mapping=aes(colour = sign))

####clusterprofiler code####
# enrichKEGG(gene = row.names(bestDf))
geneList <- data.frame(geneID = row.names(modelFeatsDf), Weight = modelFeatsDf$V8)
geneList_2<- as.data.frame(t(geneList))
genes <- c(row.names(modelFeatsDf))
weights <- as.vector(modelFeatsDf$V8)
weights <- sort(weights,decreasing = T)
weights <- as.vector(weights)
names(weights) <- genes
class(weights)
#genelis <- as.vector(geneList)
class(genelis)
# enr
# enrichKEGG(gene = weights, 
#            organism = "hsa", keyType = )

gseGO <- gseGO(geneList = weights,
               ont = "BP",
               OrgDb = org.Hs.eg.db,
               keyType = "ENSEMBL",
               maxGSSize=2000,
               eps=0,
               #pAdjustMethod="BH"
               )
gseGO@result%>%
  ggplot(aes())
p$data %>% #filter(p.adjust<0.3)%>%
  ggplot(aes(x=GeneRatio,y=forcats::fct_reorder(Description, GeneRatio)))+#, colour = p.adjust, size =Count))+
  geom_segment(aes(xend=0, yend = Description))+
  geom_point(aes(colour=p.adjust, size = Count))+
  
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE))+
  scale_size_continuous(range=c(1, 10)) +
  #scale_color_gradientn(colours = c(pal))+
  # scale_colour_gradientn(colours=c("coral", "tomato", "firebrick2", "firebrick3","slateblue4","plum4", "rosybrown", "plum3" ))+
  theme_bw(base_size = 14)+
  ylab("")

gseKEGG <- gseKEGG(geneList = weights,
                   organism = 'hsa',
                   )
p<- dotplot(gseGO)
##external validation part 3
humanPhenotype <- read.delim("genes_for_HP_0000924", col.names = c("geneID", "geneSymbol"))
externalDatabase <- read.delim("https://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt", header = F, 
                               sep = "\t")
externalDatabase$V6 <- NULL
colnames(externalDatabase)<- c("humanMarkerSymbol", "humanEntrezGeneID", "mouseMarkerSymbol", "mgiMarkerID", "mammalianPhenotypeID")
externalDatabase
#search for skeletal phenotype and 
patterns <- c("MP:0005390", "MP:0005371")

exDb <- externalDatabase %>% filter(grepl(paste(patterns, collapse = '|'), mammalianPhenotypeID))

#so now we want to know... how many of the model's predictions are known to be associated in the
# mammalian and human phenotype ontology
#but first might need to map the external genes to human ensembl IDs since the identifiers are different 
#to what I am working with now (ensembl IDs)
