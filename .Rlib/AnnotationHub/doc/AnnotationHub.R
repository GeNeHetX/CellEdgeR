## ----style, echo = FALSE, results = 'asis'--------------------------------------------------------
BiocStyle::markdown()

## ----library, message=FALSE-----------------------------------------------------------------------
library(AnnotationHub)

## ----AnnotationHub--------------------------------------------------------------------------------
ah = AnnotationHub()

## ----show-----------------------------------------------------------------------------------------
ah

## ----dataprovider---------------------------------------------------------------------------------
unique(ah$dataprovider)

## ----species--------------------------------------------------------------------------------------
head(unique(ah$species))

## ----rdataclass-----------------------------------------------------------------------------------
head(unique(ah$rdataclass))

## ----dm1------------------------------------------------------------------------------------------
dm <- query(ah, c("ChainFile", "UCSC", "Drosophila melanogaster"))
dm

## ----show2----------------------------------------------------------------------------------------
df <- mcols(dm)

## ----length---------------------------------------------------------------------------------------
length(ah)

## ----subset---------------------------------------------------------------------------------------
ahs <- query(ah, c('inparanoid8', 'ailuropoda'))
ahs

## ----BiocHubsShiny, eval=FALSE--------------------------------------------------------------------
# BiocHubsShiny::BiocHubsShiny()

## ----echo=FALSE,fig.cap="BiocHubsShiny query with terms dataprovider = 'inparanoid' and species = 'ailuropoda'"----
knitr::include_graphics("ailMel1_BiocHubsShiny.png")

## ----dm2------------------------------------------------------------------------------------------
dm
dm["AH15146"]

## ----dm3------------------------------------------------------------------------------------------
dm[["AH15146"]]

## ----show-2---------------------------------------------------------------------------------------
ah

## ----snapshot-------------------------------------------------------------------------------------
snapshotDate(ah)

## ----possibleDates--------------------------------------------------------------------------------
pd <- possibleDates(ah)
pd

## ----setdate, eval=FALSE--------------------------------------------------------------------------
# snapshotDate(ah) <- pd[1]

## ----clusterOptions1, eval=FALSE------------------------------------------------------------------
# library(AnnotationHub)
# hub <- AnnotationHub()
# gr <- hub[["AH50773"]]  ## downloaded once
# txdb <- makeTxDbFromGRanges(gr)  ## build on the fly

## ----clusterOptions2, eval=FALSE------------------------------------------------------------------
# library(AnnotationDbi)  ## if not already loaded
# txdb <- loadDb("/locationToFile/mytxdb.sqlite")

## ----eval=FALSE-----------------------------------------------------------------------------------
# proxy <- "http://my_user:my_password@myproxy:8080"
# AnnotationHub::setAnnotationHubOption("PROXY", proxy)
# ## or
# ExperimentHub::setExperimentHubOption("PROXY", proxy)

## ----eval=FALSE-----------------------------------------------------------------------------------
# ssl_opts <- list(verbose = 1L,ssl_verifypeer = 0L, ssl_verifyhost = 0L)
# 

## ----sessionInfo----------------------------------------------------------------------------------
sessionInfo()

