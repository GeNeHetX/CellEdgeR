## ----style, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(error=FALSE, warning=FALSE, message=FALSE,
    fig.retina = 0.75)
library(BiocStyle)

## ----library, echo=FALSE, results='hide'--------------------------------------
suppressPackageStartupMessages(c(
    library(SingleCellExperiment),
    library(cytomapper),
    library(imcdatasets)
))

## ----list-datasets------------------------------------------------------------
datasets <- listDatasets()
datasets <- as.data.frame(datasets)
datasets$FunctionCall <- sprintf("`%s`", datasets$FunctionCall)
knitr::kable(datasets)

## ----import-dataset-----------------------------------------------------------
sce <- IMMUcan_2022_CancerExample("sce")
sce

## ----import-images------------------------------------------------------------
images <- IMMUcan_2022_CancerExample("images")
images

## ----import-masks-------------------------------------------------------------
masks <- IMMUcan_2022_CancerExample("masks")
masks

## ----on_disk------------------------------------------------------------------
# Create temporary location
cur_path <- tempdir()

masks <- IMMUcan_2022_CancerExample(data_type = "masks", on_disk = TRUE,
    h5FilesPath = cur_path)
masks

## ----function-help------------------------------------------------------------
?IMMUcan_2022_CancerExample

## ----access-metadata, eval = FALSE--------------------------------------------
# IMMUcan_2022_CancerExample(data_type = "sce", metadata = TRUE)
# IMMUcan_2022_CancerExample(data_type = "images", metadata = TRUE)
# IMMUcan_2022_CancerExample(data_type = "masks", metadata = TRUE)

## ----usage-subset-------------------------------------------------------------
cur_images <- images[1:5]
cur_masks <- masks[1:5]

## ----usage-pixel--------------------------------------------------------------
plotPixels(
    cur_images,
    colour_by = c("CD8a", "CD68", "CDH1"),
    bcg = list(
        CD8a = c(0,4,1),
        CD68 = c(0,5,1),
        CDH1 = c(0,5,1)
    )
)

## ----usage-cell---------------------------------------------------------------
plotCells(
    cur_masks, object = sce,
    img_id = "image_number", cell_id = "cell_number",
    colour_by = c("CD8a", "CD68", "CDH1"),
    exprs_values = "exprs"
)

## ----usage-outline------------------------------------------------------------
plotPixels(
    cur_images, mask = cur_masks, object = sce,
    img_id = "image_number", cell_id = "cell_number",
    outline_by = "cell_type",
    colour_by = c("CD8a", "CD68", "CDH1"),
    bcg = list(
        CD8a  = c(0,5,1),
        CD68 = c(0,5,1),
        CDH1 = c(0,5,1)
    )
)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

