## ----echo=FALSE, results="hide"-----------------------------------------------
knitr::opts_chunk$set(error=FALSE, warning=FALSE, message=FALSE,
                        fig.retina = 0.75, crop = NULL)
library(BiocStyle)

## ----library, echo=FALSE------------------------------------------------------
library(cytomapper)

## ----quickstart-load-data-----------------------------------------------------
data(pancreasSCE)
data(pancreasImages)
data(pancreasMasks)

## ----quickstart-plotPixels----------------------------------------------------
plotPixels(image = pancreasImages, colour_by = c("H3", "CD99", "CDH"))

## ----quickstart-plotCells-2---------------------------------------------------
plotCells(mask = pancreasMasks, object = pancreasSCE,
            cell_id = "CellNb", img_id = "ImageNb", colour_by = "CD99",
            outline_by = "CellType")
plotCells(mask = pancreasMasks, object = pancreasSCE,
            cell_id = "CellNb", img_id = "ImageNb", 
            colour_by = "CellType")

## ----image-cell-id-sce--------------------------------------------------------
head(colData(pancreasSCE)[,"ImageNb"])
head(colData(pancreasSCE)[,"CellNb"])

## ----image-cell-id-cil--------------------------------------------------------
mcols(pancreasImages)[,"ImageNb"]
mcols(pancreasMasks)[,"ImageNb"]

## ----pancreasImages-----------------------------------------------------------
pancreasImages
mcols(pancreasImages)
channelNames(pancreasImages)
imageData(pancreasImages[[1]])[1:15,1:5,1]

## ----pancreasMasks------------------------------------------------------------
pancreasMasks
mcols(pancreasMasks)
imageData(pancreasMasks[[1]])[1:15,1:5]

## ----pancreasSCE--------------------------------------------------------------
pancreasSCE
names(colData(pancreasSCE))

## ----read-in-images-----------------------------------------------------------
# Read in masks
path.to.images <- system.file("extdata", package = "cytomapper")
all_masks <- loadImages(path.to.images, pattern = "_mask.tiff")
all_masks

# Read in images
all_stacks <- loadImages(path.to.images, pattern = "_imc.tiff")
all_stacks

## ----add-metadata-------------------------------------------------------------
unique(pancreasSCE$ImageNb)
mcols(all_masks)$ImageNb <- c("1", "2", "3")
mcols(all_stacks)$ImageNb <- c("1", "2", "3")

## ----image-encoding-----------------------------------------------------------
head(unique(as.numeric(all_masks[[1]])))

## ----scaleImage-1-------------------------------------------------------------
all_masks <- scaleImages(all_masks, 2^16-1)
head(unique(as.numeric(all_masks[[1]])))

## ----as.is--------------------------------------------------------------------
all_masks_2 <- loadImages(path.to.images, pattern = "_mask.tiff", as.is = TRUE)
head(unique(as.numeric(all_masks_2[[1]])))

## ----channelNames-example-----------------------------------------------------
channelNames(all_stacks) <- c("H3", "CD99", "PIN", "CD8a", "CDH")

## ----measureObjects-----------------------------------------------------------
sce <- measureObjects(all_masks, all_stacks, img_id = "ImageNb")
sce

## ----counts-sce---------------------------------------------------------------
counts(sce)[1:5, 1:5]
colData(sce)

## ----mcols--------------------------------------------------------------------
mcols(pancreasImages)
mcols(pancreasImages)$PatientID <- c("Patient1", "Patient2", "Patient3")
mcols(pancreasImages)

## ----subsetting-SimpleList----------------------------------------------------
pancreasImages[1]
pancreasImages[[1]]

## ----get-set-images-----------------------------------------------------------
cur_image <- getImages(pancreasImages, "E34_imc")
cur_image
setImages(pancreasImages, "New_image") <- cur_image
pancreasImages
mcols(pancreasImages)

## ----set-images---------------------------------------------------------------
names(cur_image) <- "Replacement"
setImages(pancreasImages, 2) <- cur_image
pancreasImages
mcols(pancreasImages)

## ----set-images-2-------------------------------------------------------------
setImages(pancreasImages, "J02_imc") <- cur_image
pancreasImages
mcols(pancreasImages)

## ----deleting-images----------------------------------------------------------
setImages(pancreasImages, c("Replacement", "New_image")) <- NULL
pancreasImages

## ----get-set-channels---------------------------------------------------------
cur_channel <- getChannels(pancreasImages, "H3")
cur_channel
channelNames(cur_channel) <- "New_H3"
setChannels(pancreasImages, 1) <- cur_channel
pancreasImages

## ----channelnames-------------------------------------------------------------
channelNames(pancreasImages)
channelNames(pancreasImages) <- c("ch1", "ch2", "ch3", "ch4", "ch5")
pancreasImages

## ----mergechannels------------------------------------------------------------
cur_channels <- getChannels(pancreasImages, 1:2)
channelNames(cur_channels) <- c("new_ch1", "new_ch2")
pancreasImages <- mergeChannels(pancreasImages, cur_channels)
pancreasImages

## ----endoapply-example--------------------------------------------------------
data("pancreasImages")

# Performing a gaussian blur
pancreasImages <- endoapply(pancreasImages, gblur, sigma = 1)

## ----normalize-default--------------------------------------------------------
data("pancreasImages")

# Default normalization
cur_images <- normalize(pancreasImages)

## ----normalize-image-wise-----------------------------------------------------
# Image-wise normalization
cur_images <- normalize(pancreasImages, separateImages = TRUE)

## ----normalize-clippingrange--------------------------------------------------
# Percentage-based clipping range
cur_images <- normalize(pancreasImages)
cur_images <- normalize(cur_images, inputRange = c(0, 0.9))
plotPixels(cur_images, colour_by = c("H3", "CD99", "CDH"))

## ----normalize-clippingrange-channelwise--------------------------------------
# Channel-wise clipping
cur_images <- normalize(pancreasImages,
                        inputRange = list(H3 = c(0, 70), CD99 = c(0, 100)))

## ----bcg-pixels---------------------------------------------------------------
data("pancreasImages")
# Increase contrast for the CD99 and CDH channel
plotPixels(pancreasImages,
            colour_by = c("H3", "CD99", "CDH"),
            bcg = list(CD99 = c(0,2,1),
                        CDH = c(0,2,1)))

## ----outline-1-pixels---------------------------------------------------------
plotPixels(pancreasImages, mask = pancreasMasks,
            object = pancreasSCE, img_id = "ImageNb",
            cell_id = "CellNb",
            colour_by = c("H3", "CD99", "CDH"),
            outline_by = "CellType")

## ----subsetting-pixels--------------------------------------------------------
cur_images <- getImages(pancreasImages, "J02_imc")
plotPixels(cur_images, colour_by = c("H3", "CD99", "CDH"))

## ----exprs_values-cells-------------------------------------------------------
plotCells(pancreasMasks, object = pancreasSCE,
            img_id = "ImageNb", cell_id = "CellNb",
            colour_by = c("CD8a", "PIN"),
            exprs_values = "exprs")

## ----colour-cells-------------------------------------------------------------
plotCells(pancreasMasks, object = pancreasSCE,
            img_id = "ImageNb", cell_id = "CellNb",
            colour_by = c("CD99", "CDH"),
            outline_by = "CellType",
            colour = list(CD99 = c("black", "red"),
                        CDH = c("black", "white"),
                        CellType = c(celltype_A = "blue",
                                    celltype_B = "green",
                                    celltype_C = "yellow")))

## ----subsetting-SCE-----------------------------------------------------------
cur_sce <- pancreasSCE[,colData(pancreasSCE)$CellType == "celltype_A"]
plotCells(pancreasMasks, object = cur_sce,
            img_id = "ImageNb", cell_id = "CellNb",
            colour_by = "CellType",
            colour = list(CellType = c(celltype_A = "red")))

## ----customization, fig.cap = "Plot customization example."-------------------
plotCells(pancreasMasks, object = pancreasSCE,
            img_id = "ImageNb", cell_id = "CellNb",
            colour_by = "CD99",
            outline_by = "CellType",
            background_colour = "white",
            missing_colour = "black",
            scale_bar = list(length = 30,
                            label = expression("30 " ~ mu * "m"),
                            cex = 2,
                            lwidth = 10,
                            colour = "cyan",
                            position = "bottomleft",
                            margin = c(5,5),
                            frame = 3),
            image_title = list(text = c("image_1", "image_2", "image_3"),
                            position = "topleft",
                            colour = "cyan",
                            margin = c(2,10),
                            font = 3,
                            cex = 2),
            legend = list(colour_by.title.font = 2,
                            colour_by.title.cex = 1.2,
                            colour_by.labels.cex = 0.7,
                            outline_by.legend.cex = 0.3,
                            margin = 10),
            margin = 2,
            thick = TRUE)

## ----return_plot--------------------------------------------------------------
cur_out <- plotPixels(pancreasImages, colour_by = c("H3", "CD99", "CDH"),
                return_plot = TRUE, return_images = TRUE, 
                display = "single")

## ----return_plot-2------------------------------------------------------------
cur_out$plot$E34_imc

## ----return_images------------------------------------------------------------
plot(cur_out$images$G01_imc)

## ----plot-integration---------------------------------------------------------
library(cowplot)
library(ggplot2)

g1 <- ggplot(mtcars) + geom_point(aes(cyl, hp))
g2 <- plotCells(pancreasMasks, object = pancreasSCE,
            img_id = "ImageNb", cell_id = "CellNb",
            colour_by = "CellType", return_plot = TRUE)
g2 <- ggdraw(g2$plot, clip = "on")

plot_grid(g1, g2)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

