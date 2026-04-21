## ----echo=FALSE, results="hide"-----------------------------------------------
knitr::opts_chunk$set(error=FALSE, warning=FALSE, message=FALSE,
                        fig.retina = 0.75, crop = NULL)
library(BiocStyle)

## ----library, echo=FALSE------------------------------------------------------
library(cytomapper)

## ----message=FALSE------------------------------------------------------------
library(HDF5Array)

# Define output directory
cur_dir <- getHDF5DumpDir()

path.to.images <- system.file("extdata", package = "cytomapper")
image.list <- loadImages(path.to.images, pattern = "mask.tiff",
                             on_disk = TRUE, h5FilesPath = cur_dir)

# Show list
image.list

# Scale images
image.list <- scaleImages(image.list, value = 2^16 - 1)
image.list$E34_mask

## -----------------------------------------------------------------------------
data("pancreasImages")

pancreasImages_onDisk <- CytoImageList(pancreasImages,
                                        on_disk = TRUE, 
                                        h5FilesPath = cur_dir)

# Image object
pancreasImages$E34_imc

# HDF5Array object
pancreasImages_onDisk$E34_imc

# Seed of HDF5Array object
seed(pancreasImages_onDisk$E34_imc)

# Size in memory
format(object.size(pancreasImages), units = "auto")
format(object.size(pancreasImages_onDisk), units = "auto")

## -----------------------------------------------------------------------------
pancreasImages_inMemory <- CytoImageList(pancreasImages_onDisk,
                                        on_disk = FALSE)

# Compare the image data to the original representation
identical(as.list(pancreasImages_inMemory), as.list(pancreasImages))

## ----sanity-check, echo=FALSE-------------------------------------------------
stopifnot(identical(as.list(pancreasImages), as.list(pancreasImages_inMemory)))

## -----------------------------------------------------------------------------
# Size of object in memory
format(object.size(pancreasImages_onDisk), units = "auto")

# Size of object on disk in kB
file.info(paste0(cur_dir, "/E34_imc.h5"))[,"size"] / 1000

pancreasImages_norm <- normalize(pancreasImages_onDisk)

seed(pancreasImages_norm$E34_imc)

# Size of object in memory
format(object.size(pancreasImages_norm), units = "auto")

# Size of object on disk in kB
file.info(paste0(cur_dir, "/E34_imc.h5"))[,"size"] / 1000

## -----------------------------------------------------------------------------
cur_Images1 <- pancreasImages_onDisk
cur_Images2 <- getChannels(pancreasImages_onDisk, 2)
channelNames(cur_Images2) <- "CD99_2"

setChannels(cur_Images1, 1) <- cur_Images2
format(object.size(cur_Images1), units = "auto")

## -----------------------------------------------------------------------------
channels1 <- getChannels(pancreasImages_onDisk, 1:2)
channels2 <- getChannels(pancreasImages_onDisk, 3:4)

dir.create(file.path(cur_dir, "test"))
cur_path_2 <- file.path(cur_dir, "test")

channels3 <- mergeChannels(channels1, channels2,
                            h5FilesPath = cur_path_2)

seed(channels3$E34_imc)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

