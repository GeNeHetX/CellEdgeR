## ----setup, echo=FALSE--------------------------------------------------------
library(knitr)
.dpi = 100
set.seed(0)
opts_chunk$set(comment=NA, fig.align="center", dpi=.dpi)
knit_hooks$set(crop=NULL)
.output = output()
switch(.output,
        html = opts_chunk$set(fig.retina=1),
        latex = opts_chunk$set(out.width=".5\\textwidth")
)

.dev = switch(.output, html="svg", latex="pdf")
options(EBImage.display = "raster")

## ----installation, eval=FALSE-------------------------------------------------
# install.packages("BiocManager")
# BiocManager::install("EBImage")

## ----library, message=FALSE---------------------------------------------------
library("EBImage")

## ----readImage----------------------------------------------------------------
f = system.file("images", "sample.png", package="EBImage")
img = readImage(f)

## ----display------------------------------------------------------------------
display(img, method="browser")

## ----display-raster, fig.width=dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
display(img, method="raster")
text(x = 20, y = 20, label = "Parrots", adj = c(0,1), col = "orange", cex = 2)

## ----dev-print, eval=FALSE----------------------------------------------------
# filename = "parrots.jpg"
# dev.print(jpeg, filename = filename , width = dim(img)[1], height = dim(img)[2])

## ----dev-print-pre, echo=FALSE, fig.show='hide', fig.width=dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
display(img, method="raster")
text(x = 20, y = 20, label = "Parrots", adj = c(0,1), col = "orange", cex = 2)
filename = "parrots.jpg"
dev.print(jpeg, filename = filename , width = dim(img)[1], height = dim(img)[2])

## ----filesize-----------------------------------------------------------------
file.info(filename)$size

## ----dev-print3, echo=FALSE, sfig.width=dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
invisible(file.remove(filename))

## ----readImageColor-pre, echo=FALSE-------------------------------------------
imgcol = readImage(system.file("images", "sample-color.png", package="EBImage"))

## ----readImageColor, eval=FALSE-----------------------------------------------
# imgcol = readImage(system.file("images", "sample-color.png", package="EBImage"))
# display(imgcol)

## ----readImageColor-post, echo=FALSE, fig.width=dim(imgcol)[1L]/.dpi, fig.height=dim(imgcol)[2L]/.dpi, dpi=.dpi/2----
display(imgcol)

## ----readImageMulti-pre, include=FALSE----------------------------------------
nuc = readImage(system.file("images", "nuclei.tif", package="EBImage"))

## ----readImageMulti, eval=FALSE-----------------------------------------------
# nuc = readImage(system.file("images", "nuclei.tif", package="EBImage"))
# display(nuc, method = "raster", all = TRUE)

## ----readImageMulti-post, echo=FALSE, fig.width=dim(nuc)[1L]/.dpi, fig.height=dim(nuc)[2L]/.dpi, dpi=.dpi----
display(nuc, method = "raster", all = TRUE)

## ----displayFrame, echo=FALSE, fig.width=dim(nuc)[1L]/.dpi, fig.height=dim(nuc)[2L]/.dpi, dpi=.dpi/2----
display(nuc, method = "raster", frame = 2)

## ----writeImage, eval=FALSE---------------------------------------------------
# writeImage(imgcol, "sample.jpeg", quality = 85)

## ----str----------------------------------------------------------------------
str(img)

## ----dim----------------------------------------------------------------------
dim(img)

## ----imageData----------------------------------------------------------------
imageData(img)[1:3, 1:6]

## ----as.array-----------------------------------------------------------------
is.Image( as.array(img) )

## ----hist, fig.width=6, fig.height=6, dev=.dev--------------------------------
hist(img)
range(img)

## ----show---------------------------------------------------------------------
img

## ----print--------------------------------------------------------------------
print(img, short=TRUE)

## ----printcol-----------------------------------------------------------------
print(imgcol, short=TRUE)

## ----numberOfFrames-----------------------------------------------------------
numberOfFrames(imgcol, type = "render")
numberOfFrames(imgcol, type = "total")

## ----nuc----------------------------------------------------------------------
nuc

## ----colorMode, fig.width=dim(imgcol)[1L]/.dpi, fig.height=dim(imgcol)[2L]/.dpi, dpi=.dpi----
colorMode(imgcol) = Grayscale
display(imgcol, all=TRUE)

## ----Image-character, fig.width=7/.dpi, fig.height=7/.dpi, dpi=10*.dpi--------
colorMat = matrix(rep(c("red","green", "#0000ff"), 25), 5, 5)
colorImg = Image(colorMat)
colorImg
display(colorImg, interpolate=FALSE)

## ----negative, fig.width=dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
img_neg = max(img) - img
display( img_neg )

## ----arithmetic, fig.width=dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
img_comb = combine(
  img,
  img + 0.3,
  img * 2,
  img ^ 0.5
)

display(img_comb, all=TRUE)

## ----cropthreshold-pre, echo=FALSE--------------------------------------------
img_crop = img[366:749, 58:441]
img_thresh = img_crop > .5

## ----cropthreshold, eval=FALSE------------------------------------------------
# img_crop = img[366:749, 58:441]
# img_thresh = img_crop > .5
# display(img_thresh)

## ----cropthreshold-post, echo=FALSE, fig.width=dim(img_thresh)[1L]/.dpi, fig.height=dim(img_thresh)[2L]/.dpi, dpi=.dpi/2----
display(img_thresh)

## ----img_thresh---------------------------------------------------------------
img_thresh

## ----transpose, fig.width=dim(img)[2L]/.dpi, fig.height=dim(img)[1L]/.dpi, dpi=.dpi/2----
img_t = transpose(img)
display( img_t )

## ----translate, fig.width=dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
img_translate = translate(img, c(100,-50))
display(img_translate)

## ----rotate-pre, echo=FALSE---------------------------------------------------
img_rotate = rotate(img, 30, bg.col = "white")

## ----rotate, eval=FALSE-------------------------------------------------------
# img_rotate = rotate(img, 30, bg.col = "white")
# display(img_rotate)

## ----rotate-post, echo=FALSE, fig.width=dim(img_rotate)[1L]/.dpi, fig.height=dim(img_rotate)[2L]/.dpi, dpi=.dpi/2----
display(img_rotate)

## ----resize-pre, echo=FALSE---------------------------------------------------
img_resize = resize(img, w=256, h=256)

## ----resize, eval=FALSE-------------------------------------------------------
# img_resize = resize(img, w=256, h=256)
# display(img_resize )

## ----resize-post, echo=FALSE, fig.width=dim(img_resize)[1L]/.dpi, fig.height=dim(img_resize)[2L]/.dpi, dpi=.dpi/2----
display(img_resize)

## ----flipflop, fig.width=2*dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
img_flip = flip(img)
img_flop = flop(img)

display(combine(img_flip, img_flop), all=TRUE)

## ----affine, fig.width=dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
m =  matrix(c(1, -.5, 128, 0, 1, 0), nrow=3, ncol=2)
img_affine = affine(img, m)
display( img_affine )

## ----makeBrush, fig.width=6, fig.height=6, dev=.dev---------------------------
w = makeBrush(size = 31, shape = 'gaussian', sigma = 5)
plot(w[(nrow(w)+1)/2, ], ylab = "w", xlab = "", cex = 0.7)

## ----lopass, fig.width=dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
img_flo = filter2(img, w)
display(img_flo)

## ----gblur, fig.width=dim(nuc)[1L]/.dpi, fig.height=dim(nuc)[2L]/.dpi, dpi=.dpi----
nuc_gblur = gblur(nuc, sigma = 5)
display(nuc_gblur, all=TRUE )

## ----highpass, fig.width=dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
fhi = matrix(1, nrow = 3, ncol = 3)
fhi[2, 2] = -8
img_fhi = filter2(img, fhi)
display(img_fhi)

## ----medianFilter, fig.width=dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
l = length(img)
n = l/10
pixels = sample(l, n)
img_noisy = img
img_noisy[pixels] = runif(n, min=0, max=1)
display(img_noisy)
img_median = medianFilter(img_noisy, 1)
display(img_median)

## ----logo-pre, echo=FALSE-----------------------------------------------------
shapes = readImage(system.file('images', 'shapes.png', package='EBImage'))
logo = shapes[110:512,1:130]

## ----logo, eval=FALSE---------------------------------------------------------
# shapes = readImage(system.file('images', 'shapes.png', package='EBImage'))
# logo = shapes[110:512,1:130]
# display(logo)

## ----logo-post, echo=FALSE, fig.width=dim(logo)[1L]/.dpi, fig.height=dim(logo)[2L]/.dpi, dpi=.dpi----
display(logo)

## ----kern, fig.width=7/.dpi, fig.height=7/.dpi, dpi=10*.dpi-------------------
kern = makeBrush(5, shape='diamond')
display(kern, interpolate=FALSE)

## ----morph, fig.width=2*dim(logo)[1L]/.dpi, fig.height=dim(logo)[2L]/.dpi, dpi=.dpi----
logo_erode= erode(logo, kern)
logo_dilate = dilate(logo, kern)

display(combine(logo_erode, logo_dilate), all=TRUE)

## ----otsu, fig.width=dim(nuc)[1L]/.dpi, fig.height=dim(nuc)[2L]/.dpi, dpi=.dpi----
threshold = otsu(nuc)
threshold
nuc_th = combine( mapply(function(frame, th) frame > th, getFrames(nuc), threshold, SIMPLIFY=FALSE) )
display(nuc_th, all=TRUE)

## ----filter2thresh, fig.width=dim(nuc)[1L]/.dpi, fig.height=dim(nuc)[2L]/.dpi, dpi=.dpi----
disc = makeBrush(31, "disc")
disc = disc / sum(disc)
offset = 0.05
nuc_bg = filter2( nuc, disc )
nuc_th = nuc > nuc_bg + offset
display(nuc_th, all=TRUE)

## ----thresh, fig.width=dim(nuc)[1L]/.dpi, fig.height=dim(nuc)[2L]/.dpi, dpi=.dpi----
display( thresh(nuc, w=15, h=15, offset=0.05), all=TRUE )

## ----bwlabel------------------------------------------------------------------
logo_label = bwlabel(logo)
table(logo_label)

## ----max_logolabel------------------------------------------------------------
max(logo_label)

## ----displaybw, fig.width=dim(logo)[1L]/.dpi, fig.height=dim(logo)[2L]/.dpi, dpi=.dpi----
display( normalize(logo_label) )

## ----colorCode, fig.width=dim(logo)[1L]/.dpi, fig.height=dim(logo)[2L]/.dpi, dpi=.dpi----
display( colorLabels(logo_label) )

## ----watershed, fig.width=dim(nuc)[1L]/.dpi, fig.height=dim(nuc)[2L]/.dpi, dpi=.dpi----
nmask = watershed( distmap(nuc_th), 2 )
display(colorLabels(nmask), all=TRUE)

## ----voronoiExample, fig.width=dim(nuc)[1L]/.dpi, fig.height=dim(nuc)[2L]/.dpi, dpi=.dpi----
voronoiExamp = propagate(seeds = nmask, x = nmask, lambda = 100)
voronoiPaint = colorLabels (voronoiExamp)
display(voronoiPaint)

## ----rmObjects, fig.width=2*dim(logo)[1L]/.dpi, fig.height=dim(logo)[2L]/.dpi, dpi=.dpi----
objects = list(
    seq.int(from = 2, to = max(logo_label), by = 2),
    seq.int(from = 1, to = max(logo_label), by = 2)
    )
logos = combine(logo_label, logo_label)
z = rmObjects(logos, objects, reenumerate=FALSE)
display(z, all=TRUE)

## ----uniqueIDs----------------------------------------------------------------
showIds = function(image) lapply(getFrames(image), function(frame) unique(as.vector(frame)))

showIds(z)

## ----reenumeratedIDs----------------------------------------------------------
showIds( reenumerate(z) )

## ----fillHull, fig.width=dim(logo)[1L]/.dpi, fig.height=dim(logo)[2L]/.dpi, dpi=.dpi----
filled_logo = fillHull(logo)
display(filled_logo)

## ----floodFill-logo, fig.width=dim(logo)[1L]/.dpi, fig.height=dim(logo)[2L]/.dpi, dpi=.dpi----
rgblogo = toRGB(logo)
points = rbind(c(50, 50), c(100, 50), c(150, 50))
colors = c("red", "green", "blue")
rgblogo = floodFill(rgblogo, points, colors)
display( rgblogo )

## ----floodFill-img, fig.width=dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
display(floodFill(img, rbind(c(200, 300), c(444, 222)), col=0.2, tolerance=0.2))

## ----paintObjects, fig.width=dim(img)[1L]/.dpi, fig.height=dim(img)[2L]/.dpi, dpi=.dpi/2----
d1 = dim(img)[1:2]
overlay = Image(dim=d1)
d2 = dim(logo_label)-1

offset = (d1-d2) %/% 2

overlay[offset[1]:(offset[1]+d2[1]), offset[2]:(offset[2]+d2[2])] = logo_label

img_logo = paintObjects(overlay, toRGB(img), col=c("red", "yellow"), opac=c(1, 0.3), thick=TRUE)

display( img_logo )

## ----load, fig.height=dim(nuc)[2L]/.dpi, fig.width=dim(nuc)[1L]/.dpi, warning=FALSE, dpi=.dpi----
nuc = readImage(system.file('images', 'nuclei.tif', package='EBImage'))
cel = readImage(system.file('images', 'cells.tif', package='EBImage'))

cells = rgbImage(green=1.5*cel, blue=nuc)
display(cells, all = TRUE)

## ----nmask, fig.width=dim(nuc)[1L]/.dpi, fig.height=dim(nuc)[2L]/.dpi, dpi=.dpi----
nmask = thresh(nuc, w=10, h=10, offset=0.05)
nmask = opening(nmask, makeBrush(5, shape='disc'))
nmask = fillHull(nmask)
nmask = bwlabel(nmask)

display(nmask, all=TRUE)

## ----ctmask, fig.width=dim(nuc)[1L]/.dpi, fig.height=dim(nuc)[2L]/.dpi, dpi=.dpi----
ctmask = opening(cel>0.1, makeBrush(5, shape='disc'))
cmask = propagate(cel, seeds=nmask, mask=ctmask)

display(ctmask, all=TRUE)

## ----res, fig.width=dim(nuc)[1L]/.dpi, fig.height=dim(nuc)[2L]/.dpi, dpi=.dpi----
segmented = paintObjects(cmask, cells, col='#ff00ff')
segmented = paintObjects(nmask, segmented, col='#ffff00')

display(segmented, all=TRUE)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

