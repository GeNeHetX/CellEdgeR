## ----style, echo=FALSE--------------------------------------------------------
knitr::opts_chunk$set(error=FALSE, warning=FALSE, message=FALSE)
library(BiocStyle)

## ----echo=FALSE, results='hide'-----------------------------------------------
marker_names <- data.frame(
    full_name = c(
        "Carbonic anhydrase IX",
        "CD3 epsilon",
        "CD8 alpha",
        "E-Cadherin",
        "cleaved-Caspase3 + cleaved-PARP",
        "Cytokeratin 5",
        "Forkhead box P3",
        "Glucose transporter 1",
        "Histone H3",
        "phospho-Histone H3 [S28]",
        "Ki-67",
        "Myeloperoxidase",
        "Programmed cell death protein 1",
        "Programmed death-ligand 1",
        "phospho-Rb [S807/S811]",
        "Smooth muscle actin",
        "Vimentin",
        "Iridium 191",
        "Iridium 193"
    ),
    short_name = c(
        "CA9",
        "CD3e",
        "CD8a",
        "CDH1",
        "cCASP3_cPARP",
        "KRT5",
        "FOXP3",
        "SLC2A1",
        "H3",
        "p_H3",
        "Ki67",
        "MPO",
        "PD_1",
        "PD_L1",
        "p_Rb",
        "SMA",
        "VIM",
        "DNA1",
        "DNA2"
    )
)

## ----echo=FALSE, results='asis'-----------------------------------------------
knitr::kable(
    marker_names,
    caption = "'full_name' and 'short_name' examples for some commonly 
        used markers"
)

## ----sessionInfo, echo=FALSE--------------------------------------------------
sessionInfo()

