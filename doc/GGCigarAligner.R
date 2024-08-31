## ----setup--------------------------------------------------------------------
library(GGCigarAligner)

## -----------------------------------------------------------------------------
df_gg <- ggBamLoader(system.file("extdata", "subset.bam", package = "GGCigarAligner"))
df_gg[1:5,]

## ----warning=FALSE, message=FALSE---------------------------------------------
result <- ggCigarAligner(df_gg, qname = "ERR188273.4711308", index = 3, my_reference = "BSgenome.Hsapiens.UCSC.hg38")
result

## -----------------------------------------------------------------------------
sessionInfo()

