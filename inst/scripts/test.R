library(SGSeq)
library(countRepeats)
pkgDir <- "~/tapscott/R_package/countRepeats"
bam_dir <- system.file("extdata", "bams", package="SGSeq")
bam_files <- list.files(bam_dir, pattern="\\.bam$", full.names=TRUE)
all.reads <- sanitizeReads(bamFiles=bam_files, cores=10, singleEnd=TRUE,
                           verbose=TRUE)
save(all.reads, file=file.path(pkgDir, "data", "all.reads.rda"))

## general examples
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
features <- exonsBy(TxDb.Hsapiens.UCSC.hg19.knownGene, by="gene")
seqlevelsStyle(features) <- "NCBI"
keeplevels <- intersect(seqlevels(features), seqlevels(all.reads[[1]]))
features <- keepSeqlevels(features, keeplevels)

se <- summarizeOverlaps.adjNH(features=features, cores=10,
                              all.reads=all.reads,
                              type="any",
                              ignore.strand=TRUE,
                              inter.feature=FALSE)
### example
data(hg19.rmsk)
seqlevelsStyle(hg19.rmsk) <- "NCBI"
keeplevels <- intersect(seqlevels(hg19.rmsk), seqlevels(all.reads[[1]]))
hg19.rmsk <- keepSeqlevels(hg19.rmsk, keeplevels)
repElements <- split(hg19.rmsk, hg19.rmsk$repName)
se.repeat <- summarizeOverlaps.adjNH(features=repElements,
                                     all.reads=all.reads,
                                     type="any",
                                     ignore.strand=TRUE,
                                     inter.feature=FALSE)
