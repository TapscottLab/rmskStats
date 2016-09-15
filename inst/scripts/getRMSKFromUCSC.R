dataDir <- "~/tapscott/R_package/hiyaRepeats/data"

getRMSKFromUCSC <- function(genome) {
    mySession = browserSession("UCSC")
    genome(mySession) <- genome
    tbl.rmsk <- getTable(ucscTableQuery(mySession, track="rmsk", table="rmsk"))

    rmsk <- GRanges(seqnames=Rle(tbl.rmsk$genoName),
                    ranges=IRanges(start=tbl.rmsk$genoStart,
                         end=tbl.rmsk$genoEnd),
                    strand=tbl.rmsk$strand)
    n <- setdiff(colnames(tbl.rmsk),
                 c("genoName", "genoStart", "genoEnd", "strand"))
    mcols(rmsk) <- tbl.rmsk[, n]
    rmsk
}

assignElementID <- function(rmsk, cores=1) {
    rmsk$dummyID <- seq(1, length(rmsk))
    repName <- split(rmsk, rmsk$repName, drop=FALSE)
    repName <- parallel::mclapply(repName, function(x) {
        idx <- seq(1, length(x))
        x$elementID <- sprintf("%s_%s_%s", seqnames(x),
                             x$repName, idx)
        x
    }, mc.cores=cores)
    tmp <- repName
    names(tmp) <- NULL
    tmp <- do.call(c, tmp)
    tmp <- tmp[order(tmp$dummyID)]
    ## TODO: take out dummy ID 
    tmp
}
#'
#' hg38 gets repeat masker track from ucsc: 4/12/2016
#'
library(rtracklayer)
library(GenomicRanges)
hg38.rmsk <- getRMSKFromUCSC(genome="hg38")
save(hg38.rmsk, file=file.path(dataDir, "hg38.rmsk.rda"))


#'
#' mm10
#'
library(rtracklayer)
library(GenomicRanges)
mySession = browserSession("UCSC")
genome(mySession) <- "mm10"
tbl.rmsk <- getTable(ucscTableQuery(mySession, track="rmsk", table="rmsk"))

mm10.rmsk <- GRanges(seqnames=Rle(tbl.rmsk$genoName),
                     ranges=IRanges(start=tbl.rmsk$genoStart,
                         end=tbl.rmsk$genoEnd),
                     strand=tbl.rmsk$strand)
n <- setdiff(colnames(tbl.rmsk), c("genoName", "genoStart", "genoEnd", "strand"))
mcols(mm10.rmsk) <- tbl.rmsk[, n]
save(mm10.rmsk, file=file.path(dataDir, "mm10.rmsk.rda"))

#'
#' mm10
#'
library(rtracklayer)
library(GenomicRanges)
mySession = browserSession("UCSC")
genome(mySession) <- "mm9"
tbl.rmsk <- getTable(ucscTableQuery(mySession, track="rmsk", table="rmsk"))

mm9.rmsk <- GRanges(seqnames=Rle(tbl.rmsk$genoName),
                     ranges=IRanges(start=tbl.rmsk$genoStart,
                         end=tbl.rmsk$genoEnd),
                     strand=tbl.rmsk$strand)
n <- setdiff(colnames(tbl.rmsk), c("genoName", "genoStart", "genoEnd", "strand"))
mcols(mm9.rmsk) <- tbl.rmsk[, n]
save(mm9.rmsk, file=file.path(dataDir, "mm9.rmsk.rda"))

#'
#' hg19
#'
library(rtracklayer)
library(GenomicRanges)
mySession = browserSession("UCSC")
genome(mySession) <- "hg19"
tbl.rmsk <- getTable(ucscTableQuery(mySession, track="rmsk", table="rmsk"))

hg19.rmsk <- GRanges(seqnames=Rle(tbl.rmsk$genoName),
                     ranges=IRanges(start=tbl.rmsk$genoStart,
                         end=tbl.rmsk$genoEnd),
                     strand=tbl.rmsk$strand)
n <- setdiff(colnames(tbl.rmsk), c("genoName", "genoStart", "genoEnd", "strand"))
mcols(hg19.rmsk) <- tbl.rmsk[, n]
save(hg19.rmsk, file=file.path(dataDir, "hg19.rmsk.rda"))
