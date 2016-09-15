sanitizeReads <- function(bamFiles, cores=1, singleEnd=TRUE,
                          verbose=TRUE) {
    ## bamFiles: BamFileLists - Rsamtools
    ## only for single ended reads
    require(BiocParallel)
    require(GenomicAlignments)
    ## are bamFiles exists?
    if (!singleEnd) warning("The package has not been tested for pair-end reads")
    
    param <- Rsamtools::ScanBamParam(tag="NH", what="seq")
    prefix <- sub(".bam", "", basename(bamFiles))

    all.reads <- BiocParallel::bplapply(bamFiles, function(x) {
        if (verbose) message("Read Gap Alignments ", x)
        
        if (singleEnd)
            gr <- as(readGAlignments(x, param=param), "GRanges")
        else
            gr <- as(readGAlignmentPairs(x, param=param), "GRanges")
        ## unique read count: keep unique reads and add count column
        tmp <- unique(gr)
        label <- paste(seqnames(gr), strand(gr), start(gr), width(gr))
        label.count <- table(label)
        tmp.label <- paste(seqnames(tmp), strand(tmp), start(tmp), width(tmp))
        tmp1 <- match(tmp.label, names(label.count))
        tmp$count <- label.count[tmp1]
        reads <- tmp
        reads
    }, BPPARAM=MulticoreParam(worker=cores))
    names(all.reads) <- prefix
    all.reads
}
