#'
#' adapted from GenomicAlignment::summarizeOverlaps. Same approach but
#' the counts are adjusted by NH: number of reported alignment.
#' Only support mode: IntersectionStrict and inter.feature:TRUE/FALSE
#' 

.countSubjectHits.2 <- function(ov, features, reads) {
    fracReadCount <- reads[queryHits(ov)]$count/reads[queryHits(ov)]$NH
    cnt <- rep(0, length(features))
    names(cnt) <- names(features)
    tmp <- tapply(fracReadCount, subjectHits(ov), sum)
    cnt[as.integer(names(tmp))] <- tmp
    as.integer(round(cnt))
}
    
countHits <- function(features, reads, type=c("any", "start", "end", "within"),
                  ignore.strand=TRUE, inter.feature=FALSE)
{   
    ov <- findOverlaps(reads, features, type=match.arg(type),
                       ignore.strand=ignore.strand)
    if (inter.feature) {
        ## Remove ambigous reads.
        reads_to_keep <- which(countQueryHits(ov) == 1L)
        ov <- ov[queryHits(ov) %in% reads_to_keep]
    }
    .countSubjectHits.2(ov, features=features, reads=reads)
}

summarizeOverlaps.adjNH <- function(features, all.reads, cores=1,
                                    type=c("any", "start", "end", "within"),
                                    ignore.strand=TRUE, inter.feature=FALSE) {
    ## features: GRangeList only
    ## all.reads: A list of GRanges or GRranges only (all.reads), the sanitize reads
    ## type: any -> union, within -> IntersectStrick
    type <- match.arg(type)
    if (any(is.na(names(all.reads))))
        stop("The list of GRanges (reads) must have names. Asigne names to the list.")
    
    counts <- BiocParallel::bplapply(all.reads, function(reads) {
        countHits(features=features, reads=reads,
                  type=type,
                  ignore.strand=ignore.strand, inter.feature=inter.feature)
    },  BPPARAM=MulticoreParam(worker=cores))

    counts <- as.matrix(do.call(cbind, counts))
    SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)),
                         rowRanges=features)
}


