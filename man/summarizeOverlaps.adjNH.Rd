\name{summarizeOverlaps.adjNH}
\alias{summarizeOverlaps.adjNH}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Permorm overlap queries between reads and repetitive genomic features
}
\description{
An extentions of 'summarizeOverlaps' by adjuesting the hit counts by the number of reported alignments (NH tag).
}
\usage{
summarizeOverlaps.adjNH(features, all.reads, cores = 1,
                        type = c("any", "start", "end", "within"),
			ignore.strand = TRUE, inter.feature = FALSE)
}

\arguments{
  \item{features}{
A GRanges or a GRangesList object of genomic regions of interest. 
}
  \item{all.reads}{
A list of GRanges consisting of unique reads with metadata column "count" and "NH". It represents the data to be counted by 'summarizeOverlaps.adjNH'. One can generate this object by using 'sanitizeReads'.
}
  \item{cores}{

}
  \item{type}{
Default 'any'. It is the pre-defined count methods. 'any' is equivalent to "Union" and "within" to "IntersectionStrick" method. 
}
  \item{ignore.strand}{
Default TURE. This function is only tested for single-ended reads. Be aware if using FALSE, the results are not gaurentee to be correct.
}
  \item{inter.feature}{
Default FALSE. A logical indicating if the counting 'type' should be aware of overlapping features. When TRUE, reads map to multiple features are droped. When FALSE, reads are retained and a count is assigned to each feature they map to.
}
}


\value{SummarizedExpriment}

\author{
Chao-Jen Wong
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{sanitizeReads}}
\code{\link{summarizeOverlaps}{GenomicAlignment}}
}
\examples{
library(SGSeq)
bam_dir <- system.file("extdata", "bams", package="SGSeq")
bam_files <- list.files(bam_dir, pattern="\\.bam$", full.names=TRUE)
all.reads <- sanitizeReads(bamFiles=bam_files, cores=2, singleEnd=TRUE,
                           verbose=TRUE)
all.reads[[1]]

data(sub.hg19.rmsk)
seqlevelsStyle(sub.hg19.rmsk) <- "NCBI"
keeplevels <- intersect(seqlevels(sub.hg19.rmsk), seqlevels(all.reads[[1]]))
sub.hg19.rmsk <- keepSeqlevels(sub.hg19.rmsk, keeplevels)
repElements <- split(sub.hg19.rmsk, sub.hg19.rmsk$repName)
se.repeat <- summarizeOverlaps.adjNH(features=repElements,
                                     all.reads=all.reads,
                                     type="any",
                                     ignore.strand=TRUE,
                                     inter.feature=FALSE)
}

