\name{sanitizeReads}
\alias{sanitizeReads}
\title{
Import Reads from BAM Files and Convert to a list of GRanges Consisting with Unique Reads with metadata "count" and "NH".
}
\description{
This step is a prerequisite of counting reads for genomic features, especially for repetitive elements. It import Gap alignments from each of the BAM files and  converts to GRanges consisting of unique reads "count" and "NH" metadata.

}
\usage{
sanitizeReads(bamFiles, cores = 1, singleEnd = TRUE, verbose = TRUE)
}

\arguments{
  \item{bamFiles}{
A vector of charactors indicating the path of the BAM files.
}
  \item{cores}{
The number of available multi-processors.
}
  \item{singleEnd}{
Default TRUE. This function only support single-ended read.
}
  \item{verbose}{
Default TRUE.
}
}

\value{ A list of GRanges with metadata column "count" and "NH"}

\references{
%% ~put references to the literature/web site here ~
}
\author{
Zizhen Yao and Chao-Jen Wong
}

\seealso{
\code{\link{summarizeOverlaps.adjNH}}
}
\examples{
library(SGSeq)
bam_dir <- system.file("extdata", "bams", package="SGSeq")
bam_files <- list.files(bam_dir, pattern="\\.bam$", full.names=TRUE)
all.reads <- sanitizeReads(bamFiles=bam_files, cores=2, singleEnd=TRUE,
                           verbose=TRUE)
}
