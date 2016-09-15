# rmskStats

## Overview
This is a standard R package providing tools to perform hit counts and enrichment/depletion analysis for repeat elements annotated by 
Repeat Master track (RMSK). The user must be familiar with class types and sequence analysis tools supported by Bioconductor such as 
the SummerizedExperiment object and DESeq2 package.

The hit counts for repeat elements (repName) is similiar to the summarizeOverlaps() function from the GenomicAlignments with extra 
adjustment w.r.t. numbers of reported multiple alignment (NH column).  *DESeq2* is used to identify the 
up/down-regulated repeat elements (repName). The enrichment/depletion analysis uses Hypergeometric test to determ the significance
of the enrichment state of a repeat family (repFamily) or class (repClass).

This package is under development and will soon be submit to Bioconductor repository.
