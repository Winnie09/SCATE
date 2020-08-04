#' Preprocess scATAC-seq
#'
#' Preprocessing scATAC-seq samples
#'
#' This function filters out scATAC-seq with low library size and transforms the reads into middle points of the reads.
#' @param input Either a character vector of locations to bam files (when type is bam) or list of GRanges object of scATAC-seq reads (when type is gr).
#' @param type Either 'bam' or 'gr'. 'bam' if the input is locations to bam files. 'gr' if the input is a list of GRanges object.
#' @param libsizefilter Numeric variable giving the minimum library size. scATAC-seq samples with library size smaller than this cutoff will be discarded.
#' @return GRanges object of list of GRanges object after preprocessing.
#' @export
#' @import GenomicAlignments GenomicRanges
#' @author Zhicheng Ji, Weiqiang Zhou, Wenpin Hou, Hongkai Ji* <whou10@@jhu.edu>
#' @examples
#' c1 <- GRanges(seqnames="chr1",IRanges(start=seq_len(100)+1e6,end=seq_len(100)+1e6))
#' c2 <- GRanges(seqnames="chr2",IRanges(start=seq_len(100)+1e6,end=seq_len(100)+1e6))
#' grl <- list(cell1=c1,cell2=c2)
#' satacprocess(grl,type='gr',libsizefilter=10)

satacprocess <- function(input,type='bam',libsizefilter=1000) {
      if (type=='bam') {
            satac <- sapply(sapply(input,readGAlignmentPairs,simplify=FALSE),GRanges,simplify=FALSE)
      } else {
            satac <- input
      }
      satac <- satac[sapply(satac,length) >= libsizefilter]
      n <- names(satac)
      satac <- lapply(satac,function(i) {
            start(i) <- end(i) <- round((start(i) + end(i))/2)
            i
      })
      names(satac) <- n
      satac
}

