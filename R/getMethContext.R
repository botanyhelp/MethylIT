#' @rdname getMethContext
#' @title Get Methylation Context from a chromosome DNA sequence
#' @description This function retrieves the methylation context from a
#'     chromosome DNA sequence in fasta format.
#' @param chr.seq DNA sequence from a chromosome in fasta format.
#' @param chromosome Chromosome name.
#' @param verbose If TRUE, prints the function log to stdout
#' @return GRanges object with three columns: 'trinucleotide', methylation
#'     context, and 'CHH' methylation subcontexts: 'CHA', 'CHC', and 'CHT'.
#'
#' @examples
#' dna <- Biostrings::DNAString(x = "CCCTAACGACCCTAACGCTACCCTAAACCTCTGAAT",
#'     start = 1, nchar = NA)
#' getMethContext(chr.seq = dna, chromosome = "1", verbose = TRUE)
#'
#'
#' @importFrom Biostrings complement
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
#'
getMethContext <- function(chr.seq, chromosome, verbose=TRUE){

   base <- strsplit(unname(as.character(chr.seq)), split="")[[1]]
   start <- grep("[C]", base)
   trinucleotide <- paste0(base[start], base[(start + 1)], base[(start + 2)])
   context <- trinucleotide

   ## ==== Positive strand ===== ##
   if (verbose) cat("*** Extracting contexts for the positive strand ...\n")
   context[grep("CG", trinucleotide)] <- "CG"
   context[grep("C[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}G", trinucleotide)] <- "CHG"
   subcontext <- context
   context[grep("C[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}",
           trinucleotide )] <- "CHH"
   subcontext[grep("C[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}A", trinucleotide)] <- "CHA"
   subcontext[grep("C[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}C", trinucleotide)] <- "CHC"
   subcontext[grep("C[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}T", trinucleotide)] <- "CHT"
   subcontext[grep("C[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}[BDEFHIJKLMNOPQRSUVWXYZ]{1}",
           trinucleotide)] <- "CHH"

   x.pos <- data.frame(seqnames=chromosome, start=start, end=start,
                     strand="+", trinucleotide=trinucleotide,
                     context=context, subcontext=subcontext)

   ## ==== negative strand ===== ##
   if (verbose) cat("*** Extracting contexts for the negative strand ...\n")
   chr.seq <- complement(chr.seq)
   base <- strsplit(unname(as.character(chr.seq)), split="")[[1]]
   start <- grep("[C]", base)
   trinucleotide <- paste0(base[start], base[(start - 1)], base[(start - 2)])
   context <- trinucleotide
   context[grep("CG", trinucleotide)] <- "CG"
   context[grep("C[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}G", trinucleotide)] <- "CHG"
   subcontext <- context
   context[grep("C[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}",
           trinucleotide)] <- "CHH"
   subcontext[grep("C[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}A", trinucleotide)] <- "CHA"
   subcontext[grep("C[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}C", trinucleotide)] <- "CHC"
   subcontext[grep("C[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}T", trinucleotide)] <- "CHT"
   subcontext[grep("C[ABCDEFHIJKLMNOPQRSTUVWXYZ]{1}[BDEFHIJKLMNOPQRSUVWXYZ]{1}",
           trinucleotide)] <- "CHH"
   x.neg <- data.frame(seqnames=chromosome, start=start, end=start, strand="-",
           trinucleotide=trinucleotide, context=context, subcontext=subcontext)
   x <- rbind(x.pos, x.neg)
   x <- makeGRangesFromDataFrame(x, keep.extra.columns=TRUE)
   return(sortBySeqnameAndStart(x))
}
