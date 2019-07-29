#' @rdname readCounts2GRangesList
#'
#' @title Read files of methylation count tables
#' @description This function is addressed to read files with methylation count
#'     table data commonly generated after the alignment of BS-seq data or found
#'     in GEO database
#' @details Read tables from files with a table methylation count data using
#'     the function fread from the package 'data.table' and and yields a list of
#'     GRanges objects with the information provided.
#'
#' @param filenames Character vector with the file names
#' @param sample.id Character vector with the names of the samples
#'     corresponding to each file
#' @param pattern Chromosome name pattern. Users working on Linux OS can
#'     specify the reading of specific lines from each file by using regular
#'     expressions.
#' @param remove Logic (TRUE). Usually the supplementary files from GEO
#'     datasets are 'gz' compressed. File datasets must be decompressed to be
#'     read. The decompressed files are removed after read if this is set
#'     'TRUE'.
#' @param columns Vector of integer numbers denoting the table columns that
#'     must be read. The numbers for 'seqnames' (chromosomes), 'start', and
#'     'end' (if different from 'start') columns must be given. The possible
#'     fields are: 'seqnames' (chromosomes),'start', 'end', 'strand',
#'     'fraction', percent' (metylation percentage), 'mC' (methylates cytosine),
#'     'uC' (non methylated cytosine), 'coverage', and 'context'
#'     (methylation context). These column headers are not required to be in the
#'     files.
#' @param chromosome.names If provided, for each GRanges object, chromosome
#'     names will be changed to those provided in 'chromosome.names' applying
#'     seqlevels(x) <- chromosome.names'. This option permits to use all the
#'     functionality of the function "seqlevels" defined from package
#'     "GenomeInfoDb", which rename, add, and reorder the seqlevels all at once
#'     (see ?seqlevels).
#' @param chromosomes If provided, it must be a character vector with the names
#'     of the chromosomes that you want to include in the final GRanges objects.
#' @param contexts If provided, it must be a character vector with the names
#'     of the contexts that you want to include in the final GRanges objects, 
#'     for example, "CG", "CHG"
#' @param seqs If provided, it must be a character vector with the names
#'     of the trinucleotide sequences that you want to include in the final 
#'     GRanges objects, for example, "CGA", "CAA"
#' @param verbose If TRUE, prints the function log to stdout
#' @param ... Additional parameters for 'fread' function from 'data.table'
#'     package
#'
#' @return A list of GRanges objects 
#'
#' @examples
#' ## Create a cov file with it's file name including "gz" (tarball extension)
#' filename <- "./file.cov"
#' gr1 <- data.frame(chr = c("chr1", "chr1"), post = c(1,2),
#'                 strand = c("+", "-"), ratio = c(0.9, 0.5),
#'                 context = c("CG", "CG"), CT = c(20, 30))
#' filename <- "./file.cov"
#' write.table(as.data.frame(gr1), file = filename,
#'             col.names = TRUE, row.names = FALSE, quote = FALSE)
#'
#' ## Read the file. It does not work. Typing mistake: "fractions"
#' LR <- try(readCounts2GRangesList(filenames = filename, remove = FALSE,
#'                             sample.id = "test",
#'                             columns = c(seqnames = 1, start = 2,
#'                                     strand = 3, fractions = 4,
#'                                     context = 5, coverage = 6)),
#'                                     silent = TRUE)
#' file.remove(filename) # Remove the file
#'
#' ## Read the file
#' ## Create a cov file with it's file name including "gz" (tarball extension)
#' filename <- "./file.cov"
#' gr1 <- data.frame(chr = c("chr1", "chr1"), post = c(1,2),
#'                 strand = c("+", "-"), ratio = c(0.9, 0.5),
#'                 context = c("CG", "CG"), CT = c(20, 30))
#' filename <- "./file.cov"
#' write.table(as.data.frame(gr1), file = filename,
#'             col.names = TRUE, row.names = FALSE, quote = FALSE)
#'
#' LR <- readCounts2GRangesList(filenames = filename, remove = TRUE,
#'                              sample.id = "test",
#'                              columns = c(seqnames = 1, start = 2,
#'                                          strand = 3, fraction = 4,
#'                                          context = 5, coverage = 6))
#'
#' ## Download supplementary files from GEO data set and store "fullpath/name"
#' ## in variable filename. The parameter 'pattern' permits us to download only
#' ## the specified filesCG, in this case, CG and CHG methylation contexts.
#'
#' filenames <- getGEOSuppFiles(GEO = "GSM881757",
#'                             pattern = "G_cytosine.txt.gz")
#'
#' ## Read the files with function 'readCounts2GRangesList'. Only lines starting
#' ## with the word 'Chr1' will be read, in acccordance with the  specification
#' ##given with parameter 'pattern'
#'
#' LR <- readCounts2GRangesList(filenames = filenames, remove = TRUE,
#'                             sample.id = c("drm2_CG", "drm2_CHG"),
#'                             columns = c(seqnames = 1, start = 2,
#'                                         mC = 4, uC = 3),
#'                             pattern = "^Chr1", verbose = TRUE)
#' file.remove(filenames) # Remove the downloaded file
#' 
#' #' filename <- "./file.cov"
#' gr1 <- data.frame(chr = c("chr1", "chr1"), post = c(1,2),
#'                 strand = c("+", "-"), ratio = c(0.9, 0.5),
#'                 context = c("CG", "CG"), CT = c(20, 30))
#' filename <- "./file.cov"
#' write.table(as.data.frame(gr1), file = filename,
#'             col.names = TRUE, row.names = FALSE, quote = FALSE)
#'
#' ## Read the file. It does not work. Typing mistake: "fractions"
#' LR <- try(readCounts2GRangesList(filenames = filename, remove = FALSE,
#'                             sample.id = "test",
#'                             columns = c(seqnames = 1, start = 2,
#'                                     strand = 3, fractions = 4,
#'                                     context = 5, coverage = 6)),
#'                                     silent = TRUE)
#' file.remove(filename) # Remove the file
#' 
#' c("seqnames", "start", "strand", "mC", "uC", "context", "seq")
#' filename <- "./file.cov"
#' gr1 <- data.frame(seqnames = rep(1:20,each=5), start = seq(10, 1000, by=10),
#'                 strand = rep(c("+", "-"),times=25,each=2), 
#'                 mC=sample(0:10,size=100,replace=T),
#'                 uC=sample(0:10,size=100,replace=T),
#'                 context = rep(c("CG","CHG","CHH"),length.out=100),
#'                 seq = rep(c("CGC", "CAG", "CTT"),length.out=100))
#' filename <- "./file.cov"
#' write.table(as.data.frame(gr1), file = filename,
#'             col.names = TRUE, row.names = FALSE, quote = FALSE)
#'
#' ## Read the file. It does not work. Typing mistake: "fractions"
#' LR <- try(readCounts2GRangesList(filenames = filename, remove = FALSE,
#'                             sample.id = "test",
#'                             columns = c(seqnames = 1, start = 2,
#'                                     strand = 3, mC = 4, uC = 5,
#'                                     context = 6, seq = 7)),
#'                                     silent = TRUE)
#' file.remove(filename) # Remove the file
#' 
#'
#' @importFrom data.table fread
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#'
#' @export
readCounts2GRangesList <- function(filenames=NULL, sample.id=NULL, 
                   pattern=NULL, remove=FALSE, columns=c(seqnames=NULL, 
                   start=NULL, end=NULL, strand=NULL,fraction=NULL, 
                   percent=NULL, mC=NULL, uC=NULL, coverage=NULL, context=NULL, 
                   seq=NULL), chromosome.names=NULL, chromosomes=NULL, 
                   contexts=NULL, seqs=NULL, verbose=TRUE, ...) {

   cn <- names(columns)

   colns <- c("seqnames", "start", "end", "strand", "fraction",
            "percent", "mC", "uC", "coverage", "context", "seq")
   
   Check <- ArgumentCheck::newArgCheck()
   if (is.null(filenames)) {
       ArgumentCheck::addError(msg="No file names provided",
                           argcheck=Check)
   }
   for (file in filenames) {
       if (!file.exists(file)) {
           ArgumentCheck::addError(msg=paste0("Unable to find: ",
                                           file), argcheck=Check)
       }
   }

   cn <- names(columns)

   if (!is.element("seqnames", cn) || !is.element("start", cn)) {
       ArgumentCheck::addError(msg=paste0("You must provide the numbers for ",
                   "'seqnames' (chromosomes names), ", "  'start' columns"),
                   argcheck=Check)
   }

   gz.ext = "[.]gz$"
   if (grepl(gz.ext, "", filenames) &&
       sum(is.element(filenames, sub(gz.ext, "", filenames))) > 0) {
       ArgumentCheck::addError(msg=paste0("File duplication. File: \n",
                   filenames[is.element(filenames, sub(gz.ext, "", filenames))],
                   "\n is also provided in compressed format '.gz' \n"),
                   argcheck=Check)
   }

   idx <- is.element(names(columns), colns)
   if (sum(idx) != length(columns)) {
       ArgumentCheck::addError(msg=paste0("Probably you have a typing mistake ",
           "in the column names: ", "'", names(columns)[!idx],"'"),
           argcheck=Check)
   }

   ArgumentCheck::finishArgCheck(Check)
   
   x_gr_list <- vector(mode = "list", length = length(filenames))
   
   for (k in 1:length(filenames)) {
       if (.Platform$OS.type == "unix" && (!is.null(pattern))) {
               x <- fread(paste0("egrep ", pattern, " ", filenames[k]),
                  select=columns, ...)
       } else {
           x <- fread (file = filenames[k], select=columns, ...)
       }

       # colnames(x) <-
       #     c("seqnames", "start", "strand", "mC", "uC", "context", "seq")
       colnames(x) <- cn
       
       if (!is.element("end", cn)) x$end <- x$start
       if (is.element("coverage", cn) && is.element("mC", cn)) {
           x$uC <- x$coverage - x$mC
       }
       if (is.element("fraction", cn) && is.element("coverage", cn)) {
           x$mC <- round(x$fraction * x$coverage)
           x$uC <- round(x$coverage) - x$mC
       }
       if (is.element("percent", cn) && is.element("coverage", cn)) {
           x$mC <- round(x$coverage * x$percent / 100)
           x$uC <- round(x$coverage) - x$mC
       }

       x_gr <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
       
       if (!is.null(chromosomes)) {
           message("*** chromosomes is NOT NULL, changing seqlevels..")
           seqlevels(x_gr, pruning.mode = "coarse") <- chromosomes
       }

       x_gr <- sortBySeqnameAndStart(x_gr)
       x_gr_list[[k]] <- x_gr
   }
   
   if (!is.null(chromosome.names)) {
       if (verbose) message("*** chromosome.names is NOT NULL..")
       if (length(unique(unlist(lapply(
           lapply(x_gr_list, seqlevels), length
       )))) == 1) {
           if (verbose) {
               message("*** seqlevels of all objects in list are equal")
           }
           if (unique(unlist(lapply(lapply(
               x_gr_list, seqlevels
           ), length))) == length(chromosome.names)) {
               if (verbose) { message(paste0("*** length of seqlevels of all ",
                            "objects in list are equal ",
                            "to length(chromosome.names), calling seqlevels.."))
               }
               for(i in 1:length(x_gr_list)){
                   if(verbose) { message(paste0("*** setting seqlevels of ",
                               " list item ",i," to new seqlevels"))
                   }
                   seqlevels(x_gr_list[[i]]) <- chromosome.names
               }
           }
       } else {
           message(paste0("*** seqlevels of all objects in list must be ",
                  "equal. They are not.  Therefore, skipping setting ",
                  "seqlevels to chromosome.names"))
       }
   } else {
       if (verbose) { message(paste0("*** chromosome.names is NULL, therefore ",
                    "not calling seqlevels.."))
       }
   }
   
   if (!is.null(chromosomes)) {
       if (verbose){
           message("*** chromosomes is NOT NULL, will subset..")
       }
       for(i in 1:length(x_gr_list)){
           if (verbose){
               message("*** subset ",i," to only contain: ", chromosomes)
           }
           x_gr_list[[i]] <- subset(x_gr_list[[i]], seqnames %in% chromosomes)
       }
   }
   
   if (!is.null(contexts)){
       if (verbose){
           message("*** contexts is NOT NULL, will subset..")
       }
       for(i in 1:length(x_gr_list)){
           if (verbose){
               message(paste0("*** subset ",i," to only contain: ", contexts))
           }
           x_gr_list[[i]] <- subset(x_gr_list[[i]], context %in% contexts)
       }
   }

   if (!is.null(seqs)){
       if (verbose){
           message("*** seqs is NOT NULL, will subset..")
       }
       for(i in 1:length(x_gr_list)){
           if (verbose){
               message(paste0("*** subset ",i," to only contain: ", seqs))
           }
           x_gr_list[[i]] <- subset(x_gr_list[[i]], seq %in% seqs)
       }
   }

   if (!is.null(sample.id)) {
       names(x_gr_list) <- sample.id
   }
   
   return(x_gr_list)
}
