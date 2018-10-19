#' Generate summarized Read File for DE analyses
#
#' @description This function will create a summarized file, decribing reads
#' from RNA-seq experiments that overlap a set of transcript features.
#' Transcript features can be described as a gtf formatted table that is
#' imported, or using a txdb. This is designed to be straightforward and
#' minimised parameter for first pass batch RNA-seq analyses.
#'
#' @param sample.table A data.frame describing samples. For paired mode it must
#' contain 3 columns, with the names "file", "group" and "pairs". The filename
#' is the name in the directory supplied with the "bam.dir" parameter below.
#' This is not required if an existing summarized file is provided. Default=NULL
#' @param bam.dir Full path to location of bam files listed in the "file" column
#'  in the sample.table provided above. This is not required if an existing
#'  summarized file is provided. Default=NULL
#' @param gtf Full path to a gtf file describing the transcript coordinates to
#' map the RNA-seq reads to. GTF file is not required if providing a
#' pre-computed summarized experiment file previously generated using
#' buildSummarized() OR a tx.db object (below). Default = NULL
#' @param tx.db An R txdb object. E.g. TxDb.Dmelanogaster.UCSC.dm3.ensGene.
#' Default=NULL
#' @param mapping.mode Options are "Union", "IntersectionStrict" and
#' "IntersectionNotEmpty". see "mode" in ?summarizeOverlaps for explanation.
#' Default = "Union"
#' @param read.format Are the reads from single-end or paired-end data? Option
#' are "paired" or "single". An option must be selected. Default = NULL
#' @param ignore.strand Ignore strand when mapping reads? see "ignore.strand" in
#' ?summarizeOverlaps for explanation. Default=FALSE
#' @param fragments When mapping.mode="paired", include reads from pairs that do
#'  not map with their corresponding pair? see "fragments" in ?summarizeOverlaps
#'   for explanation. Default=TRUE
#' @param summarized Full path to a summarized experiment file. If
#' buildSummarized() has already been performed, the output summarized file,
#' saved in "/output.log/se.R" can be used as the input (e.g. if filtering is to
#'  be done). Default=NULL
#' @param output.log Full path to directory for output of log files and saved
#' summarized experiment generated.
#' @param filter Perform filtering of low count and missing data from the
#' summarized experiment file? This uses default filtering via "filterByExpr".
#' See ?filterByExpr for further information. Default=FALSE
#' @param BamFileList.yiedsize If running into memory problems. Set the number
#' of lines to an integer value. See "yieldSize" description in ?BamFileList for
#'  an explanation.
#' @param n.cores Number of cores to utilise for reading in Bam files. Use with
#' caution as can create memory issues if BamFileList.yiedsize is not
#' parameterised. Default=1
#' @param verbose Verbosity ON/OFF. Default=FALSE
#'
#' @examples
#' ## Extract summarized following example in the vignette
#' ## load annotation tx.db for mapping reads
#' cat("The examples below will return a summarized experiment")
#' 
#' library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
#'
#' ## build a design table that lists the files and their grouping
#' file.list <- list.files(system.file("extdata", package="GenomicAlignments"),
#'                         recursive = TRUE,
#'                         pattern = "*bam$",
#'                         full = TRUE)
#'
#' ## create a sample table to be used with buildSummarized()
#' ## must be comprised of a minimum of two columns, named "file" and "group",
#' ## with one additional column: "pairing" if the data is paired
#' sample.table <- data.frame("file" = basename(file.list),
#'                            "group" = c("treat", "untreat"))
#'
#' # extract the path to the bam directory - where to search for files listed in
#' ## "sample.table"
#' bam.dir <- as.character(gsub(basename(file.list)[1], "", file.list[1]))
#' summarized.dm3 <- buildSummarized(sample.table = sample.table,
#'                                   bam.dir = bam.dir,
#'                                   tx.db = TxDb.Dmelanogaster.UCSC.dm3.ensGene,
#'                                   read.format = "paired")
#' 
#' @return A summarized experiment
#'
#' @export buildSummarized
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom Rsamtools BamFileList
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom BiocParallel register MulticoreParam
#' @importFrom edgeR filterByExpr
#' @import TxDb.Dmelanogaster.UCSC.dm3.ensGene

buildSummarized <- function(sample.table = NULL,
                            bam.dir = NULL,
                            gtf = NULL,
                            tx.db = NULL,
                            mapping.mode = "Union",
                            read.format = NULL,
                            ignore.strand = FALSE,
                            fragments = TRUE,
                            summarized = NULL,
                            output.log = NULL,
                            filter = FALSE,
                            BamFileList.yiedsize = NA_integer_,
                            n.cores = 1,
                            verbose = FALSE){

####///---- check inputs ----\\\###
if(is.null(summarized) & (is.null(read.format)))
  stop("read.format must be specified as either \"paired\", or \"single\" if
       a summarized file has not been generated.")
if(is.null(summarized) & (!is.null(read.format))){
  if(read.format != "paired" & read.format != "single"){
    stop("read.format must be specified as either \"paired\", or \"single\" if
       a summarized file has not been generated.")
  }
  # define mode for summarizeOverlaps
  if(read.format == "paired")
    singleEnd.paired <- FALSE
  if(read.format == "single")
    singleEnd.paired <- TRUE
}

if(!is.null(tx.db) & !is.null(gtf))
  warning("Both a tx.db object and path to gtf file have been provided. The path
            to the gtf file will be used in this instance.")
# be careful with more than one worker here: is extremely memory intense!
# check n.cores is integer; BamFileList.yiedsize is NA_integer_ or an integer...
is.wholenumber <-
  function(x,
           tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
if(!is.wholenumber(n.cores))
  stop("n.cores must be an integer")
  # register the number of cores to used:
  register(MulticoreParam(workers=n.cores))
if(is.null(bam.dir) & is.null(summarized))
  stop("A directory to Bam files must be provided when a summarized file is not
       provided.")
if(is.null(gtf) & is.null(tx.db) & is.null(summarized))
  stop("No summarized file is provided and missing a GTF or txdb file. Provide
        full path to gtf file or a txdb for mapping bam file reads")
if(is.null(sample.table) & is.null(summarized))
  stop("No summarized file is provided and missing sample table! Input a sample
        table using sample.table parameter")
if(is.null(output.log))
  warning("No output directory provided. The se file and sample.table will not
          be saved")
if(!(filter == TRUE || filter == FALSE))
  stop("filter can only be filter=TRUE or filter=FALSE\n")
if(filter != TRUE & is.null(sample.table))
  stop("Filtering can only be done if a sample.table table
       has been provided with groups\n")

####///---- input checks DONE ----\\\###

# create a new summarized object and save
if(is.null(summarized)){
  if(verbose){
    cat(paste("# NO summarized experiment provided", "\n", sep=""))
  }
  # check column names exist for sample.table
  if("file" %in% colnames(sample.table) == FALSE){
    stop("A sample.table must be supplied with a column labelled \"file\"")
  }
  if("group" %in% colnames(sample.table) == FALSE){
    stop("A sample.table must be supplied with a column labelled \"group\"")
  }
  bam_files <- list.files(bam.dir, full.names = FALSE, pattern = ".bam")
  # must be exact match before proceeding

  if(length(intersect(bam_files, sample.table$file)) != nrow(sample.table))
    warning(paste("There is not a matching number of files from:\n",
                  sample.table, "\n",
               "in the directory provided:", bam.dir, "\n",
               "Action: Check file names match or correct sample.table
                file/format provided.", sep=" "))

    if(!is.null(tx.db) & !is.null(gtf)){
      txdb <- makeTxDbFromGFF(gtf, format="gtf", circ_seqs=character())
  }
  if(!is.null(tx.db) & is.null(gtf)){
      txdb <- tx.db
  }
  # extract exon coordinates, by default.
  ebg <- exonsBy(txdb, by="gene")
  # Read bam files in
  bam_files <- paste(bam.dir, sample.table$file, sep="")
  bamfiles <- BamFileList(bam_files, yieldSize = BamFileList.yiedsize)
  if(verbose){
    cat(paste("# Building summarized experiment", "\n", sep=""))
  }
  # this is fine for RNA-seq, as it may be possible that only one of the pairs
  # aligns to the feature of interest due to transcript annotaiton. There is an
  # option to turn off.
  se <- summarizeOverlaps(features = ebg,
                          reads = bamfiles,
                          mode = mapping.mode,
                          singleEnd = singleEnd.paired,
                          ignore.strand = ignore.strand,
                          fragments = fragments)

  # ensure SE is labelled (important for model fits later)
  colData(se) <- S4Vectors::DataFrame(sample.table)
  colnames(se) <- sample.table$file

  # name will be the filename, not including dir.
  # save se for future use
  if(!is.null(output.log)){
    save(se, file=paste(output.log, "se.R", sep=""))
    if(verbose){
      cat(paste("# summarized experiment saved to:", "\n",
                "# ", output.log, "se.R", sep=""))
    }
  }
se.out <- se
}

# if a .se file [or] path to se.R is provided as an input
if(!is.null(summarized)){
  # either a summarized file already in R, OR read the file in
  if(is.character(summarized)==TRUE){
    # mask any existing environment variables
    attach(summarized, name = "summarized")
    if(verbose){
      cat(paste("# summarized experiment has been loaded from:", "\n",
                "# ", summarized, "\n", sep=""))
    }
    if(!exists("se"))
      stop(paste("summarized file provided has not been generated", "\n",
                 "with consensusDE. Please produce a summarized", "\n",
                 "experiment with consensusDE using buildSummarized()", "\n",
                 sep=""))
    se.out <- se
    # now detach the file to clean-up the workspace
    detach(summarized)
    # check the file format
    if(se.out@class != "RangedSummarizedExperiment")
      stop(paste("summarized file provided is not a RangedSummarizedExperiment,
                  Please produce a RangedSummarizedExperiment.", "\n",
                  sep=""))
  }
  # if already a summarized file
  if(is.character(summarized)==FALSE && summarized@class ==
     "RangedSummarizedExperiment"){
    se.out <- summarized
    if(verbose){
      cat(paste("# summarized experiment provided is as follows:", "\n",
                sep=""))
      print(se.out)
    }
  }

  sample.table <- data.frame(colData(se.out))
  # check the format of the table:
  # check column names exist for sample.table
  if("file" %in% colnames(sample.table) == FALSE){
    warning("The summarized experiment provided does not include a \"file\"
            column. This will create errors when running the DE analysis. Update
            the summarized experiment with the experimental details before
            processing to DE analysis")
  }
  if("group" %in% colnames(sample.table) == FALSE){
    warning("The summarized experiment provided does not include a \"group\"
             column. This will create errors when running the DE analysis.
             Update the summarized experiment with the experimental details
             before processing to DE analysis")
  }

  if(!is.null(output.log)){
    se <- se.out
    save(se, file=paste(output.log, "se.R", sep=""))
    if(verbose){
      cat(paste("# summarized experiment saved to:", "\n",
                "# ", output.log, "se.R", sep=""))
    }
  }
}

# report table and number of bam files (either from input, or from se file)
if(verbose){
  print(sample.table)
  cat(paste("\n", "#", nrow(sample.table), "bam file[s] selected\n\n", sep=" "))
}

# option to write sample.table to log dir
if(!is.null(output.log)){
  utils::write.table(sample.table, file=paste(output.log,
                                              "sample.table_input.tsv", sep=""),
              sep="\t", row.names=FALSE, quote=FALSE)
  if(verbose){
    cat(paste("# sample table saved to:", "\n",
              "# ", output.log, "sample.table_input.tsv", "\n", sep=""))
  }
}

# option to filter data
if(filter == TRUE){
  # filter low count data based on group assignments
  keep <- filterByExpr(assays(se.out)$counts, group=sample.table$group)
  se.out <- se.out[rownames(se.out)[keep] ,]
}
  if(verbose){
    cat("# summarizedFile ready for further analysis\n")
  }
return(se.out)
}
