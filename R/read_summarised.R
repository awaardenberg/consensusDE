#' Generate Summarised Read File for DE analyses
#
#' @description This function will create a summarised file, decribing reads
#' from RNA-seq experiments that overlap a set of transcript features.
#' Transcript features can be described as a gtf formatted table that is
#' imported, or using a txdb. This is designed to be straightforward and
#' minimised parameter for first pass batch RNA-seq analyses.
#'
#' @param sample.table A data.frame describing samples. For paired mode it must
#' contain 3 columns, with the names "file", "group" and "pairs". The filename
#' is the name in the directory supplied with the "bam.dir" parameter below.
#' This is not required if an existing summarised file is provided. Default=NULL
#' @param bam.dir Full path to location of bam files listed in the "file" column
#'  in the sample.table provided above. This is not required if an existing
#'  summarised file is provided. Default=NULL
#' @param gtf Full path to a gtf file describing the transcript coordinates to
#' map the RNA-seq reads to. GTF file is not required if providing a
#' pre-computed summarised experiment file previously generated using
#' read.summarised() OR a tx.db object (below). Default = NULL
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
#' @param summarised Full path to a summarised experiment file. If
#' read.summarised() has already been performed, the output summarised file,
#' saved in "/output.log/se.R" can be used as the input (e.g. if filtering is to
#'  be done). Default=NULL
#' @param output.log Full path to directory for output of log files and saved
#' summarised experiment generated.
#' @param filter Perform filtering of low count and missing data from the
#' summarised experiment file? This uses default filtering via "filterByExpr".
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
#' ## Extract summarised following example in the vignette
#' ## load annotation tx.db for mapping reads
#' library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
#'
#' ## build a design table that lists the files and their grouping
#' file.list <- list.files(system.file("extdata", package="GenomicAlignments"),
#'                         recursive = TRUE,
#'                         pattern = "*bam$",
#'                         full = TRUE)
#' \dontrun{
#' ## create a sample table to be used with read.summarised()
#' ## must be comprised of a minimum of two columns, named "file" and "group",
#' ## with one additional column: "pairing" if the data is paired
#' sample.table <- data.frame("file" = basename(file.list),
#'                            "group" = c("treat", "untreat"))
#'
#' # extract the path to the bam directory - where to search for files listed in
#' ## "sample.table"
#' bam.dir <- as.character(gsub(basename(file.list)[1], "", file.list[1]))
#' summarised.dm3 <- read.summarised(sample.table = sample.table,
#'                                   bam.dir = bam.dir,
#'                                   tx.db = TxDb.Dmelanogaster.UCSC.dm3.ensGene,
#'                                   read.format = "paired")
#' }
#' @return A summarised experiment
#'
#' @export read.summarised
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom Rsamtools BamFileList
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy
#' @importFrom GenomicAlignments summarizeOverlaps
#' @importFrom BiocParallel register MulticoreParam
#' @importFrom edgeR filterByExpr
# @importFrom S4Vectors DataFrame

read.summarised <- function(sample.table = NULL,
                            bam.dir = NULL,
                            gtf = NULL,
                            tx.db = NULL,
                            mapping.mode = "Union",
                            read.format = NULL,
                            ignore.strand = FALSE,
                            fragments = TRUE,
                            summarised = NULL,
                            output.log = NULL,
                            filter = FALSE,
                            BamFileList.yiedsize = NA_integer_,
                            n.cores = 1,
                            verbose = FALSE){

####///---- check inputs ----\\\###
if(is.null(summarised) & (is.null(read.format)))
  stop("read.format must be specified as either \"paired\", or \"single\" if
       a summarised file has not been generated.")
if(is.null(summarised) & (!is.null(read.format))){
  if(read.format != "paired" & read.format != "single"){
    stop("read.format must be specified as either \"paired\", or \"single\" if
       a summarised file has not been generated.")
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
if(is.null(bam.dir) & is.null(summarised))
  stop("A directory to Bam files must be provided when a summarised file is not
       provided.")
if(is.null(gtf) & is.null(tx.db) & is.null(summarised))
  stop("No summarised file is provided and missing a GTF or txdb file. Provide
        full path to gtf file or a txdb for mapping bam file reads")
if(is.null(sample.table) & is.null(summarised))
  stop("No summarised file is provided and missing sample table! Input a sample
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

# create a new summarised object and save
if(is.null(summarised)){
  if(verbose){
    cat(paste("# NO summarised experiment provided", "\n", sep=""))
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
    cat(paste("# Building summarised experiment", "\n", sep=""))
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
      cat(paste("# Summarised experiment saved to:", "\n",
                "# ", output.log, "se.R", sep=""))
    }
  }
se.out <- se
}

# if a .se file [or] path to se.R is provided as an input
if(!is.null(summarised)){
  # either a summarised file already in R, OR read the file in
  if(is.character(summarised)==TRUE){
    # mask any existing environment variables
    attach(summarised, name = "summarised")
    if(verbose){
      cat(paste("# Summarised experiment has been loaded from:", "\n",
                "# ", summarised, "\n", sep=""))
    }
    if(!exists("se"))
      stop(paste("Summarised file provided has not been generated", "\n",
                 "with consensusDE. Please produce a summarised", "\n",
                 "experiment with consensusDE using read.summarised()", "\n",
                 sep=""))
    se.out <- se
    # now detach the file to clean-up the workspace
    detach(summarised)
    # check the file format
    if(se.out@class != "RangedSummarizedExperiment")
      stop(paste("Summarised file provided is not a RangedSummarizedExperiment,
                  Please produce a RangedSummarizedExperiment.", "\n",
                  sep=""))
  }
  # if already a summarised file
  if(is.character(summarised)==FALSE && summarised@class ==
     "RangedSummarizedExperiment"){
    se.out <- summarised
    if(verbose){
      cat(paste("# Summarised experiment provided is as follows:", "\n",
                sep=""))
      print(se.out)
    }
  }

  sample.table <- data.frame(colData(se.out))
  # check the format of the table:
  # check column names exist for sample.table
  if("file" %in% colnames(sample.table) == FALSE){
    warning("The summarised experiment provided does not include a \"file\"
            column. This will create errors when running the DE analysis. Update
            the summarised experiment with the experimental details before
            processing to DE analysis")
  }
  if("group" %in% colnames(sample.table) == FALSE){
    warning("The summarised experiment provided does not include a \"group\"
             column. This will create errors when running the DE analysis.
             Update the summarised experiment with the experimental details
             before processing to DE analysis")
  }

  if(!is.null(output.log)){
    se <- se.out
    save(se, file=paste(output.log, "se.R", sep=""))
    if(verbose){
      cat(paste("# Summarised experiment saved to:", "\n",
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
    cat("# SummarisedFile ready for further analysis\n")
  }
return(se.out)
}
