#' Generate summarized Read File for DE analyses
#
#' @description This function will create a summarized file, decribing reads
#' from RNA-seq experiments that overlap a set of transcript features.
#' Transcript features can be described as a gtf formatted table that is
#' imported, or using a txdb. This is designed to be straightforward and
#' with minimised parameters for first pass batch RNA-seq analyses.
#'
#' @param sample_table A data.frame describing samples. For paired mode it must
#' contain 3 columns, with the names "file", "group" and "pairs". The filename
#' is the name in the directory supplied with the "bam_dir" parameter below.
#' This is not required if an existing summarized file is provided. Default=NULL
#' @param bam_dir Full path to location of bam files listed in the "file" column
#'  in the sample_table provided above. This is not required if an existing
#'  summarized file is provided. Default=NULL
#' @param gtf Full path to a gtf file describing the transcript coordinates to
#' map the RNA-seq reads to. GTF file is not required if providing a
#' pre-computed summarized experiment file previously generated using
#' buildSummarized() OR a tx_db object (below). Default = NULL
#' @param tx_db An R txdb object. E.g. TxDb.Dmelanogaster.UCSC.dm3.ensGene.
#' Default = NULL
#' @param mapping_mode Options are "Union", "IntersectionStrict" and
#' "IntersectionNotEmpty". see "mode" in ?summarizeOverlaps for explanation.
#' Default = "Union"
#' @param read_format Are the reads from single-end or paired-end data? Option
#' are "paired" or "single". An option must be selected. Default = NULL
#' @param ignore_strand Ignore strand when mapping reads? see "ignore_strand" in
#' ?summarizeOverlaps for explanation. Default=FALSE
#' @param fragments When mapping_mode="paired", include reads from pairs that do
#'  not map with their corresponding pair? see "fragments" in ?summarizeOverlaps
#'   for explanation. Default = TRUE
#' @param summarized Full path to a summarized experiment file. If
#' buildSummarized() has already been performed, the output summarized file,
#' saved in "/output_log/se.R" can be used as the input (e.g. if filtering is to
#'  be done). Default = NULL
#' @param output_log Full path to directory for output of log files and saved
#' summarized experiment generated.
#' @param filter Perform filtering of low count and missing data from the
#' summarized experiment file? This uses default filtering via "filterByExpr".
#' See ?filterByExpr for further information. Default=FALSE
#' @param BamFileList_yiedsize If running into memory problems. Set the number
#' of lines to an integer value. See "yieldSize" description in ?BamFileList for
#'  an explanation.
#' @param n_cores Number of cores to utilise for reading in Bam files. Use with
#' caution as can create memory issues if BamFileList_yiedsize is not
#' parameterised. Default = 1
#' @param force_build If the sample_table contains less than two replicates per
#' group, force a summarizedExperiment object to be built. Otherwise 
#' buildSummarized will halt. Default = FALSE.
#' @param verbose Verbosity ON/OFF. Default=FALSE
#'
#' @examples
#' ## Extract summarized following example in the vignette
#' ## The example below will return a summarized experiment
#' ## tx_db is obtained from TxDb.Dmelanogaster.UCSC.dm3.ensGene library
#' library(TxDb.Dmelanogaster.UCSC.dm3.ensGene)
#' ## bam files are obtained from the GenomicAlignments package
#' 
#' ## 1. Build a sample table that lists files and groupings
#' ## - obtain list of files
#' file_list <- list.files(system.file("extdata", package="GenomicAlignments"),
#'                         recursive = TRUE,
#'                         pattern = "*bam$",
#'                         full = TRUE)
#' bam_dir <- as.character(gsub(basename(file_list)[1], "", file_list[1]))
#'
#' ## - create a sample table to be used with buildSummarized() below
#' ## must be comprised of a minimum of two columns, named "file" and "group",
#' sample_table <- data.frame("file" = basename(file_list),
#'                            "group" = c("treat", "untreat"))
#'
# ## 2. Build summarized experiment, from sample table and dm3 txdb
#' summarized_dm3 <- buildSummarized(sample_table = sample_table,
#'                                   bam_dir = bam_dir,
#'                                   tx_db = TxDb.Dmelanogaster.UCSC.dm3.ensGene,
#'                                   read_format = "paired",
#'                                   force_build = TRUE)
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
#' @import Biostrings

buildSummarized <- function(sample_table = NULL,
                            bam_dir = NULL,
                            gtf = NULL,
                            tx_db = NULL,
                            mapping_mode = "Union",
                            read_format = NULL,
                            ignore_strand = FALSE,
                            fragments = TRUE,
                            summarized = NULL,
                            output_log = NULL,
                            filter = FALSE,
                            BamFileList_yiedsize = NA_integer_,
                            n_cores = 1,
                            force_build = FALSE,
                            verbose = FALSE){

####///---- check inputs ----\\\###
if(is.null(summarized) & (is.null(read_format)))
  stop("read_format must be specified as either \"paired\", or \"single\" if
       a summarized file has not been generated.")
if(is.null(summarized) & (!is.null(read_format))){
  if(read_format != "paired" & read_format != "single"){
    stop("read_format must be specified as either \"paired\", or \"single\" if
       a summarized file has not been generated.")
  }
  # define mode for summarizeOverlaps
  if(read_format == "paired")
    singleEnd_paired <- FALSE
  if(read_format == "single")
    singleEnd_paired <- TRUE
}

if(!is.null(tx_db) & !is.null(gtf)){
  warning("Both a tx_db object and path to gtf file have been provided. The path
            to the gtf file will be used in this instance.")
  tx_db <- NULL
}
      
# be careful with more than one worker here: is extremely memory intense!
# check n_cores is integer; BamFileList_yiedsize is NA_integer_ or an integer...
is_wholenumber <-
  function(x,
           tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
if(!is_wholenumber(n_cores))
  stop("n_cores must be an integer")
  # register the number of cores to used:
  register(MulticoreParam(workers=n_cores))
if(is.null(bam_dir) & is.null(summarized))
  stop("A directory to Bam files must be provided when a summarized file is not
       provided.")
if(is.null(gtf) & is.null(tx_db) & is.null(summarized))
  stop("No summarized file is provided and missing a GTF or txdb file. Provide
        full path to gtf file or a txdb for mapping bam file reads")
if(is.null(sample_table) & is.null(summarized))
  stop("No summarized file is provided and missing sample table! Input a sample
        table using sample_table parameter")
if(is.null(output_log))
  warning("No output directory provided. The se file and sample_table will not
          be saved")
if(!(filter == TRUE || filter == FALSE))
  stop("filter can only be filter = TRUE or filter = FALSE\n")

if(filter == TRUE & (is.null(sample_table) & is.null(summarized)))
  stop("Filtering can only be done if a sample_table table has been provided 
       with groups or a previously generated summarized experiment is provided")

####///---- input checks DONE ----\\\###

# create a new summarized object and save
if(is.null(summarized)){
  if(verbose){
    message("# NO summarized experiment provided")
  }
  # check column names exist for sample_table
  if("file" %in% colnames(sample_table) == FALSE){
    stop("A sample_table must be supplied with a column labelled \"file\"")
  }
  if("group" %in% colnames(sample_table) == FALSE){
    stop("A sample_table must be supplied with a column labelled \"group\"")
  }
  # check that there are minimum of two replicates in groups...
  if(force_build == FALSE && min(summary(sample_table$group)) < 2)
    stop("The sample_table provided contains groups with less than two 
         replicates")
  if(force_build == TRUE && min(summary(sample_table$group)) < 2)
    warning("The sample_table provided contains groups with less than two 
         replicates. You have selected to continue with force_build = TRUE")
  # check paired options are not matching the groups... i.e. replicated
  if("pairs" %in% colnames(sample_table) == TRUE){
    check_pairs <- paste(sample_table$group, sample_table$pairs, sep="_")
    if(length(unique(check_pairs)) != nrow(sample_table)){
      stop("pairs column in sample_table contains pairings from same group.
           Technical replication is not supported.")
    }
    if(min(summary(sample_table$pairs)) < 2)
      stop("The sample_table provided contains pairs with less than two 
           samples.")
  }
  
  bam_files <- list.files(bam_dir, full.names = FALSE, pattern = ".bam")
  # must be exact match before proceeding
  if(length(intersect(bam_files, sample_table$file)) != nrow(sample_table))
    warning(paste("There is not a matching number of files from:\n",
                  sample_table, "\n",
               "in the directory provided:", bam_dir, "\n",
               "Action: Check file names match or correct sample_table
                file/format provided.", sep=" "))
  # if a txdb is not provided, but a gtf object is:
  if(is.null(tx_db) & !is.null(gtf)){
      txdb <- makeTxDbFromGFF(gtf, format="gtf", circ_seqs=character())
  }
  # if a gtf is not provided, but a txdb object is:
  if(!is.null(tx_db) & is.null(gtf)){
      txdb <- tx_db
  }
  # extract exon coordinates, by default.
  ebg <- exonsBy(txdb, by="gene")
  # Read bam files in
  bam_files <- paste(bam_dir, sample_table$file, sep="")
  bamfiles <- BamFileList(bam_files, yieldSize = BamFileList_yiedsize)
  if(verbose){
    message("# Building summarized experiment")
  }
  # this is fine for RNA-seq, as it may be possible that only one of the pairs
  # aligns to the feature of interest due to transcript annotaiton. There is an
  # option to turn off.
  se <- summarizeOverlaps(features = ebg,
                          reads = bamfiles,
                          mode = mapping_mode,
                          singleEnd = singleEnd_paired,
                          ignore.strand = ignore_strand,
                          fragments = fragments)

  # ensure SE is labelled (important for model fits later)
  colData(se) <- S4Vectors::DataFrame(sample_table)
  colnames(se) <- sample_table$file

  # name will be the filename, not including dir.
  # save se for future use
  if(!is.null(output_log)){
    save(se, file=paste(output_log, "se.R", sep=""))
    if(verbose){
      message("# summarized experiment saved to:", 
              paste(output_log, "se.R", sep=""))
    }
  }
se_out <- se
}

# if a .se file [or] path to se.R is provided as an input
if(!is.null(summarized)){
  # either a summarized file already in R, OR read the file in
  if(is.character(summarized)==TRUE){
    # mask any existing environment variables
    attach(summarized, name = "summarized")
    if(verbose){
      message("# summarized experiment has been loaded from:", summarized)
    }
    if(!exists("se"))
      stop(paste("summarized file provided has not been generated", "\n",
                 "with consensusDE. Please produce a summarized", "\n",
                 "experiment with consensusDE using buildSummarized()", "\n",
                 sep=""))
    se_out <- se
    # now detach the file to clean-up the workspace
    detach(summarized)
    # check the file format
    if(se_out@class != "RangedSummarizedExperiment")
      stop(paste("summarized file provided is not a RangedSummarizedExperiment,
                  Please produce a RangedSummarizedExperiment.", "\n",
                  sep=""))
  }
  # if already a summarized file
  if(is.character(summarized)==FALSE && summarized@class ==
     "RangedSummarizedExperiment"){
    se_out <- summarized
    if(verbose){
      message("# summarized experiment provided is as follows:")
      message(se_out)
    }
  }

  sample_table <- data.frame(colData(se_out))
  # check the format of the table:
  # check column names exist for sample_table
  if("file" %in% colnames(sample_table) == FALSE){
    warning("The summarized experiment provided does not include a \"file\"
            column. This will create errors when running the DE analysis. Update
            the summarized experiment with the experimental details before
            processing to DE analysis")
  }
  if("group" %in% colnames(sample_table) == FALSE){
    warning("The summarized experiment provided does not include a \"group\"
             column. This will create errors when running the DE analysis.
             Update the summarized experiment with the experimental details
             before processing to DE analysis")
  }

  if(!is.null(output_log)){
    se <- se_out
    save(se, file=paste(output_log, "se.R", sep=""))
    if(verbose){
      message("# summarized experiment saved to:", 
              paste(output_log, "se.R", sep=""))
    }
  }
}

# report table and number of bam files (either from input, or from se file)
if(verbose){
  message(sample_table)
  message("#", nrow(sample_table), "bam file[s] selected")
}

# option to write sample_table to log dir
if(!is.null(output_log)){
  utils::write.table(sample_table, file=paste(output_log,
                                              "sample_table_input.tsv", sep=""),
              sep="\t", row.names=FALSE, quote=FALSE)
  if(verbose){
    message("# sample table saved to:", 
            paste(output_log, "sample_table_input.tsv", sep=""))
  }
}

# option to filter data
if(filter == TRUE){
  # filter low count data based on group assignments
  keep <- filterByExpr(assays(se_out)$counts, group=colData(se_out)$group)
  se_out <- se_out[rownames(se_out)[keep] ,]
}
  if(verbose){
    message("# summarizedFile ready for further analysis")
  }
return(se_out)
}
