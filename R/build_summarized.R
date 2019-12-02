#' Generate summarized Read File for DE analyses
#
#' @description This function will create a summarized experiment, decribing 
#' reads from RNA-seq experiments that overlap a set of transcript features.
#' Transcript features can be described as a gtf formatted table that is
#' imported, or using a txdb. The summarized experiment can be build directly
#' from bam files or by reading in counts in htseq format. This is designed to 
#' be straightforward and with minimised parameters for batch style RNA-seq 
#' analyses.
#'
#' @param sample_table A data.frame describing samples. For paired mode it must
#' at least 2 columns, "file", "group", and option additional columns, "pairs"
#' and "tech_replicate" for describing sample pairing and instances of technical
#' replicates. The filename "file" must correspong to the name of the file in 
#' the directory supplied with the "bam_dir" parameter below - or ar error will
#' be reported and buildSummarized will halt. This is not required if an 
#' existing summarized file is provided. Default = NULL
#' @param bam_dir Full path to location of bam files listed in the "file" column
#'  in the sample_table provided above. This is not required if an existing
#'  summarized file is provided. Default = NULL
#' @param htseq_dir Full path to location of htseq files listed in the "file" 
#' column in the sample_table described above. This is not required if an 
#' existing summarized file is provided. Files must end in ".txt". Default = 
#' NULL
#' @param gtf Full path to a gtf file describing the transcript coordinates to
#' map the RNA-seq reads to. GTF file is not required if providing a
#' pre-computed summarized experiment file previously generated using
#' buildSummarized() OR a tx_db object (below). Default = NULL
#' @param tx_db An R txdb object. E.g. TxDb.Dmelanogaster.UCSC.dm3.ensGene.
#' Default = NULL
#' @param technical_reps Are there technical replicates to merge counts? I.e. 
#' are there multiple technical replicates run accross multiple lanes/sequencing
#' runs. If "TRUE", unique sample names should be provided in a "tech_replicate"
#' column of the "sample_table" for identification. Options are "TRUE" or 
#' "FALSE". Default = "FALSE"
#' @param map_reads Which features to count reads by. Options are "transcript", 
#'  "exon" or "cds". This will invoke transcriptsBy(), exonsBy() or cdsBy() 
#'  respectively. Default = "transcript"
#' @param mapping_mode Options are "Union", "IntersectionStrict" and
#' "IntersectionNotEmpty". see "mode" in ?summarizeOverlaps for explanation.
#' Default = "Union"
#' @param read_format Are the reads from single-end or paired-end data? Option
#' are "paired" or "single". An option must be selected if htseq_dir is NULL and
#' read are summarized from BAM files. Default = NULL
#' @param strand_mode indicates how the reads are stranded see ?strandMode in
#' Genomic Alignments for explanation. Default = 1
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
#' See ?filterByExpr for further information. Default = FALSE
#' @param BamFileList_yieldsize If running into memory problems. Set the number
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
#' @importFrom SummarizedExperiment colData colData<- SummarizedExperiment rowRanges<-
#' @importFrom S4Vectors metadata metadata<- SimpleList
#' @importFrom Rsamtools BamFileList
#' @importFrom GenomicFeatures makeTxDbFromGFF exonsBy transcriptsBy cdsBy genes
#' @importFrom GenomicAlignments summarizeOverlaps invertStrand
#' @importFrom BiocParallel register MulticoreParam
#' @importFrom edgeR filterByExpr
#' @import TxDb.Dmelanogaster.UCSC.dm3.ensGene
#' @import Biostrings
#' @importFrom utils read.table
#' @importFrom ensembldb transcripts

buildSummarized <- function(sample_table = NULL,
                            bam_dir = NULL,
                            htseq_dir = NULL,
                            gtf = NULL,
                            tx_db = NULL,
                            technical_reps = FALSE,
                            map_reads = "transcript",
                            mapping_mode = "Union",
                            read_format = NULL,
                            strand_mode = 1,
                            fragments = FALSE,
                            summarized = NULL,
                            output_log = NULL,
                            filter = FALSE,
                            BamFileList_yieldsize = NA_integer_,
                            n_cores = 1,
                            force_build = FALSE,
                            verbose = FALSE){
  
  ####///---- check inputs ----\\\###
  if(is.null(summarized) & is.null(htseq_dir) & (is.null(read_format)))
    stop("read_format must be specified as either \"paired\", or \"single\" if
       a summarized file or htseq_dir has not been generated .")
  if(is.null(summarized)){
    if(!is.null(bam_dir)){
      if(is.null(read_format)){
        stop("read_format must be specified as either \"paired\", or \"single\" if
           a summarized file has not been generated and read counting is from
           bam files.")
      }
      if((read_format == "paired" | read_format == "single") == FALSE){
        stop("read_format must be specified as either \"paired\", or \"single\" if
         a summarized file has not been generated and read counting is from
           bam files.")
      }
      # define modes for summarizeOverlaps
      ## paired end vs single end
      if(read_format == "paired"){
        singleEnd_paired <- FALSE
      }
      if(read_format == "single"){
        singleEnd_paired <- TRUE
      }
      ## strandMode
      if(strand_mode == 0){
        ignore_strand = TRUE
        preprocess.reads = NULL
      }
      if(strand_mode == 1){
        ignore_strand = FALSE
        preprocess.reads = NULL
      }
      if(strand_mode == 2){
        ignore_strand = FALSE
        preprocess.reads = invertStrand
      }
    }
  }
  
  if(!is.null(tx_db) & !is.null(gtf)){
    warning("Both a tx_db object and path to gtf file have been provided. The path
            to the gtf file will be used in this instance.")
    tx_db <- NULL
  }
  if((technical_reps != TRUE) & (technical_reps != FALSE))
    stop("technical_reps can only be either \"TRUE\" or \"FALSE\". Please specify")
  
  # be careful with more than one worker here: is extremely memory intense!
  # check n_cores is integer; BamFileList_yiedsize is NA_integer_ or an integer...
  is_wholenumber <-
    function(x,
             tol = .Machine$double.eps^0.5) abs(x - round(x)) < tol
  if(!is_wholenumber(n_cores))
    stop("n_cores must be an integer")
  # register the number of cores to used:
  register(MulticoreParam(workers=n_cores))
  if(is.null(bam_dir) & is.null(htseq_dir) & is.null(summarized))
    stop("A directory to .bam or htseq text files must be provided when a 
        summarized file is not provided.")
  if(!is.null(bam_dir) & !is.null(htseq_dir))
    stop("EITHER a directory to .bam or htseq text files can be provided - NOT 
        BOTH. Please provide a path to only one and rerun.")
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
  
  if(!(map_reads == "transcript" || map_reads == "exon" || map_reads ==  "cds"))
    stop("map_reads must be specified as either \"transcript\", \"exon\", 
       or \"cds\"")
  
  ####///---- input checks DONE ----\\\###
  # intialise stats for metadata
  stats <- c()
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
    if(technical_reps == TRUE & ("tech_replicate" %in% colnames(sample_table) == FALSE)){
      stop("For technical_reps data, sample_table must be supplied with a column 
          labelled \"tech_replicate\"")
    }
    # ensure sample table groups are refactored
    sample_table$group <- as.factor(as.character(sample_table$group))
    # check that there are minimum of two replicates in groups...
    if(force_build == FALSE && min(summary(sample_table$group)) < 2)
      stop("The sample_table provided contains groups with less than two 
         replicates")
    if(force_build == TRUE && min(summary(sample_table$group)) < 2)
      warning("The sample_table provided contains groups with less than two 
         replicates. You have selected to continue with force_build = TRUE")
    # check paired options are not matching the groups... i.e. replicated
    
    if("pairs" %in% colnames(sample_table) == TRUE){
      # if technical replicates are to be merged
      if("tech_replicate"  %in% colnames(sample_table) == TRUE){
        check_techs <- paste(sample_table$group, sample_table$pairs, sample_table$tech_replicate, sep="_")
        check_pairs <- paste(sample_table$group, sample_table$pairs, sep="_")
        if(length(unique(check_pairs)) != length(unique(check_techs))){
          #if(length(unique(check_pairs)) != nrow(sample_table)){
          stop("pairs column in sample_table contains pairings from same group OR
             there is a problem with labelling of \"tech_replicate\" column.")
        }
      }
      # if no technical replicates are to be merged
      if("tech_replicate"  %in% colnames(sample_table) == FALSE){
        check_pairs <- paste(sample_table$group, sample_table$pairs, sep="_")
        if(length(unique(check_pairs)) != nrow(sample_table)){
          stop("pairs column in sample_table contains pairings from same group.
             Technical replication is not supported. If technical replicates are
             present set \"technical_reps\" to \"TRUE\" and ensure a 
             \"tech_replicate\" column is included in your \"sample_table\"")
        }
      }
      if(as.numeric(min(summary(unique(sample_table$pairs)))) < 2)
        stop("sample_table contains pairs with less than two samples")
    }
    
    if(!is.null(bam_dir)){
      input_files <- list.files(bam_dir, full.names = FALSE, pattern = ".bam")
    }
    if(!is.null(htseq_dir)){
      input_files <- list.files(htseq_dir, full.names = FALSE, pattern = ".txt")
    }
    # must be exact match before proceeding
    if(length(intersect(input_files, sample_table$file)) != nrow(sample_table))
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
    # optional mapping methods (only valid if reading in bam files?)
    if(!is.null(bam_dir) & map_reads == "transcript"){
      ebg <- transcriptsBy(x = txdb, by = "gene")
    }
    if(!is.null(bam_dir) & map_reads == "exon"){
      ebg <- exonsBy(x = txdb, by = "gene")
    }
    if(!is.null(bam_dir) & map_reads == "cds"){
      ebg <- cdsBy(x = txdb, by = "gene")
    }
    if(!is.null(htseq_dir)){
      if(verbose){
        message("HTseq counts selected. Txdb will be summarized at exon level.")
      }
      ebg <- exonsBy(x = txdb, by = "gene")
    }
    if(verbose){
      message("# Building summarized experiment")
    }
    if(!is.null(htseq_dir)){
      # list files
      htseq_files <- paste(htseq_dir, sample_table$file, sep="")
      # read files into tables
      htseq_files <- lapply(seq_along(htseq_files), function(i)
        read.table(htseq_files[i], 
                   col.names = c("ID", 
                                 as.character(sample_table$file[i]))))
      #  merge together...
      counts <- Reduce(function(...) merge(..., all=TRUE, by="ID"), htseq_files) 
      # some stats to keep
      stats <- data.frame(counts[grep("__", counts$ID),])
      #remove from counts
      counts <- counts[!as.character(counts$ID) %in% as.character(stats$ID),]
      rownames(counts) <- counts$ID
      counts <- counts[!colnames(counts) %in% c("ID")]
      se <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)))
      rowRanges(se) <- ebg
      #metadata(se)$stats <- stats
    }
    
    if(!is.null(bam_dir)){
      # establish bam files to read in
      bam_files <- paste(bam_dir, sample_table$file, sep="")
      
      if(read_format == "paired"){
        bamfiles <- BamFileList(bam_files, yieldSize = BamFileList_yiedsize, 
                                asMates = TRUE)
      }
       if(read_format == "single"){
         bamfiles <- BamFileList(bam_files, yieldSize = BamFileList_yiedsize)
       }
      
      se <- summarizeOverlaps(features = ebg,
                              reads = bamfiles,
                              mode = mapping_mode,
                              singleEnd = singleEnd_paired,
                              ignore.strand = ignore_strand,
                              fragments = fragments,
                              preprocess.reads = preprocess_reads)
    }
    
    # ensure SE is labelled (important for model fits later)
    colData(se) <- S4Vectors::DataFrame(sample_table)
    colnames(se) <- sample_table$file
    #colnames(se) <- make.names(sample_table$file)
    
    #will rebuild the SE if technical_reps is true
    if(technical_reps == TRUE){
      # perform sample merging here and update sample_table 
      what_to_merge <- unique(data.frame(colData(se))$tech_replicate)
      multiplex_data <- lapply(seq_along(what_to_merge), function(x)
        subset_se(se_in = se, 
                  multiplex_id = what_to_merge[x]))
      # merge together...
      counts <- Reduce(function(...) merge(..., all=TRUE, by="ID"), multiplex_data) 
      rownames(counts) <- counts$ID
      counts <- counts[!colnames(counts) %in% c("ID")]
      se <- SummarizedExperiment(assays=SimpleList(counts=as.matrix(counts)))
      rowRanges(se) <- ebg
      
      #update the sample_table $file lables
      update_sample_table <- sample_table[colnames(sample_table) %in% c("file", "group", "pairs", "tech_replicate")]
      update_sample_table$file <- update_sample_table$tech_replicate
      
      # ensure SE is labelled (important for model fits later)
      colData(se) <- S4Vectors::DataFrame(unique(update_sample_table))
      colnames(se) <- colnames(counts)
      #colnames(se) <- make.names(sample_table$file)
    }
    
    # add metadata to summarized object
    metadata(se)$gene_coords <- genes(txdb)
    metadata(se)$sample_table <- sample_table
    #metadata(se)$consensusDE_parameters <-
    metadata(se)$stats <- stats
    
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
    if(is.null(tx_db) == FALSE){
      metadata(se_out)$gene_coords <- genes(tx_db)
      message("A txdb was provided as input. Meta-data has been updated")
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
  
  # report table and number of samples (either from input, or from se file)
  if(verbose){
    #message(sample_table)
    message("#", nrow(sample_table), "sample[s] present")
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


# this function is for summing over multi-plex samples
# here multi-plexed refers to the same techical replicate samples accross
# multiple lanes
# therefore total read count is the sum of the reads
subset_se <- function(se_in = NA, 
                      multiplex_id = NA){
  se.subset <- subset(se_in, select = colData(se_in)$file %in% 
                        data.frame(colData(se_in))[data.frame(colData(se_in))$tech_replicate == multiplex_id,]$file)
  
  se.subset <- data.frame(apply(assays(se.subset)$counts, 1, sum))
  se.subset$ID <- rownames(se.subset)
  colnames(se.subset)[1] <- multiplex_id
  return(se.subset)
}
