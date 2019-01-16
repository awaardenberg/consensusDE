#' Batch - multiDE analysis of many comparisons
#
#' @description Given a summarized experiment generated using buildSummarized()
#' this function will automatically perform differential expression (DE)
#' analysis for all possible groups using 3 different methods 1) EdgeR, 2) Voom
#' and 3) DEseq2. It will also output 10x diagnostic plots automatically, if the
#' plotting options are selected (see ?diag_plots for more details).
#
#' @param summarized A "RangedSummarizedExperiment" object with included groups
#'  to be analysed. For format specifications see ?buildSummarized. E.g.
#'  accessible as "summarized$group". Groups are used to automate colouring of
#'  samples in unsupervised analyses. Default = NULL
#' @param paired Are the sample paired? If "paired" a paired statistical
#' analysis by including factors as pairs described in the "pairs" column of the
#' "RangedSummarizedExperiment" object in the model (accessible as
#' summarized$pairs). Options are "unpaired" or "paired". Default="unpaired"
#' @param intercept Optional ability to set the base term for fitting the model.
#'  This is not necessary as all pairs are computed automatically. The base
#'  term, if set, must match the name of group in "summarized$group". Default =
#'  NULL
#' @param adjust_method Method used for multiple comparison adjustment of
#' p-values. Options are: "holm", "hochberg", "hommel", "bonferroni", "BH",
#' "BY", "fdr" or "none". See ?p.adjust.methods for a details description and
#' references. Default = "BH"
#' @param ruv_correct Remove Unwanted Variation (RUV)? See ?RUVr for
#' description. Currently only RUVr, which used the residuals is enabled and one
#'  factor of variation is determined. If set to TRUE and a "plot_dir" is
#'  provided, additional plots after RUV correction and the RUV residuals will
#'   be reported. Residuals are obtained through fitting a generalised linear
#'   model (GLM) using EdgeR. Residuals are then incorporated into the
#'   SummarizedExperiment object and models for DE analysis. Options = TRUE,
#'   FALSE. Default = FALSE.
#' @param ensembl_annotate If the dataset has been mapped to ensembl transcript
#' identifiers, obtain additional annotation of the ensembl transcripts. A R 
#' Genome Wide Annotation object e.g. org.Mm.eg.db for mouse or org.Hs.eg.db for
#'  human must be provided. Default = NULL
#' @param plot_dir Full path to directory for output of plots (pdf files). See
#' ?diag_plots for more details. Default = NULL
#' @param output_voom If you wish to output the results of the Voom analysis,
#' provide a full path to directory for output of files. Default = NULL
#' @param output_edger If you wish to output the results of the EdgeR analysis,
#' provide a full path to directory for output of files. Default = NULL
#' @param output_deseq If you wish to output the results of the DEseq2 analysis,
#'  provide a full path to directory for output of files. Default = NULL
#' @param output_combined consensusDE will report the results of Voom, EdgeR
#' and DEseq2 as a combined report. If you wish to output the results of the
#' COMBINED analysis, provide a full path to directory for output of files. In
#' addition to the combined data, it will also output the raw count and
#' normalised data to the same directory. Default = NULL
#' @param verbose Verbosity ON/OFF. Default=FALSE
#'
#' @examples
#' ## Load the example data set and attach - see vignette for more details
#' ## The example below will perfrom DE analysis on all pairs of data
#' library(airway)
#' data(airway)
#' ## Name groups of the data.
#' colData(airway)$group <- colData(airway)$dex
#' ## Identify file locations
#' colData(airway)$file <- rownames(colData(airway))
#' #' ## Filter low count data:
#' airway_filter <- buildSummarized(summarized = airway,
#'                                  filter = TRUE)
#' ## for illustration, we only use random sample of 1000 transcripts
#' set.seed(1234)
#' airway_filter <- sample(airway_filter, 1000)
#' ## Run multi_de_pairs() with-out RUV correction
#' ## To run with RUV correction, use ruv_correct = TRUE
#' all_pairs_airway <- multi_de_pairs(summarized = airway_filter,
#'                                    ruv_correct = FALSE,
#'                                    paired = "unpaired")
#'
#' @return A list of all the comparisons conducted.
#' ## See vignette for more details.
#'
#' @export multi_de_pairs
#'
#' @importFrom AnnotationDbi mapIds keytypes columns
#' @importFrom S4Vectors DataFrame
#' @importFrom DESeq2 DESeqDataSet DESeq
#' @importFrom edgeR DGEList estimateDisp calcNormFactors glmFit glmLRT topTags
#' @importFrom SummarizedExperiment assays colData colData<-
#' @importFrom EDASeq betweenLaneNormalization newSeqExpressionSet
#' @importFrom RUVSeq RUVr
#' @importFrom limma voom lmFit contrasts.fit eBayes topTable makeContrasts
#' @import airway
#' @import BiocGenerics

multi_de_pairs <- function(summarized = NULL,
                           paired = "unpaired",
                           intercept = NULL,
                           adjust_method = "BH",
                           ruv_correct = FALSE,
                           ensembl_annotate = NULL,
                           plot_dir = NULL,
                           output_voom = NULL,
                           output_edger = NULL,
                           output_deseq = NULL,
                           output_combined = NULL,
                           verbose = FALSE){

####///---- check inputs ----\\\###
# check the p-value adjustment method is allowable
if((adjust_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                         "fdr", "none")) == FALSE){
  stop("adjust_method must be one of the following methods: \"holm\",
       \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\" or
       \"none\"\n")
}
if(is.null(plot_dir) && verbose == TRUE)
  message("No output directory provided for plots. Plots will not be generated")
# check the database
if(!is.null(ensembl_annotate)){
  if(inherits(ensembl_annotate, "OrgDb") == FALSE){
    stop("ensembl_annotate is not the correct format. ensembl_annotate needs to 
          be a Genome Wide Annotation object, e.g. org.Mm.eg.db for mouse or 
          org.Hs.eg.db for human from BioConductor")
  }
  # check that EMSEMBL is in there:
  if("ENSEMBL" %in% keytypes(ensembl_annotate) == FALSE){
    stop("ENSEMBL is not an available key in the ensembl_annotate database. 
          ensembl_annotate is not the correct format. ensembl_annotate needs to 
          be a Genome Wide Annotation object, e.g. org.Mm.eg.db for mouse or 
          org.Hs.eg.db for human from BioConductor")
  }
  # option to add other variables in future
  if(!identical(c("PATH", "SYMBOL", "GENENAME") %in% columns(ensembl_annotate), 
               c(TRUE, TRUE, TRUE))){
    stop("PATH, SYMBOL and GENENAME are not available columns in the 
          ensembl_annotate database provided. ensembl_annotate is not the 
          correct format. ensembl_annotate needs to be a Genome Wide Annotation 
          object, e.g. org.Mm.eg.db for mouse or org.Hs.eg.db for human from 
          BioConductor")
  }
}
  
if((ruv_correct != TRUE) & (ruv_correct != FALSE))
    stop("ruv_correct can only be either \"TRUE\" or \"FALSE\".
         Please specify\n")
if(is.null(summarized)){
  if(verbose){
    message("# NO summarized experiment provided")
  }
}

# check the format of the table
sample_table <- data.frame(colData(summarized))
if(ncol(sample_table)==2 &&
   (colnames(sample_table[1])!="file" ||
    colnames(sample_table[2])!="group" ||
    ncol(sample_table)<2))
  stop("For unpaired data, a table must be supplied with 2 columns.
       The first column must be labelled \"file\" and the second labelled
       \"group\"")
if(ncol(sample_table)==3 &&
   (colnames(sample_table[1])!="file" ||
    colnames(sample_table[2])!="group" ||
    colnames(sample_table[3])!="pairs" ||
    ncol(sample_table)<2))
  stop("For paired data, a table must be supplied using with 3 columns. The
        first column must be labelled \"file\", the second labelled \"group\"
        and third \"pairs\"")

if((paired != "paired") & (paired != "unpaired"))
  stop("paired can only be either \"paired\" or \"unpaired\". Please specify\n")

# set vaiables for plotting results
plot_this <- FALSE
write_this <- FALSE
# if the dir is not null, see plot.this to true
if(!is.null(plot_dir)){
  plot_this <- TRUE
  write_this <- TRUE
}

#########################################################
# Establish 1) design and 2) contrast matrix
design <- build_design(se = summarized,
                       pairing = paired,
                       intercept = intercept,
                       ruv = NULL)
contrast_matrix <- build_contrast_matrix(design$table, design$design)

# 1. build design.matrix
design_table <- data.frame("file" = as.character(summarized$file),
                           "group" = summarized$group)
design <- stats::model.matrix(~summarized$group)
#update the names of the matrix:
colnames(design) <- gsub("summarized\\$group", "", colnames(design))
colnames(design)[1] <- "Intercept"

# 2. normalise and do QC:
se_qc <- newSeqExpressionSet(assays(summarized)$counts,
                            phenoData = data.frame(colData(summarized)),
                            row.names = colnames(assays(summarized)$counts))
# plot the reads mapped/counts per sample
if(plot_this == TRUE){
  if(verbose){
    message("# Plotting Mapped Read counts. Will be located here: ", plot_dir)
  }
  invisible(diag_plots(se_in = se_qc,
                       write = write_this,
                       plot_dir = plot_dir,
                       mapped_reads = plot_this,
                       name = "transcript_mapped_reads"))
}
# normalise for library size (i.e. remove as an effect..)
se_qc_norm <- betweenLaneNormalization(se_qc, which="upper")
# output diagnostic plots:
if(plot_this == TRUE){
  if(verbose){
    message("# Plotting Dignostic plot. Will be located here: ", plot_dir)
  }
  invisible(diag_plots(se_in = se_qc,
                       write = write_this,
                       plot_dir = plot_dir,
                       rle = plot_this,
                       pca = plot_this,
                       hclust = plot_this,
                       density = plot_this,
                       boxplot = plot_this,
                       name = "NO_corrected_NO_normalised"))
  invisible(diag_plots(se_in = se_qc_norm,
                       write = write_this,
                       plot_dir = plot_dir,
                       rle = plot_this,
                       pca = plot_this,
                       hclust = plot_this,
                       density = plot_this,
                       boxplot = plot_this,
                       name = "NO_corrected_YES_normalised"))
}
# 3 RUV analysis
# this will basically update the design matrix and contrast matrix
if(ruv_correct==TRUE){
  if(verbose){
    message("# [OPTIONS] RUV batch correction has been selected")
  }

  # determine the intercept from the base file
  ruv_se <- ruvr_correct(se = summarized,
                         plot_dir = plot_dir,
                         design = design,
                         pairing = paired,
                         intercept = intercept,
                         plot_this = plot_this,
                         write_this = write_this,
                         verbose = verbose)
  # update to include RUVr weights
  design <- stats::model.matrix(~0 + W_1 + group, data = ruv_se[[2]]$table)
  colnames(design) <- gsub("group", "", colnames(design))
  
  # current only supports "W_1", i.e. k = 1 from RUVSeq
  design <- build_design(se = ruv_se[[1]],
                         pairing = paired,
                         intercept = intercept,
                         ruv = "W_1")
  contrast_matrix <- build_contrast_matrix(design$table, design$design)
  # update the colData in the SE
  # this is to include the W_1 in the summarized file for model
  colData(summarized) <- DataFrame(design$table)
  colnames(summarized) <- colData(summarized)$file
  # update design matrix
  design <- design$design
}

#########################################################
# Run DE analysis with or without RUV residuals in model
# will run against all three:
# add verbose
if(verbose){
  message("# Performing DEseq2 analyses")
}
deseq_res <- deseq_wrapper(summarized, design, contrast_matrix)
if(verbose){
  message("# Performing EdgeR analyses")
}
edger_res <- edger_wrapper(summarized, design, contrast_matrix)
if(verbose){
  message("# Performing voom analyses")
}
voom_res <- voom_wrapper(summarized, design, contrast_matrix)

#########################################################
# Merge results and report all comparisons done as well

# Test that dimensions of results from each test are the same!
# This is a unit test to check that the methods are still working
# the same.
if((length(voom_res$short_results)!=length(edger_res$short_results)) |
   (length(edger_res$short_results)!=length(deseq_res$short_results)))
  stop("Something has gone wrong! Cannot calculate non-DE set...\n
       This suggests a bug in either Voom/DeSeq/EdgeR...\n
       Contact a.waardenberg@gmail.com with error for debugging")

# test that names of contrasts are the same (i.e. outputs are consistent)
if(!identical(names(voom_res$short_results), names(edger_res$short_results)) &&
   !identical(names(edger_res$short_results), names(deseq_res$short_results)))
  stop("Something has gone wrong! Names of DE set are not the same...\n
       This suggests a bug in either Voom/DeSeq/EdgeR...\n
       Contact a.waardenberg@gmail.com with error for debugging")

# this is foolproof in case one of the methods isnt run
# this feature is for future versions (where a single method could be selected)
set_length <- max(length(voom_res$short_results),
                  length(edger_res$short_results),
                  length(deseq_res$short_results))

# set_length comes from above
# for future versions, this will determine the order
# i.e. what the output is anchored on, here it is edger.
# NB: average expr. is over all conditions.
if(verbose){
  message("# Creating merged/consensus results for edgeR, DEseq2 and Voom")
}
merged <- lapply(seq_len(set_length), function(i)
  merge_results(x = edger_res$short_results[[i]],
                y = deseq_res$short_results[[i]],
                z = voom_res$short_results[[i]],
                name_x = "edger",
                name_y = "deseq",
                name_z = "voom",
                adjust_method = adjust_method,
                return = "short"))

if(!is.null(ensembl_annotate)){
  if(verbose){
    message("# Annotating results")
  }
    merged <- lapply(seq_len(length(merged)), function(i)
                     annotate_ensembl(merged[[i]], merged[[i]]$ID,
                                      tx_db = ensembl_annotate))
}
names(merged) <- names(deseq_res$short_results)

if(plot_this == TRUE){
  if(verbose){
    message("# Plotting Dignostic plots for Merged Data. 
               Will be located here: ", plot_dir)
  }

  # select and change the column names (for compatability with diag_plots())
  merged_plots <- lapply(seq_len(length(merged)), function(i)
                             rename_colnames_list(data_in = merged[[i]],
                                                  vector_search = "edger_adj_p",
                                                  vector_replace = "Adj_PVal",
                                                  which_list = i)
  )
  names(merged_plots) <- names(merged)

  # This is where to run the additional plots for each of the comparisons
  if(verbose){
    message("# Producing MA plots. Will be located here: ", plot_dir)
  }
  invisible(diag_plots(se_in = NULL,
                       write = write_this,
                       plot_dir = plot_dir,
                       ma = plot_this,
                       merged_in = merged_plots,
                       name = "MA_plots"))
  if(verbose){
    message("# Producing Volcano plots. Will be located here: ", plot_dir)
  }
  invisible(diag_plots(se_in = NULL,
                       write = write_this,
                       plot_dir = plot_dir,
                       volcano = plot_this,
                       merged_in = merged_plots,
                       name = "Volcano_plots"))
  if(verbose){
    message("# Producing p-value distribution plots. 
               Will be located here: ", plot_dir)
  }
  invisible(diag_plots(se_in = NULL,
                       write = write_this,
                       plot_dir = plot_dir,
                       merged_in = merged_plots,
                       p_dist = plot_this,
                       name = "p_value_distributions"))
}

if(!is.null(output_voom) |
   !is.null(output_edger) |
   !is.null(output_deseq) |
   !is.null(output_combined)){

  if(verbose){
    message("# Writing all DE tables.")
    if(!is.null(output_voom)){
      message("# Writing Voom DE tables to: ", output_voom)
    }
    if(!is.null(output_edger)){
      message("# Writing EdgeR DE tables to: ", output_edger)
    }
    if(!is.null(output_deseq)){
      message("# Writing DEseq2 DE tables to: ", output_deseq)
    }
    if(!is.null(output_combined)){
      message("# Writing merged DE tables to: ", output_combined)
    }
  }
  # invoke write_table_wrapper
  write_table_wrapper(summarized = summarized,
                      merged = list("merged" = merged,
                                    "deseq" = deseq_res,
                                    "edger" = edger_res,
                                    "voom" = voom_res),
                      voom_dir = output_voom,
                      edger_dir = output_edger,
                      deseq_dir = output_deseq,
                      merged_dir = output_combined)
}

if(verbose){
  message("# Finished DE analysis")
}

# returns the merged results of the 3 comparisons.

if(ruv_correct == TRUE){
  return(list("merged" = merged,
              "deseq" = deseq_res,
              "edger" = edger_res,
              "voom" = voom_res,
              "summarized" = ruv_se[[1]]))
}

if(ruv_correct == FALSE){
  return(list("merged" = merged,
              "deseq" = deseq_res,
              "edger" = edger_res,
              "voom" = voom_res,
              "summarized" = summarized))
}

}

# helper function for annotating and writing tables:
# these are seperate commands for writing tables:
write_table_wrapper <- function(summarized = NULL,
                                merged = NULL,
                                voom_dir = NULL,
                                edger_dir = NULL,
                                deseq_dir = NULL,
                                merged_dir = NULL){

# format checks
  # check the merged is a list format and check the names! etc.
#########################################################
# obtain raw and normalised counts
raw <- assays(summarized)$counts
cpm_raw <- cpm(assays(summarized)$counts)
cpm_raw_log <- cpm(assays(summarized)$counts, log=TRUE)
cpm_mean <- table_means(cpm_raw, colData(summarized)$group)
cpm_mean_log <- table_means(cpm_raw_log, colData(summarized)$group)

# to normalise, push into seqexpressionset
se_new <- newSeqExpressionSet(assays(summarized)$counts,
                              phenoData = data.frame(colData(summarized)),
                              row_names=colnames(assays(summarized)$counts))
# normalise for library size
se_norm <- betweenLaneNormalization(se_new, which="upper")
norm_count <- normCounts(se_norm)
cpm_norm <- cpm(normCounts(se_norm))
cpm_norm_log <- cpm(normCounts(se_norm), log=TRUE)
cpm_norm_mean <- table_means(cpm_norm, colData(summarized)$group)
cpm_norm_log_mean <- table_means(cpm_norm_log, colData(summarized)$group)

# directories to print to
# check exists first!
if(!is.null(merged_dir)){
  utils::write.table(raw, file=paste(merged_dir, "raw_count.tsv", sep=""),
                     sep="\t")
  utils::write.table(cpm_raw, file=paste(merged_dir, "raw_cpm.tsv", sep=""),
                     sep="\t")
  utils::write.table(cpm_raw_log, file=paste(merged_dir, "raw_log_cpm.tsv",
                                             sep=""), sep="\t")
  utils::write.table(cpm_mean, file=paste(merged_dir, "raw_cpm_mean.tsv",
                                          sep=""), sep="\t", row.names=FALSE)
  utils::write.table(cpm_mean_log, file=paste(merged_dir,
                                              "raw_log_cpm_mean.tsv", sep=""),
                     sep="\t", row.names=FALSE)

  utils::write.table(norm_count, file=paste(merged_dir, "normalised_count.tsv",
                                            sep=""), sep="\t")
  utils::write.table(cpm_norm, file=paste(merged_dir, "normalised_cpm.tsv",
                                          sep=""), sep="\t")
  utils::write.table(cpm_norm_log, file=paste(merged_dir,
                                              "normalised_log_cpm.tsv", sep=""),
                     sep="\t")
  utils::write.table(cpm_norm_mean, file=paste(merged_dir,
                                               "normalised_cpm_mean.tsv",
                                               sep=""), sep="\t")
  utils::write.table(cpm_norm_log_mean, file=paste(merged_dir,
                                                   "raw_log_cpm_mean.tsv",
                                                   sep=""), sep="\t",
                     row.names=FALSE)

  # COMBINE WITH cpm_norm_mean [ normalised counts per million ] [ NOT LOGGED ]
  invisible(sapply(seq_len(length(merged$merged)), function(i)
    utils::write.table(data.frame(merge(cpm_norm_mean,
                                 merged$merged[[i]],
                                 by="ID"))[order(data.frame(merge(cpm_norm_mean,
                                                            merged$merged[[i]],
                                                            by="ID"))$rank_sum,
                                                 decreasing=FALSE), ],
                file=paste(merged_dir, names(merged$merged)[i],
                           "_combined_results.tsv", sep=""),
                sep="\t",
                row.names=FALSE)))
}

if(!is.null(deseq_dir)){
  # DESEQ:
  invisible(sapply(seq_len(length(merged$deseq$full_results)), function(i)
    utils::write.table(data.frame(merged$deseq$full_results[[i]]),
                file=paste(deseq_dir, names(merged$deseq$full_results[i]),
                           "_deseq_results.tsv", sep=""),
                sep = "\t",
                row.names = TRUE,
                col.names = NA)))
}
if(!is.null(voom_dir)){
  # VOOM
  invisible(sapply(seq_len(length(merged$voom$full_results)), function(i)
    utils::write.table(data.frame(merged$voom$full_results[[i]]),
                file=paste(voom_dir, names(merged$voom$full_results[i]),
                           "_voom_results.tsv", sep=""),
                sep = "\t",
                row.names = TRUE,
                col.names = NA)))
}

if(!is.null(edger_dir)){
  # EDGER
  invisible(sapply(seq_len(length(merged$edger$full_results)), function(i)
    utils::write.table(data.frame(merged$edger$full_results[[i]]),
                file=paste(edger_dir, names(merged$edger$full_results[i]),
                           "_edger_results.tsv", sep=""),
                sep = "\t",
                row.names = TRUE,
                col.names = NA)))
}
}
#########################################################

# helper function for renaming column names in a list
rename_colnames_list <- function(data_in = NULL,
                                 vector_search = NULL,
                                 vector_replace = NULL,
                                 which_list = NULL
                                 ){
  data_out <- data_in
  colnames(data_out) <- gsub(vector_search, vector_replace, colnames(data_out))
return(data_out)
}

# intersect function (for intersecting up to 3 vectors)
# reports intersect size and common set
# not used from version 0.05
intersect_test <- function(x=NULL, y=NULL, z=NULL){
    if(!is.null(x) & !is.null(y)){
        message("# set1 = ", length(x))
        message("# set2 = ", length(y))
        common <- intersect(x, y)
    }
    if(!is.null(x) & !is.null(y) & !is.null(z)){
        message("# set3 = ", length(z))
        common <- intersect(common, z)
    }
    message("# intersect = ", length(common))
return(common)
}

# function for obtaining non-de rows from the outputs of the XX$short_results
# from the DE analysis:
non_de_rows <- function(input_data, p_cut=0.5){
    data_out <- rownames(input_data[input_data$PValue > p_cut,])
    return(data_out)
}

# function for merging 3 DE tables in "short" format
# returns merged table with ranks for each
merge_results <- function(x, y, z,
name_x, name_y, name_z,
adjust_method="BH",
return="short"){
    x <- data.frame("ID"=rownames(x), x)
    y <- data.frame("ID"=rownames(y), y)
    z <- data.frame("ID"=rownames(z), z)

    merge_out <- merge(x, y, by="ID")
    merge_out <- merge(merge_out, z, by="ID")
    # adjust p-values
    merge_out$PValue.x <- stats::p.adjust(merge_out$PValue.x, adjust_method)
    merge_out$PValue.y <- stats::p.adjust(merge_out$PValue.y, adjust_method)
    merge_out$PValue <- stats::p.adjust(merge_out$PValue, adjust_method)
    # extract PValue columns for ranking.
    rank_this <- cbind(merge_out$PValue.x,
                       merge_out$PValue.y,
                       merge_out$PValue)
    # sum of ranks
    merge_out$sum <- rowSums(sapply(seq_len(ncol(rank_this)),
                                    function(i) rank(rank_this[,i])))
    
    
    # order results by sum of ranks
    merge_out <- merge_out[order(merge_out$sum, decreasing=FALSE),]
    
    # determine which is the max p value (or contains least probability)
    merge_out$p_max <- sapply(seq_len(nrow(merge_out)), function(i)
                        max(c(merge_out$PValue.x[i], 
                              merge_out$PValue.y[i], 
                              merge_out$PValue[i])))


    if(return=="short"){
        merge_out <- data.frame(merge_out$ID,
        merge_out$AveExpr,
        merge_out$logFC.x,
        merge_out$PValue.x,
        merge_out$PValue.y,
        merge_out$PValue,
        rank(merge_out$PValue.x),
        rank(merge_out$PValue.y),
        rank(merge_out$PValue),
        merge_out$sum,
        merge_out$p_max)
        colnames(merge_out) <- c("ID", "AveExpr", "LogFC",
        paste(name_x, "_adj_p", sep=""),
        paste(name_y, "_adj_p", sep=""),
        paste(name_z, "_adj_p", sep=""),
        paste(name_x, "_rank", sep=""),
        paste(name_y, "_rank", sep=""),
        paste(name_z, "_rank", sep=""),
        "rank_sum",
        "p_max")

    }
    return(merge_out)
}

# determine means for a matrix and vector of common columns
table_means <- function(count_matrix, matrix_names){
    colnames(count_matrix) <- matrix_names
    mean_counts <- data.frame("ID"=rownames(count_matrix),
                        sapply(seq_len(length(unique(colnames(count_matrix)))),
                               function(i)
                                 rowMeans(count_matrix[,colnames(count_matrix)
                                    %in% unique(colnames(count_matrix))[i]])))
    colnames(mean_counts) <- c("ID", paste(unique(colnames(count_matrix)),
                                           "_mean", sep=""))
    return(mean_counts)
}

# function to automate contrast matrix creation of all pairs
build_contrast_matrix <- function(table_design, design){
    # base will already be levelled - so take the first level
    base_intercept <- levels(table_design$group)[1]
    # what is remaining in the groups (other pairs)
    remaining <- setdiff(levels(table_design$group), base_intercept)
    # if more than two comparisons
    if(length(levels(table_design$group)) > 2){
        # all possible combinations (of remaining pairs):
        all_com <- utils::combn(remaining, 2)
        remaining_com <- c(sapply(seq_len(ncol(all_com)), function(i)
        paste(all_com[1,i], all_com[2,i], sep="-")))
        # names for the comparisons
        real_com <- c(sapply(seq_len(length(remaining)), function(i)
                             paste(remaining[i], base_intercept, sep="-")))
        all_c <- paste(c(remaining, remaining_com), collapse=",")
        contrast_cmd <- paste("makeContrasts(", all_c, ",levels=design)",
                              sep="")
        # evaluate command string
        contrast_matrix <- eval(parse(text=contrast_cmd))
        colnames(contrast_matrix) <- c(real_com, remaining_com)
    }
    # if exactly two comparisons
    if(length(levels(table_design$group)) == 2){
        all_c <- levels(table_design$group)[2]
        contrast_cmd <- paste("makeContrasts(", all_c, ",levels=design)",
                              sep="")
        contrast_matrix <- eval(parse(text=contrast_cmd))
        colnames(contrast_matrix) <- paste(levels(table_design$group)[2],
        levels(table_design$group)[1], sep="-")
    }
    return(contrast_matrix)
}



# function to create 1) design matrix and 2) design table
# works with W_1 from RUV only
build_design <- function(se = NULL,
                         pairing = NULL,
                         intercept = NULL,
                         ruv = NULL){
    group <- factor(se$group)
    # set the intercept term/base coeffecient for comparisons.
    if(!is.null(intercept)){
        group <- stats::relevel(group, intercept)
    }
    design_table <- data.frame("file"=as.character(se$file),
                               "group"=group)
    # default is unpaired analysis
    if(pairing=="unpaired"){
        design <- stats::model.matrix(~group)
    }
    # paired, factor pairing
    if(pairing=="paired"){
        pairs <- factor(se$pairs)
        design_table <- data.frame(design_table,
        "pairs"=se$pairs)
        design <- stats::model.matrix(~pairs + group)
    }
    # for RUV n=1
    if(!is.null(ruv)){
      if(pairing=="unpaired" && ruv=="W_1"){
          W_1 <- as.numeric(se$W_1)
          design_table <- data.frame(design_table,
          "W_1"=se$W_1)
          design <- stats::model.matrix(~W_1 + group)
      }
      if(pairing=="paired" && ruv=="W_1"){
          W_1 <- as.numeric(se$W_1)
          design_table <- data.frame(design_table,
          "W_1"=se$W_1)
          design <- stats::model.matrix(~W_1 + pairs + group)
      }
    }
    # format model.matrix
    rownames(design) <- colnames(se)
    colnames(design) <- gsub("group", "", colnames(design))
    colnames(design)[1] <- "Intercept"
return(list("design" = design, "table" = design_table))
}


# function for building DEseq compatabible contrast list from
# a contrast_matrix built using model_matrix()
DEseq_contrasts <- function(contrast_matrix, n){
    # n is the contrast of interest
    contrast_i <- contrast_matrix[,n]
    # check for positive and negative comparisons for each in the match:
    contrast_i_pos <- names(contrast_i[contrast_i==1])
    contrast_i_neg <- names(contrast_i[contrast_i==-1])
    contrast_i_name <- colnames(contrast_matrix)[n]
    # DEseq contrast lists are: nominator/denominator
    contrast_i_list <- list(c(contrast_i_pos), c(contrast_i_neg))
    return(list("contrast_list"=contrast_i_list, "contrast_name"=
                  contrast_i_name))
}


# annotation wrapper - works with "OrgDb" object from BioConductor
annotate_ensembl <- function(data_in,
                             ids,
                             tx_db){

    short_ens <- gsub("\\..*", "", ids)
    data_in$genename <- mapIds(tx_db,
                              keys = short_ens,
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first")

    data_in$symbol <- mapIds(tx_db,
                            keys = short_ens,
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
    data_in$kegg <- mapIds(tx_db,
                          keys = short_ens,
                          column="PATH",
                          keytype="ENSEMBL",
                          multiVals="first")
return(data_in)
}

# function for calling Voom analysis given se, design and contrast_matrix
voom_wrapper <- function(se, design, contrast_matrix){
    y <- DGEList(counts = assays(se)$counts,
                 group = se$group)
    y <- calcNormFactors(y)
    v <- voom(y, design = design)
    fit_voom <- lmFit(v, design = design)
    # fitting with user supplied contrast_matrix for each comparison
    all_fit <- lapply(seq_len(ncol(contrast_matrix)), function(i)
                      contrasts.fit(fit = fit_voom,
                              contrasts = t(data.frame(t(contrast_matrix[,i]),
                                  row.names = colnames(contrast_matrix)[i]))))

    all_fit <- lapply(seq_len(length(all_fit)), function(i)
    eBayes(all_fit[[i]]))

    all_results <- lapply(seq_len(length(all_fit)), function(i)
    topTable(all_fit[[i]], number=nrow(all_fit[[i]]$coefficients)))
    # short table for merging with other methods
    short_results <- lapply(seq_len(length(all_results)), function(i)
    data.frame("AveExpr"=all_results[[i]]$AveExpr,
    "logFC"=all_results[[i]]$logFC,
    "PValue"=all_results[[i]]$P.Value,
    row.names=rownames(all_results[[i]])))
    names(all_results) <- names(short_results) <- colnames(contrast_matrix)
    return(list("short_results"=short_results,
    "full_results"=all_results,
    "fitted"=all_fit,
    "contrasts"=contrast_matrix))
}

# function for calling DESeq2 analysis given se, design and contrast_matrix
edger_wrapper <- function(se, design, contrast_matrix){

    y <- DGEList(counts=assays(se)$counts,
    group=se$group)
    # normalise accross groups
    y <- calcNormFactors(y)
    # NB. there are different ways of estimating dispersion
    y_disp <- estimateDisp(y, design)
    fit_edger <- glmFit(y_disp, design)
    # fit each comparison from contrast_matrix of interest
    all_fit <- lapply(seq_len(ncol(contrast_matrix)), function(i)
    glmLRT(fit_edger, contrast=as.numeric(contrast_matrix[,i])))
    all_results <- lapply(seq_len(length(all_fit)), function(i)
    topTags(all_fit[[i]], n=nrow(all_fit[[i]]$coefficients)))
    # short table for merging with other methods
    short_results <- lapply(seq_len(length(all_results)), function(i)
    data.frame("AveExpr"=all_results[[i]]$table$logCPM,
    "logFC"=all_results[[i]]$table$logFC,
    "PValue"=all_results[[i]]$table$PValue,
    row.names=rownames(all_results[[i]]$table)))
    names(all_results) <- names(short_results) <- colnames(contrast_matrix)
    return(list("short_results"=short_results,
    "full_results"=all_results,
    "fitted"=all_fit,
    "contrasts"=contrast_matrix))
}

# function for calling DESeq2 analysis given se, design and contrast_matrix
deseq_wrapper <- function(se, design, contrast_matrix){

    dds <- DESeqDataSet(se, design = design)
    dds <- DESeq(dds)
    # extract all if the contrasts and return all results
    all_contrasts <- lapply(seq_len(ncol(contrast_matrix)),
    function(i) DEseq_contrasts(contrast_matrix, i))
    # conduct all tests
    # full_table of results
    all_results <- lapply(seq_len(length(all_contrasts)),
    function(i) results(dds, contrast=all_contrasts[[i]]$contrast_list)
    )
    names(all_results) <- sapply(seq_len(length(all_contrasts)), function(i)
    all_contrasts[[i]]$contrast_name)
    # short table for merging with other methods
    short_results <- lapply(seq_len(length(all_results)), function(i)
    data.frame("AveExpr"=log2(all_results[[i]]$baseMean),
    "logFC"=all_results[[i]]$log2FoldChange,
    "PValue"=all_results[[i]]$pvalue,
    row.names=rownames(all_results[[i]])))
    names(short_results) <- names(all_results)
    return(list("short_results"=short_results,
    "full_results"=all_results,
    "fitted"=dds,
    "contrasts"=contrast_matrix))
}

#########################################################
# Perform RUVr correction (OPTIONAL)
# - if correction with RUVr is used, update model

ruvr_correct <- function(se = NULL,
                         design = NULL,
                         pairing = "unpaired",
                         intercept = NULL,
                         plot_this = FALSE,
                         write_this = FALSE,
                         plot_dir = NULL,
                         verbose = FALSE){

    ruv <- newSeqExpressionSet(assays(se)$counts,
            phenoData = data.frame(colData(se)),
            row.names = colnames(assays(se)$counts))
    # normalise for library size (i.e. remove as an effect..)
    ruv_norm <- betweenLaneNormalization(ruv, which="upper")

    # fit GLM for estimation of W_1:
    ruv_y <- DGEList(counts=counts(ruv_norm), group=ruv_norm$group)
    ruv_y <- calcNormFactors(ruv_y, method="upperquartile")
    ruv_y <- estimateDisp(ruv_y, design)
    fit_ruv <- glmFit(ruv_y, design)
    res_ruv <- stats::residuals(fit_ruv, type="deviance")
    # >>>>this is a questionable parameter, using k=1
    # Assuming only a batch factor to remove
    ruv_se <- RUVr(ruv_norm, rownames(ruv_norm), k=1, res_ruv)
    # build design and contrast matrix for refitting model
    design <- build_design(se=ruv_se,
                           pairing=pairing,
                           intercept=intercept,
                           ruv="W_1")
    # output RLE and PCA plots
    if(plot_this == TRUE){
        if(verbose){
          message("# Plotting RUV Dignostic plots. 
                     Will be located here:", plot_dir)
        }

        invisible(diag_plots(se_in = ruv_se,
                             write = write_this,
                             plot_dir = plot_dir,
                             rle=plot_this,
                             pca=plot_this,
                             hclust=plot_this,
                             density=plot_this,
                             boxplot=plot_this,
                             name="YES_corrected_YES_normalised"))
        # plot model residuals
        invisible(diag_plots(se_in = ruv_se,
                             write = write_this,
                             plot_dir = plot_dir,
                             residuals = plot_this,
                             name = "RUVr_residuals"))
    }

    # update colData in the SE to include the W_1
    colData(se) <- S4Vectors::DataFrame(design$table)
    colnames(se) <- colData(se)$file

    return(list(se, design))
}

