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
#'  term, if set, must match the name of s group in "summarized$group". Default 
#'  = NULL
#' @param adjust_method Method used for multiple comparison adjustment of
#' p-values. Options are: "holm", "hochberg", "hommel", "bonferroni", "BH",
#' "BY", "fdr" or "none". See ?p.adjust.methods for a details description and
#' references. Default = "BH"
#' @param norm_method Methods for normalisation. Options are: "EDASeq" or 
#' "all_defaults". When "all_defaults" is selected, this will use all default 
#' normalisation methods for differential expression, EDASeq for QC, and edgeR 
#' "upperquantile" for determining RUV residuals (as per RUVSeq vignette). When 
#' "EDASeq" is selected, this will use EDASeq normalisation throughout. EDASeq 
#' normalisation method is selected using "EDASeq_method". Default = "EDASeq".
#' @param EDASeq_method Method for normalisation (applies to QC results using 
#' EDASeq and RUV when EDASeq is selected). Options are:"median","upper","full".
#' Default = "upper"
#' @param ruv_correct Remove Unwanted Variation (RUV)? See ?RUVr for
#' description. Currently only RUVr, which used the residuals is enabled and one
#' factor of variation is determined. If set to TRUE and a "plot_dir" is
#' provided, additional plots after RUV correction and the RUV residuals will
#' be reported. Residuals are obtained through fitting a generalised linear
#' model (GLM) using EdgeR. Residuals are then incorporated into the
#' SummarizedExperiment object and all models for DE analysis. Options = TRUE,
#' FALSE. Default = FALSE.
#' @param ensembl_annotate If the dataset has been mapped to ensembl transcript
#' identifiers, obtain additional annotation of the ensembl transcripts. A R 
#' Genome Wide Annotation object e.g. org.Mm.eg.db for mouse or org.Hs.eg.db for
#' human must be provided. Default = NULL
#' @param gtf_annotate Full path to a gtf file describing the transcripts. If 
#' provided will obtain gene symbols from gtf file. If a ensembl_annotate object
#' is also provided, this will extract annotations based on the symbols 
#' extracted from the GTF file. It is recommended to provide both a gtf file and
#' a tx_db for better annotation results. Default = NULL 
#' @param plot_dir Full path to directory for output of plots (pdf files). See
#' ?diag_plots for more details. Default = NULL
#' @param output_voom If you wish to output the results of the Voom analysis,
#' provide a full path to directory for output of files. Default = NULL
#' @param output_edger If you wish to output the results of the EdgeR analysis,
#' provide a full path to directory for output of files. Default = NULL
#' @param output_deseq If you wish to output the results of the DEseq2 analysis,
#' provide a full path to directory for output of files. Default = NULL
#' @param output_combined consensusDE will report the results of Voom, EdgeR
#' and DEseq2 as a combined report. If you wish to output the results of the
#' COMBINED analysis, provide a full path to directory for output of files. In
#' addition to the combined data, it will also output the raw count and
#' normalised data to the same directory. Default = NULL
#' @param legend Include legend in plots? Legend is based on group data in
#' summarized Options: TRUE, FALSE. Default = TRUE
#' @param label Include point labels in plots? Points are based on ID column
#' after DE analysis from merged results. Options: TRUE, FALSE. Default = TRUE

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
#' @importFrom S4Vectors DataFrame metadata metadata<-
#' @importFrom DESeq2 DESeqDataSet DESeq
#' @importFrom edgeR DGEList estimateDisp estimateGLMCommonDisp estimateGLMTagwiseDisp calcNormFactors glmFit glmLRT topTags
#' @importFrom SummarizedExperiment colData<- assays<-
#' @importFrom SummarizedExperiment assays colData
#' @importFrom EDASeq betweenLaneNormalization newSeqExpressionSet normCounts<-
#' @importFrom RUVSeq RUVr
#' @importFrom limma voom lmFit contrasts.fit eBayes topTable makeContrasts
#' @import airway
#' @import BiocGenerics
#' @importFrom data.table fread

multi_de_pairs <- function(summarized = NULL,
                           paired = "unpaired",
                           intercept = NULL,
                           adjust_method = "BH",
                           EDASeq_method = "upper",
                           norm_method = "EDASeq",
                           ruv_correct = FALSE,
                           ensembl_annotate = NULL,
                           gtf_annotate = NULL,
                           plot_dir = NULL,
                           output_voom = NULL,
                           output_edger = NULL,
                           output_deseq = NULL,
                           output_combined = NULL,
                           verbose = FALSE,
                           legend = TRUE,
                           label = TRUE){

####///---- check inputs ----\\\###
# check the p-value adjustment method is allowable
if((adjust_method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                         "fdr", "none")) == FALSE){
  stop("adjust_method must be one of the following methods: \"holm\",
       \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\" or
       \"none\"\n")
}
if((EDASeq_method %in% c("median","upper","full")) == FALSE){
  stop("EDASeq_method must be one of the following methods: \"median\",
       \"upper\", \"full\"")
}
if((norm_method %in% c("all_defaults","EDASeq")) == FALSE){
  stop("norm_method must be one of the following methods: \"all_defaults\",
       \"EDASeq\"")
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
    stop("ruv_correct can only be either \"TRUE\" or \"FALSE\". Please specify")
if(is.null(summarized)){
  if(verbose){
    message("# NO summarized experiment provided")
  }
}
# check the format of the table
sample_table <- data.frame(colData(summarized))
if("file" %in% colnames(sample_table) == FALSE){
  stop("For unpaired data, a table must be supplied with a column labelled 
       \"file\"")
}
if("group" %in% colnames(sample_table) == FALSE){
  stop("For unpaired data, a table must be supplied with a column labelled 
       \"group\"")
}
if("group" %in% colnames(sample_table) == FALSE){
  stop("For unpaired data, a table must be supplied with a column labelled 
       \"group\"")
}

if(!is.null(names(summary(sample_table$group)[summary(sample_table$group) < 2]))){
  warning("sample_table provided contained groups with less than 2 replicates!")
  warning(paste(names(summary(sample_table$group)[summary(sample_table$group) < 2]),
                collapse="\t"),
          "\tDo NOT contain biological replicates!")
}

if((paired != "paired") & (paired != "unpaired"))
  stop("paired can only be either \"paired\" or \"unpaired\". Please specify\n")
if(paired == "paired" & ("pairs" %in% colnames(sample_table) == FALSE)){
  stop("For paired data, a table must be supplied with a column labelled 
       \"pairs\"")
}

if(!is.null(gtf_annotate)){
  if(file.exists(gtf_annotate) == FALSE){
    warning("GTF file specified in gtf_annotate does not exist. Setting 
            gtf_annotate = NULL. Continuing without annotation")
    gtf_annotate <- NULL
  }
}

if(verbose){
  message("# [OPTIONS] Fitting models as: ", paired)
}

# set variables for plotting results
plot_this <- FALSE
write_this <- FALSE
# if dir is not null, see plot_this to TRUE
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
design_table <- design$table
design <- design$design

# set default
gene_coords <- NULL
# but if gene_coords exists as a meta-data
if(is.null(metadata(summarized)$gene_coords) == FALSE){
    gene_coords <- data.frame(metadata(summarized)$gene_coords)
    gene_coords$coords <- paste("chr",
                                gene_coords$seqnames,
                                ":",
                                gene_coords$start,
                                "-",
                                gene_coords$end,
                                sep="")
    gene_coords <- data.frame("ID" = gene_coords$gene_id,
                              "coords" = gene_coords$coords,
                              "strand" = gene_coords$strand,
                              "width" = gene_coords$width)
    gene_coords$coords <- gsub("chrchr", "chr", gene_coords$coords)
}

# normalise and do QC:
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
                       name = "transcript_mapped_reads",
                       legend = legend,
                       label = label))
}
# normalise for library size (i.e. remove as an effect..)
se_qc_norm <- betweenLaneNormalization(se_qc, which=EDASeq_method)
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
                       name = "NO_corrected_NO_normalised",
                       legend = legend,
                       label = label))
  invisible(diag_plots(se_in = se_qc_norm,
                       write = write_this,
                       plot_dir = plot_dir,
                       rle = plot_this,
                       pca = plot_this,
                       hclust = plot_this,
                       density = plot_this,
                       boxplot = plot_this,
                       name = "NO_corrected_YES_normalised",
                       legend = legend,
                       label = label))
}
# 3 RUV analysis
# this will basically update the design matrix and contrast matrix
# establish output of RUV corrected summarized object
summarized_ruv <- NULL

if(ruv_correct==TRUE){
  if(verbose){
    message("# [OPTIONS] RUV correction has been selected")
  }
  # determine the intercept from the base file
  ruv_se <- ruvr_correct(se = summarized,
                         plot_dir = plot_dir,
                         design = design,
                         pairing = paired,
                         intercept = intercept,
                         norm_method = norm_method,
                         EDASeq_method = EDASeq_method,
                         plot_this = plot_this,
                         write_this = write_this,
                         verbose = verbose,
                         legend = legend,
                         label = label)
  # update the summarized_ruv now
  summarized_ruv <- ruv_se$ruv_summarized
  # current only supports "W_1", i.e. k = 1 from RUVSeq
  design <- build_design(se = ruv_se$se,
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
deseq_res <- deseq_wrapper(summarized, 
                           design, 
                           contrast_matrix,
                           norm_method = norm_method,
                           EDASeq_method = EDASeq_method)
if(verbose){
  message("# Performing EdgeR analyses")
}
edger_res <- edger_wrapper(summarized, 
                           design, 
                           contrast_matrix,
                           norm_method = norm_method,
                           EDASeq_method = EDASeq_method)
if(verbose){
  message("# Performing voom analyses")
}
voom_res <- voom_wrapper(summarized, 
                         design, 
                         contrast_matrix,
                         norm_method = norm_method,
                         EDASeq_method = EDASeq_method)

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

if(!is.null(ensembl_annotate) & is.null(gtf_annotate)){
  if(verbose){
    message("# Annotating results from tx_db ONLY")
  }
  # now annotate each ensembl transcript
  merged <- lapply(seq_len(length(merged)), function(i)
                     annotate_ensembl(merged[[i]],
                                      tx_db = ensembl_annotate,
                                      coords = gene_coords
                                      ))
}
# select this over the ensembl
# need to check the format of the file
if(!is.null(gtf_annotate)){
  if(verbose){
    message("# Annotating results from gtf file")
  }
  gtf_table <- extract_gtf_attributes(gtf_annotate, verbose = verbose)
  # now annotate ensembl transcripts
  merged <- lapply(seq_len(length(merged)), function(i)
                    annotate_gtf(merged[[i]],
                                 gtf_table = gtf_table,
                                 tx_db = ensembl_annotate,
                                 coords = gene_coords))
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
                                                  vector_search = "p_intersect",
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
                       name = "MA_plots",
                       legend = legend,
                       label = label))
  if(verbose){
    message("# Producing Volcano plots. Will be located here: ", plot_dir)
  }
  invisible(diag_plots(se_in = NULL,
                       write = write_this,
                       plot_dir = plot_dir,
                       volcano = plot_this,
                       merged_in = merged_plots,
                       name = "Volcano_plots",
                       legend = legend,
                       label = label))
  if(verbose){
    message("# Producing p-value distribution plots. 
               Will be located here: ", plot_dir)
  }
  invisible(diag_plots(se_in = NULL,
                       write = write_this,
                       plot_dir = plot_dir,
                       merged_in = merged_plots,
                       p_dist = plot_this,
                       name = "p_value_distributions",
                       legend = legend,
                       label = label))
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
                      summarized_ruv = summarized_ruv,
                      merged = list("merged" = merged,
                                    "deseq" = deseq_res,
                                    "edger" = edger_res,
                                    "voom" = voom_res),
                      EDASeq_method = EDASeq_method,
                      norm_method = norm_method,
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
              "summarized" = summarized_ruv))
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
                                summarized_ruv = NULL,
                                EDASeq_method = EDASeq_method,
                                norm_method = norm_method,
                                merged = NULL,
                                voom_dir = NULL,
                                edger_dir = NULL,
                                deseq_dir = NULL,
                                merged_dir = NULL){

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
se_norm <- betweenLaneNormalization(se_new, which=EDASeq_method)
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
    merge_out$p_intersect <- apply(cbind(merge_out$PValue.x, 
                                         merge_out$PValue.y, 
                                         merge_out$PValue), 1, max)
    
    merge_out$p_union <- apply(cbind(merge_out$PValue.x, 
                                     merge_out$PValue.y, 
                                     merge_out$PValue), 1, min)
    
    merge_out$AveExpr_mean <- apply(cbind(merge_out$AveExpr.x, 
                                   merge_out$AveExpr.y, 
                                   merge_out$AveExpr), 1, mean)
    
    merge_out$LogFC_mean <- apply(cbind(merge_out$logFC.x, 
                                        merge_out$logFC.y, 
                                        merge_out$logFC), 1, mean)
    
    merge_out$LogFC_sd <- apply(cbind(merge_out$logFC.x, 
                                      merge_out$logFC.y, 
                                      merge_out$logFC), 1, sd)
    
    
    if(return=="short"){
        merge_out <- data.frame(merge_out$ID,
                                merge_out$AveExpr_mean,
                                merge_out$LogFC_mean,
                                merge_out$LogFC_sd,
                                merge_out$PValue.x,
                                merge_out$PValue.y,
                                merge_out$PValue,
                                rank(merge_out$PValue.x),
                                rank(merge_out$PValue.y),
                                rank(merge_out$PValue),
                                merge_out$sum,
                                merge_out$p_intersect,
                                merge_out$p_union)
        colnames(merge_out) <- c("ID", "AveExpr", "LogFC", "LogFC_sd",
                                  paste(name_x, "_adj_p", sep=""),
                                  paste(name_y, "_adj_p", sep=""),
                                  paste(name_z, "_adj_p", sep=""),
                                  paste(name_x, "_rank", sep=""),
                                  paste(name_y, "_rank", sep=""),
                                  paste(name_z, "_rank", sep=""),
                                  "rank_sum",
                                  "p_intersect",
                                  "p_union")

    }
    return(merge_out)
}

# determine means for a matrix and vector of common columns
table_means <- function(count_matrix, matrix_names){
  colnames(count_matrix) <- matrix_names
  # samples with less than two replicates (cannot take mean)
  # ensure will work with n=1
  l2 <- data.frame(count_matrix[,colnames(count_matrix) %in% names(summary(matrix_names)[summary(matrix_names) < 2])])
  colnames(l2) <- names(summary(matrix_names)[summary(matrix_names) < 2])
  # for samples with >= samples (take the mean)
  g2 <- count_matrix[,colnames(count_matrix) %in% names(summary(matrix_names)[summary(matrix_names) >= 2])]
  mean_counts <- data.frame("ID"=rownames(g2),
                            sapply(seq_len(length(unique(colnames(g2)))),
                                   function(i)
                                     rowMeans(g2[,colnames(g2)
                                                 %in% unique(colnames(g2))[i]])))
  mean_counts <- cbind(mean_counts, l2)
  colnames(mean_counts) <- c("ID", paste(unique(colnames(g2)),"_mean", sep=""),
                             paste(colnames(l2),"_mean", sep=""))
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
    design_table <- data.frame("file" = as.character(se$file),
                               "group" = group)
    # default is unpaired analysis
    if(pairing == "unpaired"){
        design <- stats::model.matrix(~group)
    }
    # paired, factor pairing
    if(pairing == "paired"){
        pairs <- factor(se$pairs)
        design_table <- data.frame(design_table,
                                  "pairs" = se$pairs)
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
                                     "W_1" = se$W_1)
          design <- stats::model.matrix(~W_1 + pairs + group)
      }
    }
    # format model.matrix
    rownames(design) <- colnames(se)
    colnames(design) <- gsub("group", "", colnames(design))
    colnames(design) <- gsub("pairs", "", colnames(design))
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
                             tx_db,
                             coords = NULL
                             ){
    short_ens <- gsub("\\..*", "", data_in$ID)
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

    # add coordinates IF PROVIDED
    # filter wont follow through with coords, so LH merge
    if(is.null(coords) == FALSE){
        data_in <- merge(data_in, coords, by="ID", all.x=TRUE)
    }
    # reorder the data.frame:
    data_in <- data_in[order(data_in$rank_sum, decreasing=FALSE),]
    
return(data_in)
}

# function for calling Voom analysis given se, design and contrast_matrix
voom_wrapper <- function(se, 
                         design, 
                         contrast_matrix,
                         norm_method,
                         EDASeq_method){
  
    if(norm_method == "all_defaults"){
      y <- DGEList(counts = assays(se)$counts,
                   group = se$group)
      y <- calcNormFactors(y)
    }
    if(norm_method == "EDASeq"){
      y <- newSeqExpressionSet(assays(se)$counts,
                               phenoData = data.frame(colData(se)),
                               row.names = colnames(assays(se)$counts))
      # normalise for library size
      y <- betweenLaneNormalization(y, which=EDASeq_method)
      
      y <- DGEList(counts=normCounts(y), group=y$group)
      y <- calcNormFactors(y, method="none")
    }

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
edger_wrapper <- function(se, 
                          design, 
                          contrast_matrix,
                          norm_method,
                          EDASeq_method){
  
    if(norm_method == "all_defaults"){
      y <- DGEList(counts = assays(se)$counts,
                   group = se$group)
      y <- calcNormFactors(y)
    }
    if(norm_method == "EDASeq"){
      y <- newSeqExpressionSet(assays(se)$counts,
                               phenoData = data.frame(colData(se)),
                               row.names = colnames(assays(se)$counts))
      # normalise for library size
      y <- betweenLaneNormalization(y, which=EDASeq_method)
      
      y <- DGEList(counts=normCounts(y), group=y$group)
      y <- calcNormFactors(y, method="none")
    }
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
                            data.frame("AveExpr" = all_results[[i]]$table$logCPM,
                                       "logFC" = all_results[[i]]$table$logFC,
                                       "PValue" = all_results[[i]]$table$PValue,
                                       row.names = rownames(all_results[[i]]$table)))
    names(all_results) <- names(short_results) <- colnames(contrast_matrix)
    return(list("short_results" = short_results,
                "full_results" = all_results,
                "fitted" = all_fit,
                "contrasts" = contrast_matrix))
}

# function for calling DESeq2 analysis given se, design and contrast_matrix
deseq_wrapper <- function(se, 
                          design, 
                          contrast_matrix,
                          norm_method,
                          EDASeq_method){
  
    if(norm_method == "all_defaults"){
      dds <- DESeqDataSet(se, design = design)
      # this is automatically performed here:
      # dds <- estimateSizeFactors(dds)
      dds <- DESeq(dds)
    }
    if(norm_method == "EDASeq"){
      y <- newSeqExpressionSet(assays(se)$counts,
                               phenoData = data.frame(colData(se)),
                               row.names = colnames(assays(se)$counts))
      # normalise for library size
      y <- betweenLaneNormalization(y, which=EDASeq_method)
      # put normalised counts back into se
      assays(se)$counts <- normCounts(y) 
      # put in dds
      dds <- DESeqDataSet(se, design = design)
      dds <- DESeq(dds)
      # set size factors to 1
      sizeFactors(dds) <- rep(1, length(sizeFactors(dds)))
    }
  
    # extract all if the contrasts and return all results
    all_contrasts <- lapply(seq_len(ncol(contrast_matrix)),
    function(i) DEseq_contrasts(contrast_matrix, i))
    # conduct all tests
    # full_table of results
    # here - cooksCutoff is disabled for comparability
    all_results <- lapply(seq_len(length(all_contrasts)),
    function(i) results(dds, 
                        contrast=all_contrasts[[i]]$contrast_list,
                        cooksCutoff = FALSE)
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
                         EDASeq_method = EDASeq_method,
                         norm_method = norm_method,
                         edgeR_norm = "upperquartile",
                         intercept = NULL,
                         plot_this = FALSE,
                         write_this = FALSE,
                         plot_dir = NULL,
                         legend = TRUE,
                         label = TRUE,
                         verbose = FALSE){

    ruv <- newSeqExpressionSet(assays(se)$counts,
            phenoData = data.frame(colData(se)),
            row.names = colnames(assays(se)$counts))
    if(norm_method == "EDASeq"){
      # normalise for library size
      ruv <- betweenLaneNormalization(ruv, which=EDASeq_method)
      ruv_y <- DGEList(counts=normCounts(ruv), group=ruv$group)
      ruv_y <- calcNormFactors(ruv_y, method="none")
    }
    if(norm_method == "all_defaults"){
      ruv_y <- DGEList(counts=counts(ruv), group=ruv$group)
      ruv_y <- calcNormFactors(ruv_y, method=edgeR_norm)
      
    }
    # update ruv norm counts depending on approach
    normCounts(ruv) <- ruv_y$counts
    
    # estimate dispersion parameters
    ruv_y <- estimateGLMCommonDisp(ruv_y, design)
    ruv_y <- estimateGLMTagwiseDisp(ruv_y, design)
    # fit
    fit_ruv <- glmFit(ruv_y, design)
    # obtain residuals
    res_ruv <- stats::residuals(fit_ruv, type="deviance")
    # k=1, assuming a single factor, refit with residuals
    ruv_se <- RUVr(ruv, rownames(ruv), k=1, res_ruv)
    
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
                             name="YES_corrected_YES_normalised",
                             legend = legend,
                             label = label))
        # plot model residuals
        invisible(diag_plots(se_in = ruv_se,
                             write = write_this,
                             plot_dir = plot_dir,
                             residuals = plot_this,
                             name = "RUVr_residuals",
                             legend = legend,
                             label = label))
    }

    # update colData in the SE to include the W_1
    colData(se) <- S4Vectors::DataFrame(design$table)
    colnames(se) <- colData(se)$file
    
    return_ruv_se <- ruv_se
    # put normalised counts back into se
    # assays(return_ruv_se)$normCounts <- normCounts(ruv_se) 

return(list("se" = se, 
            "design" = design,
            "ruv_summarized" = return_ruv_se))
}

extract_gtf_attributes <- function(gtf_path, verbose){
  if(verbose){
    message("# Reading gtf file in")
  }
  gtf_table <- fread(gtf_path, 
                     skip = 1, 
                     header = FALSE)
  colnames(gtf_table) <- c("chr",
                           "source",
                           "type",
                           "start",
                           "end",
                           "score",
                           "strand",
                           "phase",
                           "attributes")
  if(verbose){
    message("# Extracting gtf attributes")
  }
  # speed up
  gtf_attrib <- extract_attributes(gtf_attributes = gtf_table$attributes)
  gtf_table_out <- data.frame("ID" = as.character(gtf_attrib$ID),
                              "symbol" = as.character(gtf_attrib$geneid))
  if(verbose){
    message("# Cleaning up gtf attributes")
  }
  return(unique(gtf_table_out))
}

annotate_gtf <- function(data_in,
                         gtf_table,
                         tx_db,
                         coords){
  # merge with extracted gtf annotations and mapby symbol
  data_in <- merge(data_in, gtf_table, by="ID", all.x=TRUE)
  # then use the gtf symbols for annotation
  if(!is.null(tx_db)){
    data_in$genename <- mapIds(tx_db,
                               keys = as.character(data_in$symbol),
                               column="GENENAME",
                               keytype="SYMBOL",
                               multiVals="first")
    data_in$kegg <- mapIds(tx_db,
                           keys = as.character(data_in$symbol),
                           column="PATH",
                           keytype="SYMBOL",
                           multiVals="first")
  }
  # add coordinates IF PROVIDED
  if(!is.null(coords)){
    data_in <- merge(data_in, coords, by="ID", all.x=TRUE)
  }
  # reorder data.frame
  data_in <- data_in[order(data_in$rank_sum, decreasing=FALSE),]
  return(data_in)
}

# reformat to be vectorised
extract_attributes <- function(gtf_attributes){
  # option to extract additional gtf parameters
  gtf_id <- gsub("\"", "", 
                 gsub(";.*", "", 
                      gsub(".*gene_id ","", gtf_attributes)))
  gtf_gene <- gsub("\"", "", 
                   gsub(";.*", "", 
                        gsub(".*gene_name ","", gtf_attributes)))
  
  return(list("ID" = gtf_id,
              "geneid" = gtf_gene))
}
