#' Batch - multiDE analysis of many comparisons
#
#' @description Given a summarised experiment generated using read.summarised()
#' this function will automatically perform differential expression (DE)
#' analysis for all possible groups using 3 different methods 1) EdgeR, 2) Voom
#' and 3) DEseq2. It will also output 10x diagnostic plots automatically, if the
#'  plotting options are selected (see ?diag.plots for more details).
#
#' @param summarised A "RangedSummarizedExperiment" object with included groups
#'  to be analysed. For format specifications see ?read.summarised. E.g.
#'  accessible as "summarised$group". Groups are used to automate colouring of
#'  samples in unsupervised analyses. Default = NULL
#' @param paired Are the sample paired? If "paired" a paired statistical
#' analysis by including factors as pairs described in the "pairs" column of the
#' "RangedSummarizedExperiment" object in the model (accessible as
#' summarised$pairs). Options are "unpaired" or "paired". Default="unpaired"
#' @param intercept Optional ability to set the base term for fitting the model.
#'  This is not necessary as all pairs are computed automatically. The base
#'  term, if set, must match the name of group in "summarised$group". Default =
#'  NULL
#' @param adjust.method Method used for multiple comparison adjustment of
#' p-values. Options are: "holm", "hochberg", "hommel", "bonferroni", "BH",
#' "BY", "fdr" or "none". See ?p.adjust.methods for a details description and
#' references. Default = "BH"
#' @param ruv.correct Remove Unwanted Variation (RUV)? See ?RUVr for
#' description. Currently only RUVr, which used the residuals is enabled and one
#'  factor of variation is determined. If set to TRUE and a "plot.dir" is
#'  provided, additional plots after RUV correction and the RUV residuals will
#'   be reported. Residuals are obtained through fitting a generalised linear
#'   model (GLM) using EdgeR. Residuals are then incorporated into the
#'   SummarizedExperiment object and models for DE analysis. Options = TRUE,
#'   FALSE. Default = FALSE.
#' @param ensembl.annotate If the dataset has been mapped to ensembl transcript
#' identifiers, obtain additional annotation of the ensembl transcripts. If set
#' to TRUE, a tx.db must also be provided (see below). Default = FALSE
#' @param tx.db An R txdb object. E.g. TxDb.Dmelanogaster.UCSC.dm3.ensGene. Must
#'  be provided if ensembl.annotate is set to TRUE. Default = NULL
#' @param plot.dir Full path to directory for output of plots (pdf files). See
#' ?diag.plots for more details. Default = NULL
#' @param output.voom If you wish to output the results of the Voom analysis,
#' provide a full path to directory for output of files. Default = NULL
#' @param output.edger If you wish to output the results of the EdgeR analysis,
#' provide a full path to directory for output of files. Default = NULL
#' @param output.deseq If you wish to output the results of the DEseq2 analysis,
#'  provide a full path to directory for output of files. Default = NULL
#' @param output.combined consensusDE will report the results of Voom, EdgeR
#' and DEseq2 as a combined report. If you wish to output the results of the
#' COMBINED analysis, provide a full path to directory for output of files. In
#' addition to the combined data, it will also output the raw count and
#' normalised data to the same directory. Default = NULL
#' @param verbose Verbosity ON/OFF. Default=FALSE
#'
#' @examples
#' ## Load the example data set and attach - see vignette for more details
#' cat("The example below will perfrom DE analysis on all pairs of data")
#' \dontrun{
#' library(airway)
#' data("airway")
#' ## Name the groups of the data.
#' colData(airway)$group <- colData(airway)$dex
#' ## Identify the file locations
#' colData(airway)$file <- rownames(colData(airway))
#' #' ## Filter low count data:
#' airway.filter <- read.summarised(summarised = airway,
#'                                  filter = TRUE)
#' ## Run multi.de.pairs() with-out RUV correction
#' ## To run with RUV correction, use ruv.correct = TRUE
#' all.pairs.airway <- multi.de.pairs(summarised = airway.filter,
#'                                    ruv.correct = FALSE,
#'                                    paired = "unpaired")
#' }
#' @return A list of all the comparisons conducted.
#' ## See vignette for more details.
#'
#' @export multi.de.pairs
#'
#' @importFrom AnnotationDbi mapIds
#' @importFrom DESeq2 DESeqDataSet DESeq
#' @importFrom edgeR DGEList estimateDisp calcNormFactors glmFit glmLRT topTags
#' @importFrom SummarizedExperiment assays colData colData<-
#' @importFrom EDASeq betweenLaneNormalization
#' @importFrom RUVSeq RUVr
#' @importFrom limma voom lmFit contrasts.fit eBayes topTable makeContrasts

multi.de.pairs <- function(summarised = NULL,
                           paired = "unpaired",
                           intercept = NULL,
                           adjust.method = "BH",
                           ruv.correct = FALSE,
                           ensembl.annotate = FALSE,
                           tx.db = NULL,
                           plot.dir = NULL,
                           output.voom = NULL,
                           output.edger = NULL,
                           output.deseq = NULL,
                           output.combined = NULL,
                           verbose = FALSE){

####///---- check inputs ----\\\###
# check the p-value adjustment method is allowable
if((adjust.method %in% c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
                         "fdr", "none")) == FALSE){
  stop("adjust.method must be one of the following methods: \"holm\",
       \"hochberg\", \"hommel\", \"bonferroni\", \"BH\", \"BY\", \"fdr\" or
       \"none\"\n")
}

if(is.null(plot.dir) && verbose == TRUE)
  cat("No output directory provided for plots. Plots will not be generated\n")

if(ensembl.annotate == TRUE & is.null(tx.db))
  warning("ensembl.annotate was set to TRUE, but no tx.db detected. Continuing
           without annotation of data\n")
if((ruv.correct != TRUE) & (ruv.correct != FALSE))
    stop("ruv.correct can only be either \"TRUE\" or \"FALSE\".
         Please specify\n")
if(is.null(summarised)){
  if(verbose){
    cat(paste("# NO summarised experiment provided", "\n", sep=""))
  }
}

# check the format of the table
sample.table <- data.frame(colData(summarised))
if(ncol(sample.table)==2 &&
   (colnames(sample.table[1])!="file" ||
    colnames(sample.table[2])!="group" ||
    ncol(sample.table)<2))
  stop("For unpaired data, a table must be supplied with 2 columns.
       The first column must be labelled \"file\" and the second labelled
       \"group\"")
if(ncol(sample.table)==3 &&
   (colnames(sample.table[1])!="file" ||
    colnames(sample.table[2])!="group" ||
    colnames(sample.table[3])!="pairs" ||
    ncol(sample.table)<2))
  stop("For paired data, a table must be supplied using with 3 columns. The
        first column must be labelled \"file\", the second labelled \"group\"
        and third \"pairs\"")

if((paired != "paired") & (paired != "unpaired"))
  stop("paired can only be either \"paired\" or \"unpaired\". Please specify\n")

# set vaiables for plotting results
plot.this <- FALSE
write.this <- FALSE
# if the dir is not null, see plot.this to true
if(!is.null(plot.dir)){
  plot.this <- TRUE
  write.this <- TRUE
}

#########################################################
# Establish 1) design and 2) contrast matrix
design <- build.design(se = summarised,
                       pairing = paired,
                       intercept = intercept,
                       ruv = FALSE)
contrast.matrix <- build.contrast.matrix(design$table, design$design)

# 1. build design.matrix
design.table <- data.frame("file" = as.character(summarised$file),
                           "group" = summarised$group)
design <- stats::model.matrix(~summarised$group)
#update the names of the matrix:
colnames(design) <- gsub("summarised\\$group", "", colnames(design))
colnames(design)[1] <- "Intercept"

# 2. normalise and do QC:
se.qc <- newSeqExpressionSet(assays(summarised)$counts,
                            phenoData = data.frame(colData(summarised)),
                            row.names = colnames(assays(summarised)$counts))
# plot the reads mapped/counts per sample
if(plot.this == TRUE){
  if(verbose){
    cat(paste("# Plotting Mapped Read counts. Will be located here: ", plot.dir,
              "\n", sep=""))
  }
  invisible(diag.plots(se.in = se.qc,
                       write = write.this,
                       plot.dir = plot.dir,
                       mapped.reads = plot.this,
                       name = "transcript_mapped_reads"))
}
# normalise for library size (i.e. remove as an effect..)
se.qc.norm <- betweenLaneNormalization(se.qc, which="upper")
# output diagnostic plots:
if(plot.this == TRUE){
  if(verbose){
    cat(paste("# Plotting Dignostic plot. Will be located here: ", plot.dir,
              "\n", sep=""))
  }
  invisible(diag.plots(se.in = se.qc,
                       write = write.this,
                       plot.dir = plot.dir,
                       rle = plot.this,
                       pca = plot.this,
                       hclust = plot.this,
                       density = plot.this,
                       boxplot = plot.this,
                       name = "NO_corrected_NO_normalised"))
  invisible(diag.plots(se.in = se.qc.norm,
                       write = write.this,
                       plot.dir = plot.dir,
                       rle = plot.this,
                       pca = plot.this,
                       hclust = plot.this,
                       density = plot.this,
                       boxplot = plot.this,
                       name = "NO_corrected_YES_normalised"))
}
# 3 RUV analysis
# this will basically update the design matrix and contrast matrix
if(ruv.correct==TRUE){
  if(verbose){
    cat("# [OPTIONS] RUV batch correction has been selected")
  }

  # determine the intercept from the base file (this is done in design.matrix.R)
  ruv.se <- ruvr.correct(se = summarised,
                         plot.dir = plot.dir,
                         design = design,
                         pairing = paired,
                         intercept = intercept,
                         plot.this = plot.this,
                         write.this = write.this,
                         verbose = verbose)
  # update to include RUVr weights
  design <- stats::model.matrix(~0 + W_1 + group, data = ruv.se[[2]]$table)
  colnames(design) <- gsub("group", "", colnames(design))

  design <- build.design(se = ruv.se[[1]],
                         pairing = paired,
                         intercept = intercept,
                         ruv = TRUE)
  contrast.matrix <- build.contrast.matrix(design$table, design$design)
  design <- design$design
}

#########################################################
# Run DE analysis with or without RUV residuals in model
# will run against all three:
# add verbose
if(verbose){
  cat("# Performing DEseq2 analyses\n")
}
deseq.res <- deseq.wrapper(summarised, design, contrast.matrix)
if(verbose){
  cat("# Performing EdgeR analyses\n")
}
edger.res <- edger.wrapper(summarised, design, contrast.matrix)
if(verbose){
  cat("# Performing voom analyses\n")
}
voom.res <- voom.wrapper(summarised, design, contrast.matrix)

#########################################################
# Merge results and report all comparisons done as well

# Test that dimensions of results from each test are the same!
# This is a unit test to check that the methods are still working
# the same.
if((length(voom.res$short.results)!=length(edger.res$short.results)) |
   (length(edger.res$short.results)!=length(deseq.res$short.results)))
  stop("Something has gone wrong! Cannot calculate non-DE set...\n
       This suggests a bug in either Voom/DeSeq/EdgeR...\n
       Contact ashley.waardenberg@jcu.edu.au with error for debugging")

# test that names of contrasts are the same (i.e. outputs are consistent)
if(!identical(names(voom.res$short.results), names(edger.res$short.results)) &&
   !identical(names(edger.res$short.results), names(deseq.res$short.results)))
  stop("Something has gone wrong! Names of DE set are not the same...\n
       This suggests a bug in either Voom/DeSeq/EdgeR...\n
       Contact ashley.waardenberg@jcu.edu.au with error for debugging")

# this is foolproof in case one of the methods isnt run
# this feature is for future versions (where a single method could be selected)
set.length <- max(length(voom.res$short.results),
                  length(edger.res$short.results),
                  length(deseq.res$short.results))

# set.length comes from above
# for future versions, this will determine the order
# i.e. what the output is anchored on, here it is edger.
# NB: average expr. is over all conditions.
if(verbose){
  cat("# Creating merged/consensus results for edgeR, DEseq2 and Voom\n")
}
merged <- lapply(seq_len(set.length), function(i)
  merge.results(x = edger.res$short.results[[i]],
                y = deseq.res$short.results[[i]],
                z = voom.res$short.results[[i]],
                name.x = "edger",
                name.y = "deseq",
                name.z = "voom",
                adjust.method = adjust.method,
                return = "short"))

if(ensembl.annotate==TRUE & !is.null(tx.db)){
  if(verbose){
    cat("# Annotating results\n")
  }
    merged <- lapply(seq_len(length(merged)), function(i)
                     annotate.ensembl(merged[[i]], merged[[i]]$ID,
                                      tx.db = tx.db))
}
names(merged) <- names(deseq.res$short.results)

if(plot.this == TRUE){
  if(verbose){
    cat(paste("# Plotting Dignostic plots for Merged Data.
              Will be located here: ", plot.dir, "\n", sep=""))
  }

  # select and change the column names (for compatability with diag.plots())
  merged.plots <- lapply(seq_len(length(merged)), function(i)
                             rename.colnames.list(data.in = merged[[i]],
                                                  vector.search = "edger.adj.p",
                                                  vector.replace = "Adj.PVal",
                                                  which.list = i)
  )
  names(merged.plots) <- names(merged)

  # This is where to run the additional plots for each of the comparisons
  if(verbose){
    cat(paste("# Producing MA plots. Will be located here: ", plot.dir, "\n",
              sep=""))
  }
  invisible(diag.plots(se.in=NULL,
                       write=write.this,
                       plot.dir=plot.dir,
                       ma=plot.this,
                       merged.in=merged.plots,
                       name="MA_plots"))
  if(verbose){
    cat(paste("# Producing Volcano plots. Will be located here: ", plot.dir,
              "\n", sep=""))
  }
  invisible(diag.plots(se.in = NULL,
                       write = write.this,
                       plot.dir = plot.dir,
                       volcano = plot.this,
                       merged.in = merged.plots,
                       name = "Volcano_plots"))
  if(verbose){
    cat(paste("# Producing p-value distribution plots. Will be located here: ",
              plot.dir, "\n", sep=""))
  }
  invisible(diag.plots(se.in = NULL,
                       write = write.this,
                       plot.dir = plot.dir,
                       merged.in = merged.plots,
                       p.dist = plot.this,
                       name = "p_value_distributions"))
}

if(!is.null(output.voom) |
   !is.null(output.edger) |
   !is.null(output.deseq) |
   !is.null(output.combined)){

  if(verbose){
    cat(paste("# Writing all DE tables: ", "\n", sep=""))
    if(!is.null(output.voom)){
      cat(paste("# Writing Voom DE tables to: ", output.voom, "\n",
                sep=""))
    }
    if(!is.null(output.edger)){
      cat(paste("# Writing EdgeR DE tables to: ", output.edger, "\n",
                sep=""))
    }
    if(!is.null(output.deseq)){
      cat(paste("# Writing DEseq2 DE tables to: ", output.deseq, "\n",
                sep=""))
    }
    if(!is.null(output.combined)){
      cat(paste("# Writing merged DE tables to: ", output.combined, "\n",
                sep=""))
    }
  }
  # invoke write.table.wrapper
  write.table.wrapper(summarised = summarised,
                      merged = list("merged" = merged,
                                    "deseq" = deseq.res,
                                    "edger" = edger.res,
                                    "voom" = voom.res),
                      voom.dir = output.voom,
                      edger.dir = output.edger,
                      deseq.dir = output.deseq,
                      merged.dir = output.combined)
}

if(verbose){
  cat(paste("# Finished DE analysis: ", "\n", sep=""))
}

# returns the merged results of the 3 comparisons.

if(ruv.correct == TRUE){
  return(list("merged" = merged,
              "deseq" = deseq.res,
              "edger" = edger.res,
              "voom" = voom.res,
              "summarised" = ruv.se[[1]]))
}

if(ruv.correct == FALSE){
  return(list("merged" = merged,
              "deseq" = deseq.res,
              "edger" = edger.res,
              "voom" = voom.res,
              "summarised" = summarised))
}

}

# helper function for annotating and writing tables:
# these are seperate commands for writing tables:
write.table.wrapper <- function(summarised = NULL,
                                merged = NULL,
                                voom.dir = NULL,
                                edger.dir = NULL,
                                deseq.dir = NULL,
                                merged.dir = NULL){

# format checks
  # check the merged is a list format and check the names! etc.
#########################################################
# obtain raw and normalised counts
raw <- assays(summarised)$counts
cpm.raw <- cpm(assays(summarised)$counts)
cpm.raw.log <- cpm(assays(summarised)$counts, log=TRUE)
cpm.mean <- table.means(cpm.raw, colData(summarised)$group)
cpm.mean.log <- table.means(cpm.raw.log, colData(summarised)$group)

# to normalise, push into seqexpressionset
se.new <- newSeqExpressionSet(assays(summarised)$counts,
                              phenoData = data.frame(colData(summarised)),
                              row.names=colnames(assays(summarised)$counts))
# normalise for library size
se.norm <- betweenLaneNormalization(se.new, which="upper")
norm.count <- normCounts(se.norm)
cpm.norm <- cpm(normCounts(se.norm))
cpm.norm.log <- cpm(normCounts(se.norm), log=TRUE)
cpm.norm.mean <- table.means(cpm.norm, colData(summarised)$group)
cpm.norm.log.mean <- table.means(cpm.norm.log, colData(summarised)$group)

# directories to print to
# check exists first!
if(!is.null(merged.dir)){
  utils::write.table(raw, file=paste(merged.dir, "raw_count.tsv", sep=""),
                     sep="\t")
  utils::write.table(cpm.raw, file=paste(merged.dir, "raw_cpm.tsv", sep=""),
                     sep="\t")
  utils::write.table(cpm.raw.log, file=paste(merged.dir, "raw_log_cpm.tsv",
                                             sep=""), sep="\t")
  utils::write.table(cpm.mean, file=paste(merged.dir, "raw_cpm_mean.tsv",
                                          sep=""), sep="\t", row.names=FALSE)
  utils::write.table(cpm.mean.log, file=paste(merged.dir,
                                              "raw_log_cpm_mean.tsv", sep=""),
                     sep="\t", row.names=FALSE)

  utils::write.table(norm.count, file=paste(merged.dir, "normalised_count.tsv",
                                            sep=""), sep="\t")
  utils::write.table(cpm.norm, file=paste(merged.dir, "normalised_cpm.tsv",
                                          sep=""), sep="\t")
  utils::write.table(cpm.norm.log, file=paste(merged.dir,
                                              "normalised_log_cpm.tsv", sep=""),
                     sep="\t")
  utils::write.table(cpm.norm.mean, file=paste(merged.dir,
                                               "normalised_cpm_mean.tsv",
                                               sep=""), sep="\t")
  utils::write.table(cpm.norm.log.mean, file=paste(merged.dir,
                                                   "raw_log_cpm_mean.tsv",
                                                   sep=""), sep="\t",
                     row.names=FALSE)

  # COMBINE WITH cpm.norm.mean [ normalised counts per million ] [ NOT LOGGED ]
  invisible(sapply(seq_len(length(merged$merged)), function(i)
    utils::write.table(data.frame(merge(cpm.norm.mean,
                                 merged$merged[[i]],
                                 by="ID"))[order(data.frame(merge(cpm.norm.mean,
                                                            merged$merged[[i]],
                                                            by="ID"))$rank.sum,
                                                 decreasing=FALSE), ],
                file=paste(merged.dir, names(merged$merged)[i],
                           "_combined_results.tsv", sep=""),
                sep="\t",
                row.names=FALSE)))
}

if(!is.null(deseq.dir)){
  # DESEQ:
  invisible(sapply(seq_len(length(merged$deseq$full.results)), function(i)
    utils::write.table(data.frame(merged$deseq$full.results[[i]]),
                file=paste(deseq.dir, names(merged$deseq$full.results[i]),
                           "_deseq_results.tsv", sep=""),
                sep="\t",
                row.names=TRUE)))
}
if(!is.null(voom.dir)){
  # VOOM
  invisible(sapply(seq_len(length(merged$voom$full.results)), function(i)
    utils::write.table(data.frame(merged$voom$full.results[[i]]),
                file=paste(voom.dir, names(merged$voom$full.results[i]),
                           "_voom_results.tsv", sep=""),
                sep="\t",
                row.names=TRUE)))
}

if(!is.null(edger.dir)){
  # EDGER
  invisible(sapply(seq_len(length(merged$edger$full.results)), function(i)
    utils::write.table(data.frame(merged$edger$full.results[[i]]),
                file=paste(edger.dir, names(merged$edger$full.results[i]),
                           "_edger_results.tsv", sep=""),
                sep="\t",
                row.names=TRUE)))
}
}
#########################################################

# helper function for renaming column names in a list
rename.colnames.list <- function(data.in = NULL,
                                 vector.search = NULL,
                                 vector.replace = NULL,
                                 which.list = NULL
                                 ){
  data.out <- data.in
  colnames(data.out) <- gsub(vector.search, vector.replace, colnames(data.out))
return(data.out)
}

# intersect function (for intersecting up to 3 vectors)
# reports intersect size and common set
# not used in V0.05
intersect.test <- function(x=NULL, y=NULL, z=NULL){
    if(!is.null(x) & !is.null(y)){
        cat(paste("# set1 = ", length(x), "\n"))
        cat(paste("# set2 = ", length(y), "\n"))
        common <- intersect(x, y)
    }
    if(!is.null(x) & !is.null(y) & !is.null(z)){
        cat(paste("# set3 = ", length(z), "\n"))
        common <- intersect(common, z)
    }
    cat(paste("# intersect = ", length(common), "\n"))
    return(common)
}

# function for obtaining non-de rows from the outputs of the XX$short.results
# from the DE analysis:
non.de.rows <- function(input.data, p.cut=0.5){
    data.out <- rownames(input.data[input.data$PValue > p.cut,])
    return(data.out)
}

# function for merging 3 DE tables in "short" format
# returns merged table with ranks for each
merge.results <- function(x, y, z,
name.x, name.y, name.z,
adjust.method="BH",
return="short"){
    x <- data.frame("ID"=rownames(x), x)
    y <- data.frame("ID"=rownames(y), y)
    z <- data.frame("ID"=rownames(z), z)

    merge.out <- merge(x, y, by="ID")
    merge.out <- merge(merge.out, z, by="ID")
    # adjust p-values
    merge.out$PValue.x <- stats::p.adjust(merge.out$PValue.x, adjust.method)
    merge.out$PValue.y <- stats::p.adjust(merge.out$PValue.y, adjust.method)
    merge.out$PValue <- stats::p.adjust(merge.out$PValue, adjust.method)
    # extract PValue columns for ranking.
    rank.this <- cbind(merge.out$PValue.x,
    merge.out$PValue.y,
    merge.out$PValue.z)
    # sum of ranks
    merge.out$sum <- rowSums(sapply(seq_len(ncol(rank.this)),
    function(i) rank(rank.this[,i])))
    # order results by sum of ranks
    merge.out <- merge.out[order(merge.out$sum, decreasing=FALSE),]

    if(return=="short"){
        merge.out <- data.frame(merge.out$ID,
        merge.out$AveExpr,
        merge.out$logFC.x,
        merge.out$PValue.x,
        merge.out$PValue.y,
        merge.out$PValue,
        rank(merge.out$PValue.x),
        rank(merge.out$PValue.y),
        rank(merge.out$PValue),
        merge.out$sum)
        colnames(merge.out) <- c("ID", "AveExpr", "LogFC",
        paste(name.x, ".adj.p", sep=""),
        paste(name.y, ".adj.p", sep=""),
        paste(name.z, ".adj.p", sep=""),
        paste(name.x, ".rank", sep=""),
        paste(name.y, ".rank", sep=""),
        paste(name.z, ".rank", sep=""),
        "rank.sum")

    }
    return(merge.out)
}

# determine means for a matrix and vector of common columns
table.means <- function(count.matrix, matrix.names){
    colnames(count.matrix) <- matrix.names
    mean_counts <- data.frame("ID"=rownames(count.matrix),
                        sapply(seq_len(length(unique(colnames(count.matrix)))),
                               function(i)
                                 rowMeans(count.matrix[,colnames(count.matrix)
                                    %in% unique(colnames(count.matrix))[i]])))
    colnames(mean_counts) <- c("ID", paste(unique(colnames(count.matrix)),
                                           ".mean", sep=""))
    return(mean_counts)
}

# function to automate contrast matrix creation of all pairs
build.contrast.matrix <- function(table.design, design){
    # base will already be levelled - so take the first level
    base.intercept <- levels(table.design$group)[1]
    # what is remaining in the groups (other pairs)
    remaining <- setdiff(levels(table.design$group), base.intercept)
    # if more than two comparisons
    if(length(levels(table.design$group)) > 2){
        # all possible combinations (of remaining pairs):
        all.com <- utils::combn(remaining, 2)
        remaining.com <- c(sapply(seq_len(ncol(all.com)), function(i)
        paste(all.com[1,i], all.com[2,i], sep="-")))
        # names for the comparisons
        real.com <- c(sapply(seq_len(length(remaining)), function(i)
                             paste(remaining[i], base.intercept, sep="-")))
        all.c <- paste(c(remaining, remaining.com), collapse=",")
        contrast.cmd <- paste("makeContrasts(", all.c, ",levels=design)",
                              sep="")
        # evaluate command string
        contrast.matrix <- eval(parse(text=contrast.cmd))
        colnames(contrast.matrix) <- c(real.com, remaining.com)
    }
    # if exactly two comparisons
    if(length(levels(table.design$group)) == 2){
        all.c <- levels(table.design$group)[2]
        contrast.cmd <- paste("makeContrasts(", all.c, ",levels=design)",
                              sep="")
        contrast.matrix <- eval(parse(text=contrast.cmd))
        colnames(contrast.matrix) <- paste(levels(table.design$group)[2],
        levels(table.design$group)[1], sep="-")
    }
    return(contrast.matrix)
}



# function to create 1) design matrix and 2) design table
# works with W_1 from RUV only
build.design <- function(se = NULL,
                         pairing = NULL,
                         intercept = NULL,
                         ruv = NULL){
    group <- factor(se$group)
    # set the intercept term/base coeffecient for comparisons.
    if(!is.null(intercept)){
        group <- stats::relevel(group, intercept)
    }
    design.table <- data.frame("file"=as.character(se$file),
                               "group"=group)
    # default is unpaired analysis
    if(pairing=="unpaired"){
        design <- stats::model.matrix(~group)
    }
    # paired, factor pairing
    if(pairing=="paired"){
        pairs <- factor(se$pairs)
        design.table <- data.frame(design.table,
        "pairs"=se$pairs)
        design <- stats::model.matrix(~pairs + group)
    }
    # for RUV n=1
    if(pairing=="unpaired" && ruv=="W_1"){
        W_1 <- as.numeric(se$W_1)
        design.table <- data.frame(design.table,
        "W_1"=se$W_1)
        design <- stats::model.matrix(~W_1 + group)
    }
    if(pairing=="paired" && ruv=="W_1"){
        W_1 <- as.numeric(se$W_1)
        design.table <- data.frame(design.table,
        "W_1"=se$W_1)
        design <- stats::model.matrix(~W_1 + pairs + group)
    }
    # format model.matrix
    rownames(design) <- colnames(se)
    colnames(design) <- gsub("group", "", colnames(design))
    colnames(design)[1] <- "Intercept"
return(list("design" = design, "table" = design.table))
}


# function for building DEseq compatabible contrast list from
# a contrast.matrix built using model.matrix()
DEseq.contrasts <- function(contrast.matrix, n){
    # n is the contrast of interest
    contrast.i <- contrast.matrix[,n]
    # check for positive and negative comparisons for each in the match:
    contrast.i.pos <- names(contrast.i[contrast.i==1])
    contrast.i.neg <- names(contrast.i[contrast.i==-1])
    contrast.i.name <- colnames(contrast.matrix)[n]
    # DEseq contrast lists are: nominator/denominator
    contrast.i.list <- list(c(contrast.i.pos), c(contrast.i.neg))
    return(list("contrast.list"=contrast.i.list, "contrast.name"=
                  contrast.i.name))
}


# annotation wrapper
annotate.ensembl <- function(data.in,
                             ids,
                             tx.db){

    short_ens <- gsub("\\..*", "", ids)
    data.in$genename <- mapIds(tx.db,
                              keys = short_ens,
                              column="GENENAME",
                              keytype="ENSEMBL",
                              multiVals="first")

    data.in$symbol <- mapIds(tx.db,
                            keys = short_ens,
                            column="SYMBOL",
                            keytype="ENSEMBL",
                            multiVals="first")
    data.in$kegg <- mapIds(tx.db,
                          keys = short_ens,
                          column="PATH",
                          keytype="ENSEMBL",
                          multiVals="first")
return(data.in)
}

# function for calling Voom analysis given se, design and contrast.matrix
voom.wrapper <- function(se, design, contrast.matrix){
    y <- DGEList(counts = assays(se)$counts,
                 group = se$group)
    y <- calcNormFactors(y)
    v <- voom(y, design = design)
    fit.voom <- lmFit(v, design = design)
    # fitting with user supplied contrast.matrix for each comparison
    all.fit <- lapply(seq_len(ncol(contrast.matrix)), function(i)
                      contrasts.fit(fit = fit.voom,
                              contrasts = t(data.frame(t(contrast.matrix[,i]),
                                  row.names = colnames(contrast.matrix)[i]))))

    all.fit <- lapply(seq_len(length(all.fit)), function(i)
    eBayes(all.fit[[i]]))

    all.results <- lapply(seq_len(length(all.fit)), function(i)
    topTable(all.fit[[i]], number=nrow(all.fit[[i]]$coefficients)))
    # short table for merging with other methods
    short.results <- lapply(seq_len(length(all.results)), function(i)
    data.frame("AveExpr"=all.results[[i]]$AveExpr,
    "logFC"=all.results[[i]]$logFC,
    "PValue"=all.results[[i]]$P.Value,
    row.names=rownames(all.results[[i]])))
    names(all.results) <- names(short.results) <- colnames(contrast.matrix)
    return(list("short.results"=short.results,
    "full.results"=all.results,
    "fitted"=all.fit,
    "contrasts"=contrast.matrix))
}

# function for calling DESeq2 analysis given se, design and contrast.matrix
edger.wrapper <- function(se, design, contrast.matrix){

    y <- DGEList(counts=assays(se)$counts,
    group=se$group)
    # normalise accross groups
    y <- calcNormFactors(y)
    # NB. there are different ways of estimating dispersion
    y.disp <- estimateDisp(y, design)
    fit.edger <- glmFit(y.disp, design)
    # fit each comparison from contrast.matrix of interest
    all.fit <- lapply(seq_len(ncol(contrast.matrix)), function(i)
    glmLRT(fit.edger, contrast=as.numeric(contrast.matrix[,i])))
    all.results <- lapply(seq_len(length(all.fit)), function(i)
    topTags(all.fit[[i]], n=nrow(all.fit[[i]]$coefficients)))
    # short table for merging with other methods
    short.results <- lapply(seq_len(length(all.results)), function(i)
    data.frame("AveExpr"=all.results[[i]]$table$logCPM,
    "logFC"=all.results[[i]]$table$logFC,
    "PValue"=all.results[[i]]$table$PValue,
    row.names=rownames(all.results[[i]]$table)))
    names(all.results) <- names(short.results) <- colnames(contrast.matrix)
    return(list("short.results"=short.results,
    "full.results"=all.results,
    "fitted"=all.fit,
    "contrasts"=contrast.matrix))
}

# function for calling DESeq2 analysis given se, design and contrast.matrix
deseq.wrapper <- function(se, design, contrast.matrix){

    dds <- DESeqDataSet(se, design = design)
    dds <- DESeq(dds)
    # extract all if the contrasts and return all results
    all.contrasts <- lapply(seq_len(ncol(contrast.matrix)),
    function(i) DEseq.contrasts(contrast.matrix, i))
    # conduct all tests
    # full.table of results
    all.results <- lapply(seq_len(length(all.contrasts)),
    function(i) results(dds, contrast=all.contrasts[[i]]$contrast.list)
    )
    names(all.results) <- sapply(seq_len(length(all.contrasts)), function(i)
    all.contrasts[[i]]$contrast.name)
    # short table for merging with other methods
    short.results <- lapply(seq_len(length(all.results)), function(i)
    data.frame("AveExpr"=log2(all.results[[i]]$baseMean),
    "logFC"=all.results[[i]]$log2FoldChange,
    "PValue"=all.results[[i]]$pvalue,
    row.names=rownames(all.results[[i]])))
    names(short.results) <- names(all.results)
    return(list("short.results"=short.results,
    "full.results"=all.results,
    "fitted"=dds,
    "contrasts"=contrast.matrix))
}

#########################################################
# Perform RUVr correction (OPTIONAL)
# - if correction with RUVr is used, update model

ruvr.correct <- function(se = NULL,
                         design = NULL,
                         pairing = "unpaired",
                         intercept = NULL,
                         plot.this = FALSE,
                         write.this = FALSE,
                         plot.dir = NULL,
                         verbose = FALSE){

    ruv <- newSeqExpressionSet(assays(se)$counts,
            phenoData = data.frame(colData(se)),
            row.names = colnames(assays(se)$counts))
    # normalise for library size (i.e. remove as an effect..)
    ruv.norm <- betweenLaneNormalization(ruv, which="upper")

    # fit GLM for estimation of W_1:
    ruv.y <- DGEList(counts=counts(ruv.norm), group=ruv.norm$group)
    ruv.y <- calcNormFactors(ruv.y, method="upperquartile")
    ruv.y <- estimateDisp(ruv.y, design)
    fit.ruv <- glmFit(ruv.y, design)
    res.ruv <- stats::residuals(fit.ruv, type="deviance")
    # >>>>this is a questionable parameter, using k=1
    # Assuming only a batch factor to remove
    ruv.se <- RUVr(ruv.norm, rownames(ruv.norm), k=1, res.ruv)
    # build design and contrast matrix for refitting model
    design <- build.design(se=ruv.se,
                           pairing=pairing,
                           intercept=intercept,
                           ruv="W_1")
    #contrast.matrix <- build.contrast.matrix(desig, design)

    # output RLE and PCA plots
    if(plot.this == TRUE){
        if(verbose){
            cat(paste("# Plotting RUV Dignostic plots. Will be located here: ",
                      plot.dir, "\n", sep=""))
        }

        invisible(diag.plots(se.in = ruv.se,
                             write = write.this,
                             plot.dir = plot.dir,
                             rle=plot.this,
                             pca=plot.this,
                             hclust=plot.this,
                             density=plot.this,
                             boxplot=plot.this,
                             name="YES_corrected_YES_normalised"))
        # plot model residuals
        invisible(diag.plots(se.in = ruv.se,
                             write = write.this,
                             plot.dir = plot.dir,
                             residuals = plot.this,
                             name = "RUVr_residuals"))
    }

    # update colData in the SE to include the W_1
    colData(se) <- S4Vectors::DataFrame(design$table)
    colnames(se) <- colData(se)$file

    return(list(se, design))
}

