#' QC/diagnostic plotting
#
#' @description Wrappers for a series of plots to be used as diagnostics in
#' RNA-seq analyses. Currently 10 plots are possible using this function: 1)
#' Mapped reads, 2) Relative Log Expression (RLE), 3) Principle Component
#' Analyis (PCA), 4) Residuals from a batch correction model, e.g. RUVseq, 5)
#' Hierarchical clustering, 6) Densitiy distributions, 7) Boxplots, 8) MA plots,
#'  9) Volcano Plots and 10) P-value distribution plots. Plots 1 to 6 utilise a
#'  "SeqExpressionSet" object for extracting information to plot. Plots 8-10
#'  utilised a simple list class, containing all the data.frames of each
#'  comparison performed. See descriptions of each in the parameter options
#'  below and for format specification. See vignette for more information and
#'  examples.
#
#' @param se_in A "SeqExpressionSet" object or "RangedSummarizedExperiment"
#' generated using "buildSummarized()". If the input is a "SeqExpressionSet",
#' ensure that it included groups to be analysed. E.g. accessible as
#' "se_in$group. Groupings are used to automate colouring of samples in
#' unsupervised analyses. Default = NULL
#' @param merged_in A data.frame that contains the merged results which are
#' included in the outputs from multi_de_pairs(). These contain the ouputs from
#' the pair-wise comparisons which allows plotting of MA, Volcano and p-value
#' distributions. Where the outputs of multi_de_pairs() are to be used as inputs
#'  into diag_plots(), use multi_de_pairs()$merged as inputs. See example below.
#'   Default = NULL
#' @param write Write the results to a pdf file? Options: TRUE, FALSE. This is
#' to be used together with "plot_dir" and "write" parameters (below). Will
#' report an error and halt if is TRUE and "plot_dir" and "write" are NULL.
#' Default = FALSE
#' @param plot_dir If "write" is TRUE, where to write the files to? The
#' directory must already exist. E.g. "/path/to/my/pretty/plots/". Default =
#' NULL
#' @param legend Include legend in plots? Legend is based on group data in
#' se_in. Options: TRUE, FALSE. Default = FALSE
#' @param label Include point labels in plots? Points are based on ID column
#' from merged_in. Options: TRUE, FALSE. Default = FALSE
#' @param name If "write" is TRUE, what to name the plot? The file name will
#' always be preceded with "QC_" and end in ".pdf". E.g.
#' name="very_pretty_plots" will produce a file named "QC_very_pretty_plots.pdf"
#'  in "/path/to/my/pretty/plots/". Default = NULL
#' @param mapped_reads Plot mapped reads per sample as a barchart. Requires
#' se_in to be a "SeqExpressionSet" and utilise "group" meta-data for colouring.
#'  Options: TRUE, FALSE. Default = FALSE
#' @param rle Plot Relative Log Expressio (RLE) of samples for assessment of
#' sample quality. See ?plotRLE for further details. Requires se_in to be a
#' "SeqExpressionSet"and utilise "group" meta-data for colouring. Options: TRUE,
#' FALSE. Default = FALSE
#' @param pca Perform unsupervised Principle Component Analysis (PCA) and plot
#' results. By default performs Singular Value Decomposition. Requires se_in to
#' be a "SeqExpressionSet" and utilise "group" meta-data for colouring. Options:
#'  TRUE, FALSE. Default = FALSE
#' @param residuals If RUV-seq has been applied to dataset, plot the residuals
#' identified in the model. Only works for one set of residuals. Data is also
#' accessible using pData(se_in)$W_1. Requires se_in to be a "SeqExpressionSet"
#' and utilise "group" meta-data for colouring. Options: TRUE, FALSE. Default =
#' FALSE
#' @param hclust Performs unsupervised hierarchical clustering of samples.
#' Colours sample below plot according to group and numbered by inputs. Requires
#'  se_in to be a "SeqExpressionSet" and utilise "group" meta-data for
#'  colouring. Options: TRUE, FALSE. Default = FALSE
#' @param density Plot density distributions of log2(count-per-million). Will
#' automatically extract normalised counts over non-normalised counts is
#' available in "SeqExpressionSet". Requires se_in to be a "SeqExpressionSet"
#' and utilise "group" meta-data for colouring. Options: TRUE, FALSE. Default =
#'  FALSE
#' @param boxplot Boxplot of density distributions of log2(count-per-million).
#' Will automatically extract normalised counts over non-normalised counts is
#' available in "SeqExpressionSet". Requires se_in to be a "SeqExpressionSet"
#' and utilise "group" meta-data for colouring. Options: TRUE, FALSE. Default =
#' FALSE
#' @param ma Plot Mean versus. Log2 Fold-Change of comparison. Requires a
#' data.frame as input to "merged_in" with the following column names "ID",
#' "AvExpr", "Log2FC" and "Adj_PVal".The data frame should be sorted, as the top
#'  10 in the table are also  plotted. Options: TRUE,  FALSE. Default = FALSE
#' @param volcano Volcano plot of Log2 Fold-Change and significance of
#' comparison. Requires a data.frame as input to "merged_in" with the following
#' column names "ID", "AvExpr", "Log2FC" and "Adj_PVal". The data frame should
#' be sorted, as the top 10 in the table are also plotted. Options: TRUE,
#' FALSE. Default = FALSE
#' @param p_dist P-value distribution plot. Requires a data.frame as input to
#' "merged_in" with the following column names "ID", "AvExpr", "Log2FC" and
#' "Adj_PVal". The data frame should be sorted, as the top 10 in the table are
#' also plotted. Options: TRUE,  FALSE. Default = FALSE
#'
#' @examples
#' ## Load the example data set and attach
#' ## The example below will display a PCA plot before normalisation
#' data("airway")
#' ## Name the groups of the data.
#' colData(airway)$group <- colData(airway)$dex
#' ## Identify the file locations
#' colData(airway)$file <- rownames(colData(airway))
#' ## Filter low count data:
#' airway_filter <- buildSummarized(summarized = airway,
#'                                  filter = TRUE)
#' ## for illustration, use random sample of 1000 transcripts
#' set.seed(1234)
#' airway_filter <- sample(airway_filter, 1000)
#' ## The following is example code to perform a PCA plot
#' ## see vignette for more details of displaying each plot
#' ## diag_plots(se_in = airway_filter,
#' ##            name = "airway example data",
#' ##            pca = TRUE)
#'               
#' @return Returns pretty plots
#'
#' @export diag_plots
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pcaMethods prep pca R2cum scores
#' @importFrom DESeq2 counts results
#' @importFrom EDASeq normCounts newSeqExpressionSet plotRLE
#' @importFrom edgeR cpm
#' @importFrom SummarizedExperiment assays colData
#' @importFrom Biobase pData
#' @importFrom dendextend colored_bars set %>%
#' @importFrom stats as.dendrogram
#' @import airway

diag_plots <- function(se_in = NULL,
                       merged_in = NULL,
                       write = FALSE,
                       plot_dir = NULL,
                       legend = TRUE,
                       label = TRUE,
                       name = NULL,
                       mapped_reads = FALSE,
                       rle = FALSE,
                       pca = FALSE,
                       residuals = FALSE,
                       hclust = FALSE,
                       density = FALSE,
                       boxplot = FALSE,
                       ma = FALSE,
                       volcano = FALSE,
                       p_dist = FALSE
                       ){

####///---- check inputs ----\\\###

# format checks of input data types:
if(!is.null(se_in) && se_in@class == "RangedSummarizedExperiment"){
    se_in <- newSeqExpressionSet(assays(se_in)$counts,
                                phenoData = data.frame(colData(se_in)),
                                row.names = colnames(assays(se_in)$counts))
    }

if(!is.null(se_in) && se_in@class != "SeqExpressionSet")
  stop(paste("summarized file provided is not a SeqExpressionSet.", "\n",
              "Please produce a SeqExpressionSet.", "\n", sep=""))

if(!is.null(merged_in) && !inherits(merged_in, "list"))
    stop(paste("merged_in is not a list. If you want to plot with one comparison
                only, put the single dataframe into a list as follows:\n
                my_list <- list(\"name_of_comparison\" = merged_in)\n", sep=""))

# check that all are data.frames and all contain variables of interest:
df_check <- sapply(seq_len(length(merged_in)), function(i)
                   is.data.frame(merged_in[[i]]))
if(sum(df_check==TRUE) != length(merged_in))
    stop(paste("merged_in contains", sum(df_check==FALSE), "slots that are not
               \"data.frame\" objects.\n"))

df_col_check <- sapply(seq_len(length(merged_in)), function(i)
                       identical(c("ID", "AveExpr", "Adj_PVal") %in%
                       colnames(merged_in[[i]]), c(TRUE, TRUE, TRUE)))
if(sum(df_col_check==TRUE) != length(merged_in))
  stop(paste("merged_in contains", sum(df_col_check==FALSE), "slots that are not
             \"data.frame\" objects with column names containing \"ID\",
             \"AveExpr\", and \"Adj_PVal\".\n"))

# if writing, make sure other fields are not empty:
if(write==TRUE && (is.null(plot_dir) | is.null(name)))
  stop(paste("When write = TRUE, a path to the plot_dir and a name for the file
              must be provided\n"))


# check RUV-residuals exist:
if(residuals==TRUE && is.null(pData(se_in)$W_1)){
    warning(paste("Cannot plot residuals when pData(se_in)$W_1 is empty\n"))
  # this will continue the other plots
  residuals <- FALSE
}
####///---- FINISH CHECK of inputs ----\\\###

  # establish colours - "Set3" is up to 12 distinct colours
  # must have a minimum of 3, this set-up as follows:
  colors <- brewer.pal(max(length(unique(se_in$group)), 3), "Set3")

  # remove the QC_?
  if(write==TRUE){
    grDevices::pdf(file=paste(plot_dir, "QC_", name, ".pdf", sep=""))
  }

  if(mapped_reads==TRUE){
    read_counts <- colSums(counts(se_in))
    #par(mar=c(15,3,2,2)+1)
    graphics::barplot(read_counts,
            col=colors[se_in$group],
            #names_arg=se_in$file,
            names.arg=seq_len(length(se_in$file)),
            las=1,
            ylab="mapped reads",
            cex.names=0.5,
            cex=0.5)
    if(legend == TRUE){
      legend("topright",
             c(as.character(unique(se_in$group))),
             col=colors,
             pch = c(rep(19, length(unique(se_in$group)))),
             title = "SAMPLE GROUPS", inset = .02, cex=0.5)
    }
  }

  if(rle==TRUE){
    plotRLE(se_in, outline=FALSE,
            ylim=c(-2, 2),
            xlab = "file number (see legend)",
            ylab = "Relative Log Expression",
            col=colors[se_in$group],
            names=seq_len(length(se_in$file)),
            las=1, cex.axis=1, main=paste("RLE-", name))
    if(legend == TRUE){
      legend("topright",
             c(as.character(unique(se_in$group)),
               paste(seq_len(length(se_in$file)), se_in$file, sep="-")),
             col=c(colors[1:length(unique(se_in$group))],
                   rep("black", length(se_in$file))),
             pch = c(rep(19, length(unique(se_in$group))),
                     rep(0, length(se_in$file))),
             title = "SAMPLE GROUPS", inset = .02, cex=0.5)
    }
  }

  if(pca==TRUE){
    # this will plot 3 PCAs with and without labels
    plot_PCA_wrapper(se_in=se_in,
                     title=name,
                     comp1=1, comp2=2,
                     legend=legend,
                     label=label,
                     colors=colors)
  }

  if(hclust==TRUE){
    plot_hclust_wrapper(se_in=se_in,
                        title=name,
                        colors=colors,
                        name=name,
                        legend=legend)
  }

  if(density==TRUE){
    plot_density_wrapper(se_in=se_in,
                         title=name,
                         colors=colors,
                         legend=legend)
  }

  if(boxplot==TRUE){
    boxplot_wrapper(se_in=se_in,
                    title=name,
                    colors=colors,
                    legend=legend)
  }

  if(residuals==TRUE){
    # RUV residuals from GLM
    ruv_res <- pData(se_in)$W_1
    graphics::barplot(ruv_res,
            col=colors[se_in$group],
            names.arg=seq_len(length(se_in$file)),
            xlab = "file number (see legend)",
            las=1, ylab="residuals", cex.names=1)

    if(legend == TRUE){
      legend("topright",
             c(as.character(unique(se_in$group)),
               paste(seq_len(length(se_in$file)), se_in$file, sep="-")),
             col=c(colors[1:length(unique(se_in$group))],
                   rep("black", length(se_in$file))),
             pch = c(rep(19, length(unique(se_in$group))),
                     rep(0, length(se_in$file))),
             title = "SAMPLE GROUPS", inset = .02, cex=0.5)
    }
  }

  if(ma==TRUE){
    # will be a list of all the pairwise comparisons
    # could have option here, for if not a list...
    sapply(seq_len(length(merged_in)), function(i)
           plot_ma_wrapper(merged_in[[i]],
                           names(merged_in[i]),
                           label=label,
                           legend=legend))
  }

  if(volcano==TRUE){
    sapply(seq_len(length(merged_in)), function(i)
          plot_volcano_wrapper(merged_in[[i]],
                               names(merged_in[i]),
                               label=label,
                               legend=legend))
  }

  if(p_dist==TRUE){
    # distribution of p-values
    sapply(seq_len(length(merged_in)), function(i)
      barplot_pval_wrapper(merged_in[[i]],
                           names(merged_in[i])))
  }

  if(write==TRUE){
    grDevices::dev.off()
  }
}

# function to check if data has been normalised
# will take normalised data if it is available
check_normalise <- function(se_in=NULL){
  is_normalised <- TRUE
  # if it hasn't been normalised, all data will be NAs
  if((sum(is.na(normCounts(se_in))) == (nrow(se_in)*ncol(se_in)))==TRUE){
    is_normalised <- FALSE
  }
  if(is_normalised == TRUE){
    data_in <- cpm(normCounts(se_in), log=TRUE)
  }
  if(is_normalised == FALSE){
    data_in <- cpm(counts(se_in), log=TRUE)
  }
return(data_in)
}

barplot_pval_wrapper <- function(merged_in=NULL,
                                 names_merged_in=NULL){
  # create 10 bins from 0-1
  bin_counts <- sapply(seq_len(10), function(i)
    nrow(merged_in[merged_in$Adj_PVal <= i/10 &
                     merged_in$Adj_PVal > (i-1)/10,]))
  # plot bin counts
  graphics::barplot(bin_counts,
          las=1,
          ylab="frequency",
          names.arg=c(seq_len(10))/10,
          xlab="adj_p",
          main=names_merged_in)
}

plot_volcano_wrapper <- function(merged_in=NULL,
                                 names_merged_in=NULL,
                                 legend = TRUE,
                                 label = TRUE){

  # set of transcripts to plot at different p_cuts
  top10 <- merged_in[1:10,]
  p_01 <- merged_in[merged_in$Adj_PVal <= 0.01,]
  p_05 <- merged_in[merged_in$Adj_PVal > 0.01 & merged_in$Adj_PVal <= 0.05 ,]
  p_rest <- merged_in[merged_in$Adj_PVal > 0.05,]

  # from exprs
  graphics::plot(p_rest$LogFC,
       -log(p_rest$Adj_PVal, 10),
       main=names_merged_in,
       ylab="-log10(adj_p)",
       xlab="Log2FC",
       xlim=range(c(p_rest$LogFC, p_01$LogFC, p_05$LogFC, top10$LogFC)),
       ylim=range(c(-log(p_rest$Adj_PVal, 10),
                    -log(p_01$Adj_PVal, 10),
                    -log(p_05$Adj_PVal, 10),
                    -log(top10$Adj_PVal, 10))),
       cex=0.4, pch=16
  )
  graphics::points(top10$LogFC,
         -log(top10$Adj_PVal, 10),
         col="red",
         cex=1.5,
         pch=16)
  graphics::points(p_05$LogFC,
         -log(p_05$Adj_PVal, 10),
         col="lightgreen",
         cex=0.7,
         pch=16)
  graphics::points(p_01$LogFC,
         -log(p_01$Adj_PVal, 10),
         col="blue",
         cex=1,
         pch=16)
  if(legend == TRUE){
    legend("topright", c("top10", "p <= 0.01", "p <= 0.05"),
           col=c("red", "blue", "lightgreen"),
           pch = c(rep(16, 3)),
           title = "p-val cutoffs", inset = .02, cex=0.7)
  }
  if(label == TRUE){
    # label top 10 points
    graphics::text(top10$LogFC, -log(top10$Adj_PVal, 10), top10$ID, cex=0.7,
                   pos=4, col="black")
  }
}

plot_ma_wrapper <- function(merged_in = NULL,
                            names_merged_in = NULL,
                            legend = TRUE,
                            label = TRUE){

  # set of transcripts to plot at different p_cuts
  top10 <- merged_in[1:10,]
  p_01 <- merged_in[merged_in$Adj_PVal <= 0.01,]
  p_05 <- merged_in[merged_in$Adj_PVal > 0.01 & merged_in$Adj_PVal <= 0.05 ,]
  p_rest <- merged_in[merged_in$Adj_PVal > 0.05,]

  # from exprs
  graphics::plot(p_rest$AveExpr,
       p_rest$LogFC,
       main=names_merged_in,
       ylab="Log2FC",
       xlab="Average Expression",
       xlim=range(c(p_rest$AveExpr, p_01$AveExpr, p_05$AveExpr, top10$AveExpr)),
       ylim=range(c(p_rest$LogFC, p_01$LogFC, p_05$LogFC, top10$LogFC)),
       cex=0.4, pch=16
  )
  graphics::points(top10$AveExpr,
         top10$LogFC,
         col="red",
         cex=1.5,
         pch=16)
  graphics::points(p_05$AveExpr,
         p_05$LogFC,
         col="lightgreen",
         cex=0.7,
         pch=16)
  graphics::points(p_01$AveExpr,
         p_01$LogFC,
         col="blue",
         cex=1,
         pch=16)

  #abline(h=0, col="lightgrey")

  if(legend==TRUE){
    legend("topright", c("top10", "p <= 0.01", "p <= 0.05"),
           col=c("red", "blue", "lightgreen"),
           pch = c(rep(16, 3)),
           title = "p-val cutoffs", inset = .02, cex=0.7)

  }
  if(label==TRUE){
    # label top 10 points
    graphics::text(top10$AveExpr, top10$LogFC, top10$ID, cex=0.7, pos=4,
                   col="black")
  }
}

boxplot_wrapper <- function(se_in=NULL,
                            title=NA,
                            colors=NA,
                            legend=TRUE){

  data_in <- check_normalise(se_in=se_in)

  graphics::boxplot(data_in, col=colors[se_in$group],
          names=c(seq_len(ncol(data_in))), las=2,
          ylab="log(cpm)",
          xlab="Sample number (see legend)",
          cex=0.1,
          main=paste("Boxplot - ", title, sep=""))

  if(legend == TRUE){
    legend("topright",
           c(as.character(unique(se_in$group)),
             paste(seq_len(length(se_in$file)), se_in$file, sep="-")),
           col=c(colors[1:length(unique(se_in$group))],
                 rep("black", length(se_in$file))),
           pch = c(rep(19, length(unique(se_in$group))),
                   rep(0, length(se_in$file))),
           title = "SAMPLE GROUPS", inset = .02, cex=0.5)
  }
}

plot_density_wrapper <- function(se_in=NULL,
                                 title=NA,
                                 colors=NA,
                                 legend=TRUE){

  data_in <- check_normalise(se_in=se_in)
  # colours
  colors_density <- colors[se_in$group]
  # determine densities
  # everything:
  dens <- stats::density(data_in)

  all_dens <- lapply(seq_len(ncol(data_in)), function(i)
    stats::density(data_in[,i]))


  y_range <- range(sapply(seq_len(length(all_dens)), function(i)
                   all_dens[[i]]$y))
  # to include everything in the range:
  y_range <- range(c(y_range, dens$y))
  x_range <- range(sapply(seq_len(length(all_dens)), function(i)
                   all_dens[[i]]$x))

  # initialise plot
  graphics::plot(all_dens[[1]], xlim = x_range, ylim = y_range,
       main=paste("Density - ", title, sep=""),
       col = colors_density[1],
       xlab="log(cpm)")
  # add samples
  sapply(2:length(all_dens), function(i)
    graphics::lines(all_dens[[i]], col = colors_density[i])
  )
  # add total density of all samples
  graphics::lines(dens, col ="black", lwd=2)
  # add legend
  if(legend==TRUE){
    legend("topright", c(as.character(unique(se_in$group)), "ALL SAMPLES"),
           col=c(colors, "black"),
           pch = c(rep(19, length(unique(se_in$group)))),
           title = "SAMPLE GROUPS", inset = .02, cex=0.5)
  }
}

plot_hclust_wrapper <- function(se_in=NULL,
                                title=NA,
                                colors=NA,
                                name=name,
                                legend=TRUE
                                ){

  data_in <- check_normalise(se_in=se_in)
  data_in <- t(data_in)
  # relabel hclust by sample number
  rownames(data_in) <- c(seq_len(nrow(data_in)))
  dd <- stats::dist(scale(data_in), method = "euclidean")
  hc <- stats::hclust(dd, method = "ward.D2")
  # labelling of nodes, by samples
  colors_hclust <- colors[se_in$group]
  # adding coloured bars
  the_bars <- colors_hclust
  colors_hclust <- sort(colors_hclust)[hc$order]

  # GENERATE DENDROGRAM
  dend <- data_in %>% scale %>% stats::dist(method = "euclidean") %>%
    stats::hclust(method = "ward.D2") %>% as.dendrogram
  dend %>% set("labels_col", value=c(colors_hclust))
  dend %>% graphics::plot(main = name,
                          sub="euclidean + ward(see legend for sample numbers)")
  colored_bars(colors = the_bars, dend = dend, sort_by_labels_order = TRUE)

  if(legend == TRUE){
    legend("topright",
           c(as.character(unique(se_in$group)),
             paste(seq_len(length(se_in$file)), se_in$file, sep="-")),
           col=c(colors[1:length(unique(se_in$group))],
                 rep("black", length(se_in$file))),
           pch = c(rep(19, length(unique(se_in$group))),
                   rep(0, length(se_in$file))),
           title = "SAMPLE GROUPS", inset = .02, cex=0.5)
  }

}

plot_PCA_wrapper <- function(se_in=NULL,
                             title=NA,
                             comp1=1,
                             comp2=2,
                             legend=TRUE,
                             label=TRUE,
                             colors=NA){

  data_in <- check_normalise(se_in=se_in)
  # data transformations
  md <- prep(t(data_in), scale = "none", centre = FALSE)
  pc <- pca(md, method="svd", center=FALSE, nPcs=ncol(data_in))
  var_3 <- R2cum(pc)[3] # accumulated variance
  pc_1 <- round(pc@R2[comp1]*100, 2)
  pc_2 <- round(pc@R2[comp2]*100, 2)

  pc_scores <- as.data.frame(scores(pc))
  pc_scores <- data.frame(pc_scores, "group"=se_in$group, "file"=se_in$file)

  # initialise plot:
  # -2 is remove the group and file name variable.
  graphics::plot(1, type="n", xlim=c(min(pc_scores[comp1])-5,
                           max(pc_scores[comp1])+5),
                    ylim=c(min(pc_scores[comp2])-5,
                           max(pc_scores[comp2])+5),
       axes=TRUE,
       xlab=paste("PC", comp1, " - ", pc_1, "%", sep=""),
       ylab=paste("PC", comp2, " - ", pc_2, "%", sep=""),
       main=paste("PCA - ", title, sep="")
       )

  # add labels:
  graphics::abline(h = 0, v = 0, col = "gray", lty = 2)
  # add legend:
  if(legend==TRUE){
    legend("bottomright", c(as.character(unique(pc_scores$group))),
                            col=colors,
                          pch = c(rep(19, length(unique(pc_scores$group)))),
           title = "SAMPLE GROUPS", inset = .02, cex=0.5)
  }
  if(label==TRUE){
    graphics::text(pc_scores[,comp1], pc_scores[,comp2], pc_scores$file,
                   cex=0.5, pos=3, col="black")
  }
  graphics::points(pc_scores[,comp1], pc_scores[,comp2], cex = 1,
                   col = colors[se_in$group], pch=19)
}

