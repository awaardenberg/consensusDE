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
#' @param se.in A "SeqExpressionSet" object or "RangedSummarizedExperiment"
#' generated using "read.summarised()". If the input is a "SeqExpressionSet",
#' ensure that it included groups to be analysed. E.g. accessible as
#' "se.in$group. Groupings are used to automate colouring of samples in
#' unsupervised analyses. Default = NULL
#' @param merged.in A data.frame that contains the merged results which are
#' included in the outputs from multi.de.pairs(). These contain the ouputs from
#' the pair-wise comparisons which allows plotting of MA, Volcano and p-value
#' distributions. Where the outputs of multi.de.pairs() are to be used as inputs
#'  into diag.plots(), use multi.de.pairs()$merged as inputs. See example below.
#'   Default = NULL
#' @param write Write the results to a pdf file? Options: TRUE, FALSE. This is
#' to be used together with "plot.dir" and "write" parameters (below). Will
#' report an error and halt if is TRUE and "plot.dir" and "write" are NULL.
#' Default = FALSE
#' @param plot.dir If "write" is TRUE, where to write the files to? The
#' directory must already exist. E.g. "/path/to/my/pretty/plots/". Default =
#' NULL
#' @param legend Include legend in plots? Legend is based on group data in
#' se.in. Options: TRUE, FALSE. Default = FALSE
#' @param label Include point labels in plots? Points are based on ID column
#' from merged.in. Options: TRUE, FALSE. Default = FALSE
#' @param name If "write" is TRUE, what to name the plot? The file name will
#' always be preceded with "QC_" and end in ".pdf". E.g.
#' name="very_pretty_plots" will produce a file named "QC_very_pretty_plots.pdf"
#'  in "/path/to/my/pretty/plots/".Default = NULL
#' @param mapped.reads Plot mapped reads per sample as a barchart. Requires
#' se.in to be a "SeqExpressionSet" and utilise "group" meta-data for colouring.
#'  Options: TRUE, FALSE. Default = FALSE
#' @param rle Plot Relative Log Expressio (RLE) of samples for assessment of
#' sample quality. See ?plotRLE for further details. Requires se.in to be a
#' "SeqExpressionSet"and utilise "group" meta-data for colouring. Options: TRUE,
#' FALSE. Default = FALSE
#' @param pca Perform unsupervised Principle Component Analysis (PCA) and plot
#' results. By default performs Singular Value Decomposition. Requires se.in to
#' be a "SeqExpressionSet" and utilise "group" meta-data for colouring. Options:
#'  TRUE, FALSE. Default = FALSE
#' @param residuals If RUV-seq has been applied to dataset, plot the residuals
#' identified in the model. Only works for one set of residuals. Data is also
#' accessible using pData(se.in)$W_1. Requires se.in to be a "SeqExpressionSet"
#' and utilise "group" meta-data for colouring. Options: TRUE, FALSE. Default =
#' FALSE
#' @param hclust Performs unsupervised hierarchical clustering of samples.
#' Colours sample below plot according to group and numbered by inputs. Requires
#'  se.in to be a "SeqExpressionSet" and utilise "group" meta-data for
#'  colouring. Options: TRUE, FALSE. Default = FALSE
#' @param density Plot density distributions of log2(count-per-million). Will
#' automatically extract normalised counts over non-normalised counts is
#' available in "SeqExpressionSet". Requires se.in to be a "SeqExpressionSet"
#' and utilise "group" meta-data for colouring. Options: TRUE, FALSE. Default =
#'  FALSE
#' @param boxplot Boxplot of density distributions of log2(count-per-million).
#' Will automatically extract normalised counts over non-normalised counts is
#' available in "SeqExpressionSet". Requires se.in to be a "SeqExpressionSet"
#' and utilise "group" meta-data for colouring. Options: TRUE, FALSE. Default =
#' FALSE
#' @param ma Plot Mean versus. Log2 Fold-Change of comparison. Requires a
#' data.frame as input to "merged.in" with the following column names "ID",
#' "AvExpr", "Log2FC" and "Adj.PVal".The data frame should be sorted, as the top
#'  10 in the table are also  plotted. Options: TRUE,  FALSE. Default = FALSE
#' @param volcano Volcano plot of Log2 Fold-Change and significance of
#' comparison. Requires a data.frame as input to "merged.in" with the following
#' column names "ID", "AvExpr", "Log2FC" and "Adj.PVal". The data frame should
#' be sorted, as the top 10 in the table are also plotted. Options: TRUE,
#' FALSE. Default = FALSE
#' @param p.dist P-value distribution plot. Requires a data.frame as input to
#' "merged.in" with the following column names "ID", "AvExpr", "Log2FC" and
#' "Adj.PVal". The data frame should be sorted, as the top 10 in the table are
#' also plotted. Options: TRUE,  FALSE. Default = FALSE
#'
#' @examples
#' ## Load the example data set and attach
#' library(airway)
#' data("airway")
#' ## Name the groups of the data.
#' colData(airway)$group <- colData(airway)$dex
#' ## Identify the file locations
#' colData(airway)$file <- rownames(colData(airway))
#' ## Filter low count data:
#' airway.filter <- read.summarised(summarised = airway,
#'                                  filter = TRUE)
#' ## Below we will perform a PCA plot
#' ## see vignette for more details of displaying each plot
#'diag.plots(se.in = airway.filter,
#'           name = "airway example data",
#'           pca = TRUE)
#'
#' @return Returns pretty plots.
#'
#' @export diag.plots
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

diag.plots <- function(se.in = NULL,
                       merged.in = NULL,
                       write = FALSE,
                       plot.dir = NULL,
                       legend = TRUE,
                       label = TRUE,
                       name = NULL,
                       mapped.reads = FALSE,
                       rle = FALSE,
                       pca = FALSE,
                       residuals = FALSE,
                       hclust = FALSE,
                       density = FALSE,
                       boxplot = FALSE,
                       ma = FALSE,
                       volcano = FALSE,
                       p.dist = FALSE
                       ){

####///---- check inputs ----\\\###

# format checks of input data types:
if(!is.null(se.in) && se.in@class == "RangedSummarizedExperiment"){
    se.in <- newSeqExpressionSet(assays(se.in)$counts,
                                phenoData = data.frame(colData(se.in)),
                                row.names = colnames(assays(se.in)$counts))
    }

if(!is.null(se.in) && se.in@class != "SeqExpressionSet")
  stop(paste("Summarised file provided is not a SeqExpressionSet.", "\n",
              "Please produce a SeqExpressionSet.", "\n", sep=""))

if(!is.null(merged.in) && !inherits(merged.in, "list"))
    stop(paste("merged.in is not a list. If you want to plot with one comparison
                only, put the single dataframe into a list as follows:\n
                my.list <- list(\"name.of.comparison\" = merged.in)\n", sep=""))

# check that all are data.frames and all contain variables of interest:
df.check <- sapply(seq_len(length(merged.in)), function(i)
                   is.data.frame(merged.in[[i]]))
if(sum(df.check==TRUE) != length(merged.in))
    stop(paste("merged.in contains", sum(df.check==FALSE), "slots that are not
               \"data.frame\" objects.\n"))

df.col.check <- sapply(seq_len(length(merged.in)), function(i)
                       identical(c("ID", "AveExpr", "Adj.PVal") %in%
                       colnames(merged.in[[i]]), c(TRUE, TRUE, TRUE)))
if(sum(df.col.check==TRUE) != length(merged.in))
  stop(paste("merged.in contains", sum(df.col.check==FALSE), "slots that are not
             \"data.frame\" objects with column names containing \"ID\",
             \"AveExpr\", and \"Adj.PVal\".\n"))

# if writing, make sure other fields are not empty:
if(write==TRUE && (is.null(plot.dir) | is.null(name)))
  stop(paste("When write = TRUE, a path to the plot.dir and a name for the file
              must be provided\n"))


# check RUV-residuals exist:
if(residuals==TRUE && is.null(pData(se.in)$W_1)){
    warning(paste("Cannot plot residuals when pData(se.in)$W_1 is empty\n"))
  # this will continue the other plots
  residuals <- FALSE
}
####///---- FINISH CHECK of inputs ----\\\###

  # establish colours - "Set3" is up to 12 distinct colours
  # must have a minimum of 3, this set-up as follows:
  colors <- brewer.pal(max(length(unique(se.in$group)), 3), "Set3")

  # remove the QC_?
  if(write==TRUE){
    grDevices::pdf(file=paste(plot.dir, "QC_", name, ".pdf", sep=""))
  }

  if(mapped.reads==TRUE){
    read.counts <- colSums(counts(se.in))
    #par(mar=c(15,3,2,2)+1)
    graphics::barplot(read.counts,
            col=colors[se.in$group],
            #names.arg=se.in$file,
            names.arg=seq_len(length(se.in$file)),
            las=1,
            ylab="mapped reads",
            cex.names=0.5,
            cex=0.5)
    if(legend == TRUE){
      legend("topright",
             c(as.character(unique(se.in$group))),
             col=colors,
             pch = c(rep(19, length(unique(se.in$group)))),
             title = "SAMPLE GROUPS", inset = .02, cex=0.5)
    }
  }

  if(rle==TRUE){
    plotRLE(se.in, outline=FALSE,
            ylim=c(-2, 2),
            xlab = "file number (see legend)",
            ylab = "Relative Log Expression",
            col=colors[se.in$group],
            names=seq_len(length(se.in$file)),
            las=1, cex.axis=1, main=paste("RLE-", name))
    if(legend == TRUE){
      legend("topright",
             c(as.character(unique(se.in$group)),
               paste(seq_len(length(se.in$file)), se.in$file, sep="-")),
             col=c(colors[1:length(unique(se.in$group))],
                   rep("black", length(se.in$file))),
             pch = c(rep(19, length(unique(se.in$group))),
                     rep(0, length(se.in$file))),
             title = "SAMPLE GROUPS", inset = .02, cex=0.5)
    }
  }

  if(pca==TRUE){
    # this will plot 3 PCAs with and without labels
    plot.PCA.wrapper(se.in=se.in,
                     title=name,
                     comp1=1, comp2=2,
                     legend=legend,
                     label=label,
                     colors=colors)
  }

  if(hclust==TRUE){
    plot.hclust.wrapper(se.in=se.in,
                        title=name,
                        colors=colors,
                        name=name,
                        legend=legend)
  }

  if(density==TRUE){
    plot.density.wrapper(se.in=se.in,
                         title=name,
                         colors=colors,
                         legend=legend)
  }

  if(boxplot==TRUE){
    boxplot.wrapper(se.in=se.in,
                    title=name,
                    colors=colors,
                    legend=legend)
  }

  if(residuals==TRUE){
    # RUV residuals from GLM
    ruv.res <- pData(se.in)$W_1
    graphics::barplot(ruv.res,
            col=colors[se.in$group],
            names.arg=seq_len(length(se.in$file)),
            xlab = "file number (see legend)",
            las=1, ylab="residuals", cex.names=1)

    if(legend == TRUE){
      legend("topright",
             c(as.character(unique(se.in$group)),
               paste(seq_len(length(se.in$file)), se.in$file, sep="-")),
             col=c(colors[1:length(unique(se.in$group))],
                   rep("black", length(se.in$file))),
             pch = c(rep(19, length(unique(se.in$group))),
                     rep(0, length(se.in$file))),
             title = "SAMPLE GROUPS", inset = .02, cex=0.5)
    }
  }

  if(ma==TRUE){
    # will be a list of all the pairwise comparisons
    # could have option here, for if not a list...
    sapply(seq_len(length(merged.in)), function(i)
           plot.ma.wrapper(merged.in[[i]],
                           names(merged.in[i]),
                           label=label,
                           legend=legend))
  }

  if(volcano==TRUE){
    sapply(seq_len(length(merged.in)), function(i)
          plot.volcano.wrapper(merged.in[[i]],
                               names(merged.in[i]),
                               label=label,
                               legend=legend))
  }

  if(p.dist==TRUE){
    # distribution of p-values
    sapply(seq_len(length(merged.in)), function(i)
      barplot.pval.wrapper(merged.in[[i]],
                           names(merged.in[i])))
  }

  if(write==TRUE){
    grDevices::dev.off()
  }
}

# function to check if data has been normalised
# will take normalised data if it is available
check.normalise <- function(se.in=NULL){
  is.normalised <- TRUE
  # if it hasn't been normalised, all data will be NAs
  if((sum(is.na(normCounts(se.in))) == (nrow(se.in)*ncol(se.in)))==TRUE){
    is.normalised <- FALSE
  }
  if(is.normalised == TRUE){
    data.in <- cpm(normCounts(se.in), log=TRUE)
  }
  if(is.normalised == FALSE){
    data.in <- cpm(counts(se.in), log=TRUE)
  }
return(data.in)
}

barplot.pval.wrapper <- function(merged.in=NULL,
                                 names.merged.in=NULL){
  # create 10 bins from 0-1
  bin.counts <- sapply(seq_len(10), function(i)
    nrow(merged.in[merged.in$Adj.PVal <= i/10 &
                     merged.in$Adj.PVal > (i-1)/10,]))
  # plot bin counts
  graphics::barplot(bin.counts,
          las=1,
          ylab="frequency",
          names.arg=c(seq_len(10))/10,
          xlab="adj.p",
          main=names.merged.in)
}

plot.volcano.wrapper <- function(merged.in=NULL,
                                 names.merged.in=NULL,
                                 legend = TRUE,
                                 label = TRUE){

  # set of transcripts to plot at different p.cuts
  top10 <- merged.in[1:10,]
  p.01 <- merged.in[merged.in$Adj.PVal <= 0.01,]
  p.05 <- merged.in[merged.in$Adj.PVal > 0.01 & merged.in$Adj.PVal <= 0.05 ,]
  p.rest <- merged.in[merged.in$Adj.PVal > 0.05,]

  # from exprs
  graphics::plot(p.rest$LogFC,
       -log(p.rest$Adj.PVal, 10),
       main=names.merged.in,
       ylab="-log10(adj.p)",
       xlab="Log2FC",
       xlim=range(c(p.rest$LogFC, p.01$LogFC, p.05$LogFC, top10$LogFC)),
       ylim=range(c(-log(p.rest$Adj.PVal, 10),
                    -log(p.01$Adj.PVal, 10),
                    -log(p.05$Adj.PVal, 10),
                    -log(top10$Adj.PVal, 10))),
       cex=0.4, pch=16
  )
  graphics::points(top10$LogFC,
         -log(top10$Adj.PVal, 10),
         col="red",
         cex=1.5,
         pch=16)
  graphics::points(p.05$LogFC,
         -log(p.05$Adj.PVal, 10),
         col="lightgreen",
         cex=0.7,
         pch=16)
  graphics::points(p.01$LogFC,
         -log(p.01$Adj.PVal, 10),
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
    graphics::text(top10$LogFC, -log(top10$Adj.PVal, 10), top10$ID, cex=0.7,
                   pos=4, col="black")
  }
}

plot.ma.wrapper <- function(merged.in = NULL,
                            names.merged.in = NULL,
                            legend = TRUE,
                            label = TRUE){

  # set of transcripts to plot at different p.cuts
  top10 <- merged.in[1:10,]
  p.01 <- merged.in[merged.in$Adj.PVal <= 0.01,]
  p.05 <- merged.in[merged.in$Adj.PVal > 0.01 & merged.in$Adj.PVal <= 0.05 ,]
  p.rest <- merged.in[merged.in$Adj.PVal > 0.05,]

  # from exprs
  graphics::plot(p.rest$AveExpr,
       p.rest$LogFC,
       main=names.merged.in,
       ylab="Log2FC",
       xlab="Average Expression",
       xlim=range(c(p.rest$AveExpr, p.01$AveExpr, p.05$AveExpr, top10$AveExpr)),
       ylim=range(c(p.rest$LogFC, p.01$LogFC, p.05$LogFC, top10$LogFC)),
       cex=0.4, pch=16
  )
  graphics::points(top10$AveExpr,
         top10$LogFC,
         col="red",
         cex=1.5,
         pch=16)
  graphics::points(p.05$AveExpr,
         p.05$LogFC,
         col="lightgreen",
         cex=0.7,
         pch=16)
  graphics::points(p.01$AveExpr,
         p.01$LogFC,
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

boxplot.wrapper <- function(se.in=NULL,
                            title=NA,
                            colors=NA,
                            legend=TRUE){

  data.in <- check.normalise(se.in=se.in)

  graphics::boxplot(data.in, col=colors[se.in$group],
          names=c(seq_len(ncol(data.in))), las=2,
          ylab="log(cpm)",
          xlab="Sample number (see legend)",
          cex=0.1,
          main=paste("Boxplot - ", title, sep=""))

  if(legend == TRUE){
    legend("topright",
           c(as.character(unique(se.in$group)),
             paste(seq_len(length(se.in$file)), se.in$file, sep="-")),
           col=c(colors[1:length(unique(se.in$group))],
                 rep("black", length(se.in$file))),
           pch = c(rep(19, length(unique(se.in$group))),
                   rep(0, length(se.in$file))),
           title = "SAMPLE GROUPS", inset = .02, cex=0.5)
  }
}

plot.density.wrapper <- function(se.in=NULL,
                                 title=NA,
                                 colors=NA,
                                 legend=TRUE){

  data.in <- check.normalise(se.in=se.in)
  # colours
  colors.density <- colors[se.in$group]
  # determine densities
  # everything:
  dens <- stats::density(data.in)

  all.dens <- lapply(seq_len(ncol(data.in)), function(i)
    stats::density(data.in[,i]))


  y.range <- range(sapply(seq_len(length(all.dens)), function(i)
                   all.dens[[i]]$y))
  # to include everything in the range:
  y.range <- range(c(y.range, dens$y))
  x.range <- range(sapply(seq_len(length(all.dens)), function(i)
                   all.dens[[i]]$x))

  # initialise plot
  graphics::plot(all.dens[[1]], xlim = x.range, ylim = y.range,
       main=paste("Density - ", title, sep=""),
       col = colors.density[1],
       xlab="log(cpm)")
  # add samples
  sapply(2:length(all.dens), function(i)
    graphics::lines(all.dens[[i]], col = colors.density[i])
  )
  # add total density of all samples
  graphics::lines(dens, col ="black", lwd=2)
  # add legend
  if(legend==TRUE){
    legend("topright", c(as.character(unique(se.in$group)), "ALL SAMPLES"),
           col=c(colors, "black"),
           pch = c(rep(19, length(unique(se.in$group)))),
           title = "SAMPLE GROUPS", inset = .02, cex=0.5)
  }
}

plot.hclust.wrapper <- function(se.in=NULL,
                                title=NA,
                                colors=NA,
                                name=name,
                                legend=TRUE
                                ){

  data.in <- check.normalise(se.in=se.in)
  data.in <- t(data.in)
  # relabel hclust by sample number
  rownames(data.in) <- c(seq_len(nrow(data.in)))
  dd <- stats::dist(scale(data.in), method = "euclidean")
  hc <- stats::hclust(dd, method = "ward.D2")
  # labelling of nodes, by samples
  colors.hclust <- colors[se.in$group]
  # adding coloured bars
  the_bars <- colors.hclust
  colors.hclust <- sort(colors.hclust)[hc$order]

  # GENERATE DENDROGRAM
  dend <- data.in %>% scale %>% stats::dist(method = "euclidean") %>%
    stats::hclust(method = "ward.D2") %>% as.dendrogram
  dend %>% set("labels_col", value=c(colors.hclust))
  dend %>% graphics::plot(main = name,
                          sub="euclidean + ward(see legend for sample numbers)")
  colored_bars(colors = the_bars, dend = dend, sort_by_labels_order = TRUE)

  if(legend == TRUE){
    legend("topright",
           c(as.character(unique(se.in$group)),
             paste(seq_len(length(se.in$file)), se.in$file, sep="-")),
           col=c(colors[1:length(unique(se.in$group))],
                 rep("black", length(se.in$file))),
           pch = c(rep(19, length(unique(se.in$group))),
                   rep(0, length(se.in$file))),
           title = "SAMPLE GROUPS", inset = .02, cex=0.5)
  }

}

plot.PCA.wrapper <- function(se.in=NULL,
                             title=NA,
                             comp1=1,
                             comp2=2,
                             legend=TRUE,
                             label=TRUE,
                             colors=NA){

  data.in <- check.normalise(se.in=se.in)
  # data transformations
  md <- prep(t(data.in), scale = "none", centre = FALSE)
  pc <- pca(md, method="svd", center=FALSE, nPcs=ncol(data.in))
  var_3 <- R2cum(pc)[3] # accumulated variance
  pc.1 <- round(pc@R2[comp1]*100, 2)
  pc.2 <- round(pc@R2[comp2]*100, 2)

  pc.scores <- as.data.frame(scores(pc))
  pc.scores <- data.frame(pc.scores, "group"=se.in$group, "file"=se.in$file)

  # initialise plot:
  # -2 is remove the group and file name variable.
  graphics::plot(1, type="n", xlim=c(min(pc.scores[comp1])-5,
                           max(pc.scores[comp1])+5),
                    ylim=c(min(pc.scores[comp2])-5,
                           max(pc.scores[comp2])+5),
       axes=TRUE,
       xlab=paste("PC", comp1, " - ", pc.1, "%", sep=""),
       ylab=paste("PC", comp2, " - ", pc.2, "%", sep=""),
       main=paste("PCA - ", title, sep="")
       )

  # add labels:
  graphics::abline(h = 0, v = 0, col = "gray", lty = 2)
  # add legend:
  if(legend==TRUE){
    legend("bottomright", c(as.character(unique(pc.scores$group))),
                            col=colors,
                          pch = c(rep(19, length(unique(pc.scores$group)))),
           title = "SAMPLE GROUPS", inset = .02, cex=0.5)
  }
  if(label==TRUE){
    graphics::text(pc.scores[,comp1], pc.scores[,comp2], pc.scores$file,
                   cex=0.5, pos=3, col="black")
  }
  graphics::points(pc.scores[,comp1], pc.scores[,comp2], cex = 1,
                   col = colors[se.in$group], pch=19)
}

