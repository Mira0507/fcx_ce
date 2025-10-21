#' Create a subchunk calling DT::datatable(dataframe) or plotly::ggplotly(plot)
#'
#' @param name string. subchunk name
#' @param input string. 'df' for data.frame or 'plot' for ggplot object
#' @param width integer. plot width
#' @param height integer. plot height
#' @param ggplotly logical. TRUE for interactive plotting using plotly::ggplot,
#'                          FALSE for non-interactive plotting
#' @return string to print a subchunk
subchunkify <- function(name, input='df', width=7, height=5, ggplotly=TRUE) {
    if (input == 'df') {
        t_deparsed <- paste0("DT::datatable(df)")
        more <- ""
    } else if (input == 'plot') {
        t_deparsed <- paste0("plotly::ggplotly(p)")
        if (!ggplotly) {
            t_deparsed <- paste0("print(p)")
        }
        more <- paste0(", fig.width=", width, ", fig.height=", height)
    } else {
        stop('Incorrect argument: input')
    }
    sub_chunk <- paste0("```{r sub_chunk_", name, ", results='asis', echo=FALSE", more, "}",
        "\n",
        t_deparsed,
        "\n\n```\n\n\n")

    cat(knitr::knit(text = sub_chunk, quiet=TRUE))
}


#' Print a Markdown link for a plot
#'
#' @param fn filepath
#' @return String ready to be printed
link_plot <- function(fn){
    cat(paste('\n\n- Download Plot: [', fn, '](', fn, ')\n\n'))
}

#' Print a Markdown link for a table
#'
#' @param fn filepath
#' @return String ready to be printed
link_table <- function(fn){
    cat(paste('\n\n- Download Table: [', fn, '](', fn, ')\n\n'))
}

#' Create sashimi plots
#'
#' @param cluster_to_plot String intron cluster (e.g. "clu_1_+")
#' @param main_title String plot title
#' @param exons_table data.frame exon table consisting of chr, start, end, strand, and gene_name columns
#' @param meta data.frame metadata frame consisting of sample and group in the first two columns
#' @param cluster_ids character vector consisting of cluster IDs
#' @param counts data.frame cluster-by-sample matrix for the number of extracted junctions
#' @param introns data.frame intron data frame consisting of clusterID, gene, ensemblID, chr, start, end,
#'                           verdict, deltapsi, and transcripts
#' @param snp_pos String coordinate of SNPs
#' @return grid of ggplot2 objects
#'
#' NOTE: This is a modified version of `make_cluster_plot` provided by the `LeafCutter` package
make_cluster_plot_ms <- function (cluster_to_plot,
                                  main_title = NA,
                                  exons_table = NULL,
                                  meta = NULL,
                                  cluster_ids = NULL,
                                  counts = NULL,
                                  introns = NULL,
                                  snp_pos = NA) {
    if (is.null(cluster_to_plot)) {
        print("no cluster selected!")
    }
    group_names <- levels(meta$group)
    y <- t(counts[cluster_ids == cluster_to_plot, ])
    x <- levels(meta$group)
    length_transform <- function(g) {
        log(g + 1)
    }

    introns_to_plot <- introns[introns$clusterID == cluster_to_plot, 
        ]
    summary_func <- colSums
    legend_title <- "Mean counts"
    alphabet <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", 
        "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", 
        "u", "v", "w", "x", "y", "z")
    junction_colour <- "red"
    cryptic_colour <- "pink"

    intron_meta <- get_intron_meta(colnames(y))
    intron_meta$verdict <- introns_to_plot$verdict[match(paste(intron_meta$start, 
        intron_meta$end), paste(introns_to_plot$start, introns_to_plot$end))]
    intron_meta$dPSI <- introns_to_plot$deltapsi[match(paste(intron_meta$start, 
        intron_meta$end), paste(introns_to_plot$start, introns_to_plot$end))]

    ranks <- alphabet[1:nrow(intron_meta)]
    absPSI <- intron_meta$dPSI[order(abs(intron_meta$dPSI), decreasing = TRUE)]
    intron_meta$rank <- ranks[match(intron_meta$dPSI, absPSI)]

    if (all(!grepl("chr", intron_meta$chr))) {
        intron_meta$chr <- paste0("chr", as.character(intron_meta$chr))
    }
    if (!is.null(exons_table)) {
        if (all(!grepl("chr", exons_table$chr))) {
            exons_table$chr <- paste0("chr", as.character(exons_table$chr))
        }
    }
    new_theme_empty <- theme_bw(base_size = 15)
    new_theme_empty$panel.background = element_rect(fill = "white", 
        colour = "white")
    new_theme_empty$line <- element_blank()
    new_theme_empty$rect <- element_blank()
    new_theme_empty$strip.text <- element_blank()
    new_theme_empty$axis.text <- element_blank()
    groups <- sort(unique(x))
    temp_max <- mclapply(groups, function(tis) {
            intron_meta$counts <- summary_func(y[tis == x, , drop = F])
        }) %>%
    unlist() %>%
    max()
    max_log <- 0.5 * ceiling(2 * log10(1 + temp_max))

    breaks <- if (max_log <= 2.5) {
        seq(0, max_log, by = 0.5)
    } else {
        seq(0, ceiling(max_log), by = 1)
    }

    limits <- c(0, max_log)
    intron_meta$id = as.factor(1:nrow(intron_meta))
    temp <- intron_meta[, c("id", "start", "end")]
    m <- reshape2::melt(temp, id.vars = "id")
    s <- unique(m$value)
    if (!is.na(snp_pos)) {
        s <- c(s, snp_pos)
    }
    s <- sort(s)
    d <- s[2:length(s)] - s[1:length(s) - 1]
    trans_d <- length_transform(d)
    coords <- c(0, cumsum(trans_d))
    names(coords) <- s
    snp_coord <- coords[as.character(snp_pos)]
    total_length <- sum(trans_d)
    my_xlim <- c(-0.05 * total_length, 1.05 * total_length)
    first_plot <- T
    min_height <- 0
    max_height <- 0
    curv <- 0.1
    min_exon_length <- 0.5
    maxratio <- 0
    minratio <- 1
    yFactor <- 0.65
    yConstant <- -0.25
    labelTextSize <- 3.5
    curveMax <- 10
    curveExponent <- 2
    yOffset <- 0
    centreLineWidth <- 3
    mainPalette <- c(junction_colour, cryptic_colour)
    summary_func <- function(a) apply(sweep(a, 1, rowSums(a), 
        "/"), 2, function(g) mean(g, na.rm = T))
    plots <- lapply(groups, function(tis) {
        intron_meta$counts <- summary_func(y[tis == x, , drop = F])
        maxratio <- max(c(max(intron_meta$counts/sum(intron_meta$counts)), maxratio))
        minratio <- min(c(min(intron_meta$counts/sum(intron_meta$counts)), minratio))
        list(maxratio=maxratio, minratio=minratio)
        })

    last_group <- groups[length(groups)]
    plots <- list()
    for (fancyVar in 1:length(groups)) {
        intron_meta$counts <- summary_func(y[groups[fancyVar] == x, , drop = F])
        intron_meta$prop <- intron_meta$counts
        group_sample_size <- sum(groups[fancyVar] == x)
        elist <- lapply(1:nrow(intron_meta), function(i) {
            if (i%%2 == 1) 
              return(NULL)
            start <- coords[as.character(intron_meta$start[i])]
            end <- coords[as.character(intron_meta$end[i])]
            l <- end - start
            edge <- data.frame(start, end)
            edge$startv <- intron_meta$start[i]
            edge$endv <- intron_meta$end[i]
            edge$start <- start
            edge$end <- end
            edge$log10counts = intron_meta$counts[i] + 1
            edge$label <- paste0(format(intron_meta$prop[i], 
              digits = 2, scientific = FALSE), "^", intron_meta$rank[i])
            edge$clu <- intron_meta$clu[i]
            edge$Group <- i
            edge$xtext <- start + l/2
            edge$ytext <- -((l^(yFactor)/2) + yConstant)
            edge$verdict <- ifelse(intron_meta$verdict[i] == 
              "annotated", yes = "annotated", no = "cryptic")
            edge
        })
        allEdges <- do.call(rbind, elist)
        eplist <- lapply(1:nrow(intron_meta), function(i) {
            if (i%%2 == 0) 
              return(NULL)
            start <- coords[as.character(intron_meta$start[i])]
            end <- coords[as.character(intron_meta$end[i])]
            l <- end - start
            edge <- data.frame(start, end)
            edge$startv <- intron_meta$start[i]
            edge$endv <- intron_meta$end[i]
            edge$start <- start
            edge$end <- end
            edge$log10counts = intron_meta$counts[i] + 1
            edge$label = paste0(format(intron_meta$prop[i], 
              digits = 2, scientific = FALSE), "^", intron_meta$rank[i])
            edge$clu <- intron_meta$clu[i]
            edge$Group <- i
            edge$xtext <- start + l/2
            edge$ytext <- l^(yFactor)/2 + yConstant
            edge$SIZE <- intron_meta$prop[i] + 1
            edge$verdict <- ifelse(intron_meta$verdict[i] == 
              "annotated", yes = "annotated", no = "cryptic")
            edge
        })
        allEdgesP <- do.call(rbind, eplist)
  
        if (all(is.na(main_title)) | !first_plot) {
            new_theme_empty$plot.title <- element_blank()
        }
        first_plot <- FALSE
        MAXcounts <- max(c(max(with(allEdgesP, 1)), max(with(allEdges, 1))))
        YLIMP <- max(allEdgesP$ytext) + 0.25 * max(allEdgesP$ytext)
        YLIMN <- min(allEdges$ytext) + 0.25 * min(allEdges$ytext)
        g <- ggplot() + geom_curve(data = allEdgesP, aes(x = start, 
            xend = xtext, y = yOffset, yend = ytext, group = Group, 
            colour = verdict, size = curveMax * (log10counts - 
                1)^curveExponent), angle = 90, curvature = -curv, 
            lineend = "round") + geom_curve(data = allEdgesP, 
            aes(x = xtext, xend = end, y = ytext, yend = yOffset, 
                group = Group, colour = verdict, size = curveMax * 
                  (log10counts - 1)^curveExponent), angle = 90, 
            curvature = -curv, lineend = "round") + geom_curve(data = allEdges, 
            aes(x = start, xend = xtext, y = -yOffset, yend = ytext, 
                group = Group, colour = verdict, size = curveMax * 
                  (log10counts - 1)^curveExponent), angle = 90, 
            curvature = curv, lineend = "round") + geom_curve(data = allEdges, 
            aes(x = xtext, xend = end, y = ytext, yend = -yOffset, 
                group = Group, colour = verdict, size = curveMax * 
                  (log10counts - 1)^curveExponent), angle = 90, 
            curvature = curv, lineend = "round") + new_theme_empty + 
            ylab(paste0(groups[fancyVar], " (n=", group_sample_size, 
                ")")) + xlab("") + xlim(my_xlim) + ggtitle(paste0(groups[fancyVar], 
            " (n=", group_sample_size, ")")) + geom_hline(yintercept = 0, 
            size = centreLineWidth, colour = "white") + geom_hline(yintercept = 0, 
            alpha = 0.9, size = 1) + geom_label(data = allEdgesP, 
            aes(x = xtext, y = 0.95 * ytext, label = label), 
            size = labelTextSize, label.size = NA, parse = TRUE, 
            fill = "white", colour = "black", label.r = unit(0.3, 
                "lines"), label.padding = unit(0.3, "lines")) + 
            geom_label(data = allEdges, aes(x = xtext, y = 0.95 * 
                ytext, label = label), size = labelTextSize, 
                label.size = NA, parse = TRUE, fill = "white", 
                colour = "black", label.r = unit(0.3, "lines"), 
                label.padding = unit(0.3, "lines")) + ylim(YLIMN, 
            YLIMP) + scale_size_continuous(limits = c(0, 10), 
            guide = "none")
        if (!is.na(snp_coord)) {
            df <- data.frame(x = snp_coord, xend = snp_coord, 
                y = 0, yend = max_height * 1.1)
            g <- g + geom_segment(data = df, aes(x = x, y = y, 
                xend = xend, yend = yend))
        }
        plots[[fancyVar]] <- g
    }
    df <- data.frame(x = coords,
                     xend = total_length * (s - min(s))/(max(s) - min(s)),
                     y = 0,
                     yend = min_height)
    if (!is.null(exons_table)) {
        exons_chr <- exons_table[exons_table$chr == intron_meta$chr[1], ]
        stopifnot(nrow(exons_chr) > 0)
        exons_here <- exons_chr[(min(s) <= exons_chr$start & 
            exons_chr$start <= max(s)) | (min(s) <= exons_chr$end & 
            exons_chr$end <= max(s)), ]
        if (nrow(exons_here) > 0) {
            exons_here <- unique(exons_here[(exons_here$end %in% 
                intron_meta$start | exons_here$start %in% intron_meta$end) & 
                (exons_here$end - exons_here$start <= 500 | exons_here$end == 
                  min(intron_meta$start) | exons_here$start == 
                  max(intron_meta$end)), ])
        }
        if (nrow(exons_here) > 0) {
            exons_here$gene_name <- factor(exons_here$gene_name)
            n_genes <- seq(1, length(levels(exons_here$gene_name)))
            gene_name_df <- data.frame(x = 0.2 * total_length, 
                y = YLIMN, label = rev(levels(exons_here$gene_name)))
            gene_name_df$N <- table(exons_here$gene_name)[gene_name_df$label]
            gene_name_df <- gene_name_df[order(gene_name_df$N, 
                decreasing = TRUE), ]
            cbbPalette <- c("#000000", "#E69F00", "#56B4E9", 
                "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
            cbbPalette <- cbbPalette[1:length(gene_name_df$label)]
            names(cbbPalette) <- gene_name_df$label
            mainPalette <- c(cbbPalette, junction_colour, cryptic_colour)
            names(mainPalette)[(length(mainPalette) - 1):length(mainPalette)] <- c("annotated", 
                "cryptic")
            invert_mapping <- function(pos) {
                if (pos %in% s) 
                  coords[as.character(pos)]
                else if (pos < min(s)) 
                  my_xlim[1]
                else if (pos > max(s)) 
                  my_xlim[2]
                else {
                  w = which(pos < s[2:length(s)] & pos > s[1:(length(s) - 
                    1)])
                  stopifnot(length(w) == 1)
                  coords[w] + (coords[w + 1] - coords[w]) * (pos - 
                    s[w])/(s[w + 1] - s[w])
                }
            }
            exon_df <- data.frame(x = sapply(exons_here$start, 
                invert_mapping), xend = sapply(exons_here$end, 
                invert_mapping), y = 0, yend = 0, label = exons_here$gene_name)
            exon_df[(exon_df$xend - exon_df$x) < min_exon_length, 
                ]$xend <- exon_df[(exon_df$xend - exon_df$x) < 
                min_exon_length, ]$x + min_exon_length
            exon_df <- exon_df[!duplicated(paste(exon_df$x, exon_df$xend)), ]
            principal_gene <- gene_name_df$label[1]
            gene_strand <- unique(exons_here[exons_here$gene_name == principal_gene, ]$strand)
            if (length(gene_strand) > 1 & !is.na(gene_strand)) {
                gene_strand <- NULL
            } else {
                exon_intervals <- intervals::Intervals(matrix(data = c(exon_df$x, 
                  exon_df$xend), ncol = 2))
                intron_intervals <- intervals::interval_complement(exon_intervals)
                intron_intervals <- intron_intervals[2:(nrow(intron_intervals) - 
                  1), ]
                strand_df <- as.data.frame(intron_intervals)
                strand_df <- strand_df[(strand_df$V2 - strand_df$V1) > 
                  0.025 * total_length, ]
                strand_df$midpoint <- strand_df$V1 + (strand_df$V2 - 
                  strand_df$V1)/2
                group <- c()
                for (i in 1:nrow(strand_df)) {
                  group <- c(group, rep(i, 2))
                }
                if (gene_strand == "+") {
                  strand_df <- data.frame(x = c(rbind(strand_df$V1, 
                    strand_df$midpoint)), group = group, y = 0)
                }
                if (gene_strand == "-") {
                  strand_df <- data.frame(x = c(rbind(strand_df$midpoint, 
                    strand_df$V2)), group = group, y = 0)
                }
                if (gene_strand == "+") {
                  strand_pos <- "last"
                }
                if (gene_strand == "-") {
                  strand_pos <- "first"
                }
                for (i in 1:length(plots)) {
                  plots[[i]] <- plots[[i]] + geom_line(data = strand_df, 
                    aes(x = x, y = y, group = group), colour = "black", 
                    size = 1, arrow = arrow(ends = strand_pos, 
                      type = "open", angle = 30, length = unit(0.1, 
                        units = "inches")))
                }
            }
            for (i in 1:length(plots)) {
                plots[[i]] <- plots[[i]] + geom_segment(data = exon_df, 
                  aes(x = x, y = y, xend = xend, yend = yend, 
                    colour = label), alpha = 1, size = 6) + geom_segment(data = exon_df, 
                  aes(x = x, xend = x + 0.01, y = y, yend = yend), 
                  colour = "white", size = 6, alpha = 1) + geom_segment(data = exon_df, 
                  aes(x = xend - 0.01, xend = xend, y = y, yend = yend), 
                  colour = "white", size = 6, alpha = 1)
            }
        }
    }
    if (all(!is.na(main_title))) {
        if (length(main_title) > 1) {
            plots[[1]] <- plots[[1]] + ggtitle(main_title[1], 
                subtitle = main_title[2]) + theme(plot.title = element_text(face = "bold.italic", 
                colour = "black", size = 15, hjust = 0.45), plot.subtitle = element_text(hjust = 0.45))
        }
        else {
            plots[[1]] = plots[[1]] + ggtitle(main_title)
        }
    }
    plots[[1]] <- plots[[1]] + scale_colour_manual("", values = mainPalette) + 
        guides(colour = FALSE)
    plots[[2]] <- plots[[2]] + scale_colour_manual("", values = mainPalette) + 
        theme(legend.position = "bottom", legend.justification = "right")
    if (is.na(snp_pos)) {
        gridExtra::grid.arrange(plots[[1]], plots[[2]], ncol = 1)
    }
}

