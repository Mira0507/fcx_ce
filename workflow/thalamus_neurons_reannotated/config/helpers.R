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

#' Summarize the number of samples (N) from input metadata
#'
#' @param meta_df data.frame metadata
#' @return data.frame summary data frame for the number of samples
summarize_N <- function(meta_df) {

    # Count the number of samples
    n_df <- meta_df %>%
        group_by_at(vars(subset_col, disease_col)) %>%
        summarize(N=n()) %>%
        spread(disease_col, N)

    # Replace missing values with zero
    n_df[is.na(n_df)] <- 0
    # Add a new column for all N
    n_df[['all_N']] <- rowSums(n_df[, colnames(n_df)[-1]])

    return(n_df)
}

#' Specify the plot format displaying counts
#'
#' @param input_df data.frame input data frame for plotting
#' @param dcol string column name specifying disease from the input data frame
#' @param scol string column name specifying sample ID
#' @param fcol string column name specifying sample IDs for sporadic FTD
#' @param ylabel string y-label
#' @return ggplot output beeswarm plot
count_plot <- function(input_df, dcol, scol, fcol, ylabel) {

    ggplot(input_df,
        aes_string(x=dcol, y='counts', fill=dcol)) +
        geom_violin(trim=TRUE) +
        theme_bw() +
        geom_quasirandom(data=input_df,
                         aes_string(x=dcol, y='counts', color=scol, shape=fcol),
                         width=0.2, size=3) +
        stat_summary(fun="median",
             geom="crossbar",
             width=0.1,
             color="black") +
        ylab(ylabel) +
        theme(axis.title.x=element_blank(),
        axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
        scale_fill_manual(values=color_map[[dcol]]) +
        scale_color_manual(values=color_map[[scol]])
}

#' Specify the plot format displaying ratios
#'
#' @param input_df data.frame input data frame for plotting
#' @param dcol string column name specifying disease from the input data frame
#' @param scol string column name specifying sample ID
#' @param fcol string column name specifying sample IDs for sporadic FTD
#' @param ylabel string y-label
#' @return ggplot output beeswarm plot
ratio_plot <- function(input_df, dcol, scol, fcol, ylabel) {

    ggplot(input_df, aes_string(x=dcol, y='ratio', fill=dcol)) +
        geom_violin(trim=TRUE) +
        geom_quasirandom(data=input_df,
                         aes_string(x=dcol, y='ratio', color=scol, shape=fcol),
                         width=0.2, size=3) +
           stat_summary(fun="median",
                        geom="crossbar",
                        width=0.1,
                        color="black") +
           theme_bw() +
           ylab(ylabel) +
           theme(axis.title.x=element_blank(),
                 axis.text.x=element_text(angle=90, vjust=0.5, hjust=1)) +
           scale_fill_manual(values=color_map[[dcol]]) +
           scale_color_manual(values=color_map[[scol]])
}

#' Extract the min and max values for ratio
#'
#' @param input_list list input nested list storing ratio data
#' @return vector min and max ratios
set_min_max <- function(input_list) {

    # Extract all ratios
    vec <- ratio_list %>%
        map(~.x %>%
            bind_rows() %>%
            pull(ratio)) %>%
            unlist()

    # Specify the max ratio value
    c(min=min(vec), max=max(vec))
}

#' Clean input label
#'
#' @param input_label string 
#' @return string cleaned group label
clean_label <- function(input_label) {
    ifelse(str_detect(input_label, "\\:"),
           str_replace(input_label, "\\:", "_"), input_label
    )
}

#' Calculate p-values using t-test or wilcox test
#'
#' @param dataframe data.frame input data frame
#' @param response_col string column name for the response variable
#' @param explanatory_col string column name for the explanatory variable
#' @param split_col string column name to subset rows 
#' @param stat string "ttest" or "wilcox"
#' @return p-value
calc_pval <- function(dataframe, response_col, explanatory_col, split_col, stat="ttest") {

    # Specify a design formula
    formula_str <- paste(response_col, "~", explanatory_col)
    sapply(unique(dataframe[[split_col]]), function(j) {
        df_i <- dataframe[dataframe[[split_col]] == j,]
        tryCatch(
            { 
                if (stat == "ttest") {
                    t.test(as.formula(formula_str), data=df_i)$p.value
                } else if (stat == "wilcox") {
                    wilcox.test(as.formula(formula_str), data=df_i, exact=FALSE)$p.value
                } else {
                    stop("Wrong statistics!")
                }
            },
            error = function(e) { 
                NA 
            })
        }
    )
}

#' Clean p-values and add FDR
#'
#' @param dataframe data.frame input data frame
#' @dcol string column name for disease
#' @control string disease condition
#' @padj_method string method for p-value adjustment
#' @return data frame where the output raw and adjusted p-values are added
clean_pval <- function(dataframe, dcol, control="Control", padj_method="BH") {
    # Replace p-value with 1 for rows from control group
    dataframe[['pvalue']] <- ifelse(dataframe[[dcol]] == control,
                                    1,
                                    dataframe[['pvalue']])
    # Add a column for FDR
    dataframe[['FDR']] <- p.adjust(dataframe[['pvalue']],
                                   method=padj_method)
    return(dataframe)
}

#' Subset and clean a count data frame
#'
#' @param dataframe data.frame input data frame
#' @param disease_label string for disease label
#' @return data frame with the number of junctions counted per sample within the contrast
subset_clean_counts <- function(dataframe, disease_label) {
    dataframe[dataframe[[disease_col]] %in% c('Control', disease_label),] %>%
        group_by_at(vars(sample_col,
                         "exon_detection",
                         "gene_spliced",
                         disease_col)) %>%
        summarize(counts=sum(counts))
}

#' Run a count division
#'
#' @param dataframe input data frame
#' @param denom string "nonce" (divide by non-CE junctions) or "all" (divide by all junctions)
#' @return output data frame with calculated ratios
calculate_ratios <- function(dataframe, denom="nonce") {

    # Prep input for a boxplot
    dataframe <- dataframe %>%
        pivot_wider(names_from=exon_detection, values_from=counts) %>%
        # Impute missing values with zero
        mutate_if(is.numeric, function(x) ifelse(is.na(x), 0, x))

    # Add a column for all counts
    dataframe[["all junctions"]] <- rowSums(dataframe[, e_detection])

    # Remove samples if they have zero counts for all exon_detection groups
    dataframe <- dataframe[rowSums(dataframe[, e_detection]) > 0,] %>%
        # Add a pseudocount of nonzero_min to all numeric entries
        mutate_if(is.numeric, function(x) x + pseudo_p * nonzero_min) %>%
        as.data.frame() %>%
        unique()

    # Add a new column calculating division
    denom_col <- ifelse(denom == "nonce",
                        "junction detected without CE",
                        "all junctions")
    dataframe[['ratio']] <- dataframe[["CE detected"]] / dataframe[[denom_col]]
    # Add a new column as a placeholder in case not needing to split the data frame
    dataframe[['metric']] <- 'Ratio'

    return(dataframe)
}


