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

