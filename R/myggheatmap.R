#' Generate a ggplot2 heatmap with row and column dendrograms
#' modified from https://github.com/tleonardi/ggheatmap
#' @param dataMatrix A data.frame containing the input data.
#' @param orderCol Reorder the columns (default=T)
#' @param orderRow Reorder the rows (default=T)
#' @param dendroLineSize Size of the dendrogram lines (default=0.5)
#' @param revCol reverse the column order
#' @param revRow revers the row order
#' @param fontSize Font size (default=20)
#' @param vlines Add vertical lines (default = FALSE)
#' @param vlinecol Color of vertical lines if specified
#' @param vlinesize Size of vertical colors if specified
#' @param colorPalette Color palette (default='Spectral')
#' @param scaleName Name of the colorscale (default='value')
#' @param distMethod Distance method (default='euclidean', see ?dist)
#' @param clustMethod Clustering method (default='complete', see ?hclust)
#' @param revColors Invert color scale
#' @param title Title of the plot (default=NULL) 
#' @examples myggheatmap(mtcars)
#' @importFrom magrittr %>%
#' @import ggplot2
#' @export 
myggheatmap <- function(dataMatrix, orderCol = T, orderRow = T, points = F,
                        dendroLineSize = 0.5,
                        revCol = F,
                        revRow = F,
                        rowLabels = NA,
                        colLabels = NA,
                        fontSize = 20,
                        vlines = FALSE,
                        vlinecol = NULL,
                        vlinesize = NULL,
                        viridis = FALSE,
                        colorPalette = "default",
                        scaleName = "value",
                        distMethod = "euclidean",
                        clustMethod = "complete",
                        revColors=F,
                        title = NULL, ...) {
    
    data_m <- tibble::rownames_to_column(dataMatrix) %>% reshape2::melt()
    
    # Cluster rows
    if (orderRow) {
        dd.row <- as.dendrogram(hclust(dist(dataMatrix, method = distMethod),
                                       method = clustMethod))
        
        if (revRow) dd.row = rev(dd.row)
        
        row.ord <- order.dendrogram(dd.row)
        ordered_row_names <- row.names(dataMatrix[row.ord, , drop = FALSE])
        data_m$rowname <- factor(data_m$rowname, levels = ordered_row_names)
        
        
    }
    
    if (!is.na(rowLabels)) levels(data_m$rowname ) = rowLabels[levels(data_m$rowname )]
    
    # Cluster columns
    if (orderCol) {
        dd.col <- as.dendrogram(hclust(dist(t(dataMatrix), method = distMethod), 
                                       method = clustMethod))
        
        if (revCol) dd.col = rev(dd.col)
        
        col.ord <- order.dendrogram(dd.col)
        ordered_col_names <- colnames(dataMatrix[, col.ord])
        data_m$variable <- factor(data_m$variable, levels = ordered_col_names)
    }
    
    if (!is.na(colLabels)) levels(data_m$variable ) = colLabels[levels(data_m$variable )]
    
    if (colorPalette == "default") {
        if (viridis) {
            colorPalette = "viridis"
        } else {
            colorPalette = "Spectral"
        }
    }
    
    
    heat_plot <- ggplot2::ggplot(data_m, ggplot2::aes(x = variable, y = rowname))
    
    if (points == TRUE) {
        heat_plot = heat_plot +
            ggplot2::geom_point(aes(color = value, size = abs(value) ) ) +
            guides(size = "none")
        
        if (viridis) {
            heat_plot = heat_plot +
                viridis::scale_color_viridis(option = colorPalette,name = scaleName)
        } else {
            heat_plot = heat_plot +
                ggplot2::scale_color_distiller(palette = colorPalette, name = scaleName,
                                               direction=ifelse(revColors,1,-1), ...)
        }
        
    } else {
        heat_plot = heat_plot +
            ggplot2::geom_tile(aes(fill = value) )
        ## ggplot2::geom_vline(xintercept=1:nlevels(data_m$variable) + 0.5,
        ##                     color = "gray40" ,size=0.6)+
        if (viridis) {
            heat_plot = heat_plot +
                viridis::scale_fill_viridis(option = colorPalette, name = scaleName)
        } else {
            heat_plot = heat_plot +
                ggplot2::scale_fill_distiller(palette = colorPalette, name = scaleName,
                                              direction=ifelse(revColors,1,-1), ... )
        }
    }
    
    heat_plot = heat_plot +  
        ggplot2::theme_minimal() + 
        ggplot2::theme(axis.line = ggplot2::element_line(size = 0),
                       text = ggplot2::element_text(size = fontSize),
                       axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
                       panel.border = element_rect(color = "gray20", fill = NA) ) + 
        ggplot2::scale_y_discrete(position = "right") + 
        ggplot2::xlab("") + 
        ggplot2::ylab("")
    
    if (! is.null(title)) {
        heat_plot = heat_plot + ggtitle(title)
    }
    
    if (vlines) {
        if (is.null(vlinecol)) vlinecol = "gray40"
        if (is.null(vlinesize)) vlinesize = 0.8
        
        segment_data = data.frame(xstarts = 1:nlevels(data_m$variable) + 0.5,
                                  xends = 1:nlevels(data_m$variable) + 0.5,
                                  ystarts = 0.5, yends = nlevels(data_m$rowname)+ 0.5)
        
        heat_plot = heat_plot +
            ggplot2::geom_segment(aes(x = xstarts, y = ystarts, xend = xends, yend = yends),
                                  color = vlinecol, size = vlinesize,
                                  data = segment_data)
    }
    final_plot <- heat_plot
    
    if (orderRow) {
        dendro_data_row <- ggdendro::dendro_data(dd.row, type = "rectangle")
        dendro_row <- cowplot::axis_canvas(heat_plot, axis = "y", coord_flip = TRUE) + 
            ggplot2::geom_segment(data = ggdendro::segment(dendro_data_row),
                                  ggplot2::aes(y = -y, 
                                               x = x, xend = xend, yend = -yend), size = dendroLineSize) + 
            ggplot2::coord_flip()
        final_plot <- cowplot::insert_yaxis_grob(final_plot, dendro_row, grid::unit(0.2, 
                                                                                    "null"), position = "left")
    }
    
    if (orderCol) {
        dendro_data_col <- ggdendro::dendro_data(dd.col, type = "rectangle")
        dendro_col <- cowplot::axis_canvas(heat_plot, axis = "x") + 
            ggplot2::geom_segment(data = ggdendro::segment(dendro_data_col),
                                  ggplot2::aes(x = x, y = y, xend = xend, yend = yend), size = dendroLineSize)
        final_plot <- cowplot::insert_xaxis_grob(final_plot, dendro_col, grid::unit(0.2, 
                                                                                    "null"), position = "top")
    }
    
    cowplot::ggdraw(final_plot)
    
}
