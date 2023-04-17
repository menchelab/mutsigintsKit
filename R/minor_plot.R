#' expands the panel with a given factor
#' @param ggout Input ggplot file
#' @param title Title of the output file
#' @param expand.factor Expand the output with the given factor
#' @param expand.factor.y Expand factor for the y axis if they are different
#' 
#' @export

minor_plot = function(ggout, expand.factor, expand.factor.y = NULL, title = NULL) {
    
    bb = ggplot_build(ggout)
    ylims = bb$layout$panel_scales_y[[1]]$range$range
    xlims = bb$layout$panel_scales_x[[1]]$range$range
    
    xwidth  = abs(xlims[2] - xlims[1])
    ywidth = abs(ylims[2] - ylims[1])
    
    xlims.out = c(xlims[1] - expand.factor * xwidth,  xlims[2] + expand.factor * xwidth)
    
    if(! is.null(expand.factor.y) ) {
        ylims.out = c(ylims[1] - expand.factor.y * ywidth,  ylims[2] + expand.factor.y * ywidth)
    }else {
        ylims.out = c(ylims[1] - expand.factor * ywidth,  ylims[2] + expand.factor * ywidth)
    }

    if (! is.null(title)) {
        pout = ggout + ggtitle(title)
    } else {
        pout = ggout
    }
    
    if (! is.null(xlims) & ! is.null(ylims)) {
        pout = pout + 
            coord_cartesian(xlim = xlims.out, ylim = ylims.out)
    } else {
        pout = pout + coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))
    }
    return(pout)
}
