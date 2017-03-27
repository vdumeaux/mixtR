#' @title Break points and color vector for heatmaps.
#' @description
#' \code{heatmap.color.scheme} computes break points and color vector for heatmaps. Can also plot
#' the HCL space for easier customization of colors.
#'
#' @param col.palette either a palette name or index, or a full specification of an hcl palette.
#' See the 'supported.palettes' definition in the code for examples.
#' @param low.breaks  breakpoints for the low color
#' @param high.breaks breakpoints for the low color.
#' @param masked.color color to be used for the gap between the low and high colors.
#' @param plot if TRUE, will plot the HCL color space for the chosen color palette.
#' @param gp graphical parameters passed to grid.points and grid.lines when 'plot' is TRUE.
#'
#' @import grid
#'
#' @return a list with three fields:
#' \item{breaks}{vector with break points, size is one more that 'col'}
#' \item{col}{vector with colors, size one less than 'breaks'}
#' \item{masked.idx}{idx into 'col' and 'breaks' pointing to the color and
#' the left break point of the masked region (if such a region exists).}
#' @keywords internal
#'
heatmap.color.scheme <- function(low.breaks=seq(-3, 0, length.out=21),
                                 high.breaks=seq(0, 3, length.out=21),
                                 masked.color="black",
                                 col.palette=c("blue.red"),
                                 plot=FALSE, gp=gpar())
{
    stopifnot(length(low.breaks) > 0 || length(high.breaks) > 0)
    stopifnot(length(low.breaks) == 0 || is.numeric(low.breaks) && length(low.breaks) >= 2)
    stopifnot(length(high.breaks) == 0 || is.numeric(high.breaks) && length(high.breaks) >= 2)
    stopifnot(length(low.breaks) == 0 || all(diff(low.breaks) > 0))
    stopifnot(length(high.breaks) == 0 || all(diff(high.breaks) > 0))
    stopifnot(length(low.breaks) == 0 || length(high.breaks) == 0 || utils::tail(low.breaks, 1) <= utils::head(high.breaks, 1))

    supported.palettes <- list("blue.red"=list(h=c(260, 0), c=c(0, 80), l=c(85, 30), power=1.5),
                               "blue.red.2"=list(h=c(260, 0), c=c(0, 100), l=c(90, 50), power=1.5),
                               "blue.orange.1"=list(h=c(240, 60), c=c(0, 77), l=c(98, 72), power=1.5))

    if (is.character(col.palette) || is.numeric(col.palette))
        pal <- supported.palettes[[col.palette]]
    else
        pal <- col.palette

    if (!c("h", "c", "l", "power") %in% names(pal) ||
        !is.numeric(pal$h) || !is.numeric(pal$c) || !is.numeric(pal$l) || !is.numeric(pal$power) ||
        length(pal$h) != 2 || length(pal$c) != 2 || length(pal$l) != 2 || length(pal$power) != 1)
        stop("invalid 'col.palette' specification")

    color.points <- function(b) {
        if (length(b) == 2)
        {
            ret <- 1
        }
        else if (length(b) == 3)
        {
            ret <- c(0, 1)
        }
        else
        {
            mid.b <- b[-c(1, length(b))]
            ret <- c(b[1], mid.b[-length(mid.b)] + diff(mid.b)/2, b[length(b)])
            ret <- ret - ret[1]
            ret <- ret / utils::tail(ret, 1)
        }
        ret
    }

    if (length(low.breaks) > 0)
    {
        low.color.points <- color.points(low.breaks)^pal$power
        low.L <- pal$l[1] + diff(pal$l) * low.color.points
        low.C <- pal$c[1] + diff(pal$c) * low.color.points
        low.colors <- rev(colorspace::hex(colorspace::polarLUV(L=low.L, C=low.C, H=pal$h[1]), fixup=TRUE))
    }

    if (length(high.breaks) > 0)
    {
        high.color.points <- color.points(high.breaks)^pal$power
        high.L <- pal$l[1] + diff(pal$l) * high.color.points
        high.C <- pal$c[1] + diff(pal$c) * high.color.points
        high.colors <- colorspace::hex(colorspace::polarLUV(L=high.L, C=high.C, H=pal$h[2]), fixup=TRUE)
    }


    scheme <- list()
    if (length(low.breaks) == 0)
    {
        scheme$breaks <- high.breaks
        scheme$col <- c(high.colors)
    }
    else if (length(high.breaks) == 0)
    {
        scheme$breaks <- low.breaks
        scheme$col <- c(low.colors)

    }
    else if (utils::tail(low.breaks, 1) == utils::head(high.breaks, 1))
    {
        scheme$breaks <- c(low.breaks, high.breaks[-1])
        scheme$col <- c(low.colors, high.colors)
    }
    else
    {
        scheme$breaks <- c(low.breaks, high.breaks)
        scheme$col <- c(low.colors, masked.color, high.colors)
        scheme$masked.idx <- length(low.colors) + 1
    }

    if (isTRUE(plot))
    {
        plot.h <- function(h, pixels)
        {
            c <- seq(0,100,length=pixels)[-1];
            l <- seq(0,100,length=pixels)[-1];
            xy <- expand.grid(c=c, l=l);
            rect.col <- grDevices::hcl(h, xy$c, xy$l, fixup=FALSE)
            grid.rect(xy$c/100, xy$l/100, width=1/length(c), height=1/length(l),
                      gp=gpar(fill=rect.col, col=rect.col, linejoin="mitre", lwd=0.01), just=c("right", "top"),
                      default.units="native")
            grid.text(h, 0.8, 0.1, default.units="native")
        }

        pushViewport(viewport(layout=grid.layout(3, 2, heights=c(4, 1, 4))))

        pushViewport(viewport(layout.pos.row=1, layout.pos.col=1, xscale=c(1, 0), clip="off"))
        plot.h(pal$h[1], 100)
        grid.lines(pal$c/100, pal$l/100, default.units="native", gp=gp)
        grid.points(low.C / 100, low.L / 100, default.units="native", gp=gp)
        popViewport()

        pushViewport(viewport(layout.pos.row=1, layout.pos.col=2))
        plot.h(pal$h[2], 100)
        grid.lines(pal$c/100, pal$l/100, default.units="native", gp=gp)
        grid.points(high.C / 100, high.L / 100, default.units="native", gp=gp)
        popViewport()

        pushViewport(viewport(layout.pos.row=2, layout.pos.col=1:2, xscale=range(scheme$breaks)))
        grid.rect()
        grid.rect(scheme$breaks[-1], 0, width=diff(scheme$breaks), height=1,
                  gp=gpar(fill=scheme$col, col=scheme$col, linejoin="mitre", lwd=0.01),
                  just=c("right", "bottom"),
                  default.units="native")
        popViewport()

        pushViewport(viewport(layout.pos.row=3, layout.pos.col=1:2, xscale=range(scheme$breaks)))
        pushViewport(viewport(layout=grid.layout(4, 9)))

        rc <- expand.grid(c=1:9, r=1:4)
        h <- seq(0, 359, 10)
        for (i in 1:length(h))
        {
            pushViewport(viewport(layout.pos.row=rc[i, "r"], layout.pos.col=rc[i, "c"]))
            plot.h(h[i], 80)
            popViewport()
        }

        popViewport()
        popViewport()

        popViewport()
    }

    return(scheme)
}


#' @title Function for checking color schemes
#' @name heatmap.is.valid.color.scheme
#'
#' @description
#' \code{heatmap.is.valid.color.scheme} internal function for checking color schemes.
#' @keywords internal
heatmap.is.valid.color.scheme <- function(color.scheme)
{
    is.list(color.scheme) &&
    (c("breaks", "col") %in% names(color.scheme)) &&
    is.numeric(color.scheme$breaks) &&
    (length(color.scheme$breaks) == length(color.scheme$col) + 1) &&
    (is.null(color.scheme$masked.idx) || color.scheme$masked.idx %in% 1:length(color.scheme$col))
}


#' @title Plot heatmap
#' @description
#' \code{heatmap.map} plots a heatmap of 'x' in the current viewport. Returns the color
#' scheme used.
#' @param x a numeric matrix with the heatmap values
#' @param color.scheme  see the function heatmap.color.scheme() for details
#' @param na.color color used for NA values in 'x'
#' @param col.grid color used for grid lines in heatmap
 #'
#' @import grid stats
#'
#' @return The color scheme invisibly.
#' @keywords internal
heatmap.map <- function(x, color.scheme=heatmap.color.scheme(), na.color = "grey", col.grid=NA, use.raster=FALSE)
{
    if (length(dim(x)) != 2 || !is.numeric(x))
        stop("'x' must be a numeric matrix")

    bc <- color.scheme
    stopifnot(heatmap.is.valid.color.scheme(bc))

    ## make sure no values are outside the range ['min.val', 'max.val']
    min.val <- utils::head(bc$breaks, 1)
    max.val <- utils::tail(bc$breaks, 1)
    x[x < min.val] <- min.val
    x[x > max.val] <- max.val

    ## convert numbers to colors according to our breaks
    x.colors <- findInterval(x, bc$breaks, all.inside=TRUE)

    ## add na.color at the end of the color scheme
    bc$col <- c(bc$col, na.color)
    x.colors[is.na(x.colors)] <- length(bc$col)

    ## this viewport is set up so that origo is on the top left and the
    ## ranges of the "native" units are number of columns and rows of x.
    pushViewport(viewport(xscale=c(0, ncol(x)), yscale=c(nrow(x),0)))

    ## now draw the heatmap
    if (isTRUE(use.raster))
    {
        colmatrix <- bc$col[x.colors]
        dim(colmatrix) <- dim(x)
        grid.raster(colmatrix, interpolate=FALSE, width=1, height=1)
    }
    else
    {
        grid.rect(col(x), row(x), gp=gpar(col=bc$col[x.colors], fill=bc$col[x.colors], linejoin="mitre", lwd=0.01, col=col.grid),
                  width=1, height=1, just=c("right", "top"), default.units="native")
    }

    popViewport()

    return(invisible(color.scheme))
}


#' Plot heatmap key
#'
#' @description
#' \code{heatmap.key} plots a heatmap key in the current viewport.
#'
#' @param x numeric matrix that was used to plot a heatmap
#' @param color.scheme the color scheme used to plot the heatmap
#' @param  key.min should be less than the least break point of the color scheme
#' @param key.max should be greater than the greates break point of the color scheme
#' @param key.line.color  the color used to draw the histogram
#'
#' @import grid stats
#'
#' @return invisible(NULL)
#' @keywords internal
heatmap.key <- function(x, color.scheme=heatmap.color.scheme(),
                        key.min = utils::head(color.scheme$breaks, 1) - 1,
                        key.max = utils::tail(color.scheme$breaks, 1) + 1,
                        key.line.color = "cyan")
{
    if (!is.numeric(x))
        stop("'x' must be a numeric matrix")
    stopifnot(heatmap.is.valid.color.scheme(color.scheme))

    bc <- color.scheme

    min.val <- utils::head(bc$breaks, 1)
    max.val <- utils::tail(bc$breaks, 1)
    key.min <- min(key.min, min.val)
    key.max <- max(key.max, max.val)

    ## add breakpoints and colors for the masked region, from min.val to
    ## key.min, and from max.val to key.max. This is so that we can get
    ## a good histogram showing in the key.
    break.width <- min(diff(bc$breaks))

    if (!is.null(bc$masked.idx)) # extend masked region
    {
        masked.color <- bc$col[bc$masked.idx]
        a <- bc$breaks[bc$masked.idx]
        b <- bc$breaks[bc$masked.idx + 1]
        n <- round((b - a) / break.width)
        n <- max(1, n)
        new.breaks <- seq(a, b, length.out = n + 1)[-c(1, n+1)]
        bc$breaks <- c(bc$breaks[1:bc$masked.idx],
                       new.breaks,
                       bc$breaks[(bc$masked.idx + 1):length(bc$breaks)])
        bc$col <- c(bc$col[1:bc$masked.idx],
                    rep(masked.color, n - 1),
                    bc$col[(bc$masked.idx + 1):length(bc$col)])
    }
    if (key.min < min.val) # extend region between min.val and key.min
    {
        n <- round((min.val - key.min) / break.width)
        n <- max(1, n)
        new.breaks <- seq(key.min, min.val, length.out = n + 1)[-(n+1)]
        bc$breaks <- c(new.breaks, bc$breaks)
        bc$col <- c(rep(bc$col[1], n), bc$col)
    }
    if (key.max > max.val) # extend region between max.val and key.max
    {
        n <- round((key.max - max.val) / break.width)
        n <- max(1, n)
        new.breaks <- seq(max.val, key.max, length.out = n + 1)[-1]
        bc$breaks <- c(bc$breaks, new.breaks)
        bc$col <- c(bc$col, rep(bc$col[length(bc$col)], n))
    }

    ## draw the key

    ## we will leave 30% of the area below the key histogram for drawing
    ## an axis, which needs about 2 lines of text.
    axis.text <- do.call(paste, as.list(grid.pretty(c(key.min, key.max))))
    height <- convertHeight(unit(2, "lines"), "points", valueOnly = TRUE)/convertHeight(unit(1, "npc"), "points", valueOnly = TRUE)
    axis.cex <- min((as.double(convertHeight(unit(height, "native"), "lines")) / 2),
                    (1/as.double(convertWidth(unit(1, "strwidth", axis.text), "npc"))))*0.8
    height <- height * axis.cex

    pushViewport(viewport(0.5, 0.9, width=0.8, height=1-height, xscale=range(bc$breaks), yscale=c(1, 0),
                          just=c("centre", "top")))
    grid.rect(bc$breaks[-1], 0,
              width=diff(bc$breaks), height=1, gp=gpar(col=NA, fill=bc$col, lty=0),
              default.units="native", just=c("right", "bottom"))


    if (!all(is.na(x)))
    {
        h <- graphics::hist(x, plot = FALSE,
                  breaks = c(min(x, na.rm=TRUE), bc$breaks, max(x, na.rm = TRUE)))
        pushViewport(viewport(xscale=c(key.min, key.max), yscale=c(0, max(h$counts)*1.05)))
        breaks <- h$breaks[-c(1, length(h$breaks))]
        counts <- h$counts[-c(1, length(h$counts))]
        grid.lines(rep(breaks, each=2), c(0, rep(counts, each=2), 0),
                   gp=gpar(col=key.line.color), default.units="native")
        grid.xaxis(gp=gpar(cex=axis.cex, lex=axis.cex))
        popViewport()
    }
    else
    {
        pushViewport(viewport(xscale=c(key.min, key.max)))
        grid.xaxis(gp=gpar(cex=axis.cex, lex=axis.cex))
        popViewport()
    }
    popViewport()

    return(invisible(NULL))
}


#' Plot clinical variables
#'
#' @description
#' \code{heatmap.clinical} plots clinical variables as boxes in the current frame.  Assumes 'x'
#' is a matrix with color information in each cell. Colors are usually
#' given using character strings such as "white", "red", etc. See
#' 'colors()' for a list of color names. NAs are plotted with color as
#' given by 'na.color'.
#' @param x matrix with colors in each cell, or NAs.
#' Will also accept data frames (with rows representing patients)
#' @param na.color color to be used for missing values, i.e., NAs
#'
#' @import grid
#'
#' @return NULL
#' @keywords internal
heatmap.clinical <- function(x, na.color = "grey")
{
    if (length(dim(x)) != 2)
        stop("'x' must be a color matrix")

    if (is.data.frame(x))
        x <- t(as.matrix(x))

    x[is.na(x)] <- na.color

    pushViewport(viewport(xscale=c(0, ncol(x)), yscale=c(nrow(x), 0)))

    grid.rect(col(x), row(x), width=1, height=1, gp=gpar(col=x, fill=x, linejoin="mitre", lwd=0.01),
              just=c("right", "top"), default.units="native")
    grid.rect(ncol(x), 1:nrow(x), width=ncol(x), height=1, gp=gpar(col="black", fill=NA),
              just=c("right", "top"), default.units="native")

    popViewport()

    return(invisible(NULL))
}


#' Flexible clustering function
#'
#' @description
#' \code{heatmap.ddr} is a very flexible function for clustering. It allows multiple ways of
#' specifying how the dendrogram is to be computed and it knows about
#' the cor() and hclust() functions.
#' @param exprs a matrix-like object whose columns are to be clustered.
#' May be omitted if a dendrogram can be obtained using only the other arguments.
#' @param dist.elem specifies how to compute the distance matrix. May be
#' omitted if a dendrogram can be obtained from 'clust.elem' alone.
#' Possible specifications are: An object that can be coerced to a distance matrix
#' (using as.dist()), a function that takes 'exprs' and
#' produces a distance matrix, or a character string specifying
#' how a distance matrix is to be computed from 'exprs'.
#' The acceptable character strings are "cor" and any method specification that can be
#' passed on to the dist() function. If 'dist.elem' is "cor", then 1-correlation
#' (between the columns of 'exprs') is used as distance.
#' @param clust.elem specifies how to compute the dendrogram. Possible
#' specifications are:
#' An object that can be coerced to "dendrogram" (using as.dendrogram()), a function
#' that takes a distance matrix and produces an object that can be coerced to a dendrogram,
#' a character string specifying how to compute, from a distance
#' matrix, an object that can be coerced to a dendrogram, or a character
#' string specifying how to compute a dendrogram from a distance matrix.
#' passed to the dist() function and used only when 'dist.elem' is "minkowski"
#' (see help on the dist() function for details).
#'
#' @import stats
#'
#' @keywords internal
heatmap.ddr <- function(exprs=NULL, dist.elem="cor", clust.elem="average", p=2)
{
    if (class(clust.elem) != "dendrogram")
    {
        ret <- try(as.dendrogram(clust.elem), silent=TRUE)
        if (class(ret) == "try-error")
        {
            if (mode(dist.elem) == "numeric" && class(dist.elem) == "matrix")
            {
                dist.mat <- as.dist(dist.elem)
            }
            else if (is.null(exprs))
            {
                stop("heatmap.ddr: 'exprs' missing")
            }
            else if (mode(dist.elem) == "function")
            {
                dist.mat <- as.dist(dist.elem(exprs))
            }
            else if (mode(dist.elem) == "character" && length(dist.elem) == 1 && dist.elem == "cor")
            {
                cors <- cor(exprs, use="pairwise.complete.obs")
                cors[is.na(cors)] <- 0
                dist.mat <- as.dist(1 - cors)
            }
            else if (mode(dist.elem) == "character" && length(dist.elem) == 1)
            {
                dist.mat <- dist(t(exprs), method=dist.elem, p=p)
            }
            else
            {
                stop("heatmap.ddr: bad 'dist.elem'")
            }

            if (mode(clust.elem) == "function")
            {
                ret <- as.dendrogram(clust.elem(dist.mat))
            }
            else if (mode(clust.elem) == "character" && length(clust.elem) == 1)
            {
                ret <- as.dendrogram(hclust(dist.mat, method=clust.elem))
            }
            else
            {
                stop("heatmap.ddr: bad 'dist.elem'")
            }
        }
    }
    else
    {
        ret <- clust.elem
    }

    return(ret)
}


#' Return a cex value for text plotted in a viewport
#'
#' @description
#' \code{heatmap.labels.cex}computes a cex value that can be used when plotting text in the
#' current viewport. The assumption is that 'labels' is a character
#' vector containing text strings that are meant to be drawn in the
#' current viewport (either horizontally or vertically). This function
#' will compute and return a cex value that can be passed to
#' grid.text() such that the strings in 'labels' will fit in the
#' current viewport when drawn on separate lines.
#'
#' @param labels character vector
#' @param type either "row.labels" or "col.labels" indicating whether the
#'         labels are intended to be drawn horizontally or vertically.
#'
#' @return a cex value to be used when drawing the strings in 'labels' in
#' the current viewport
#' @keywords internal
#' @import grid
#'
heatmap.labels.cex <- function(labels, type=c("row.labels", "col.labels"))
{
    type <- match.arg(type)

    rw <- convertWidth(max(do.call(unit.c, lapply(labels, function(label) {unit(1, "strwidth", label)}))), "npc")
    rw <- convertWidth((rw + unit(2, "strwidth", "X")), "npc")
    rh <- convertHeight(unit(1, "strheight", "X"), "npc") * 1.8 * length(labels)
    if (type == "row.labels")
    {
        cex <- min(1/as.double(convertWidth(rh, "npc")), 1/as.double(rw))
    }
    else
    {
        rw <- convertUnit(rw, "npc", "x", "dimension", "y", "dimension")
        rh <- convertUnit(rh, "npc", "y", "dimension", "x", "dimension")
        cex <- min(1/as.double(convertWidth(rh, "npc")), 1/as.double(rw))
    }
    return(cex)
}


#' Plot labels in a viewport
#'
#' @description
#' \code{heatmap.labels} plots each string in 'labels' on a separate line in the current
#' viewport. The size of the text and spacing between strings are
#' adjusted so that the text will fill the viewport as much as possible.
#'
#' @param labels character vector with strings to be drawn in the current
#'         viewport
#' @param type either "row.labels" or "col.lables", determining if the text
#'         is to be drawn horizontally or vertically.
#' @param just the justification of the text.
#' @param rot rotation of the text in angles.
#'
#' @return invisible(NULL)
#' @keywords internal
#' @import grid
heatmap.labels <- function(labels, type=c("row.labels", "col.labels"), just=c("left", "right", "centre"), rot=0)
{
    type <- match.arg(type)
    just <- match.arg(just)

    if (type == "row.labels")
    {
        pushViewport(viewport(xscale=c(0, 1), yscale=c(1, 0)))
        if (just == "left")
            x <- unit(1, "strwidth", "X")
        else if (just == "right")
            x <- unit(1, "npc") - unit(1, "strwidth", "X")
        else
            x <- 0.5
        y <- ((1:length(labels)) - 0.5) /length(labels)
        grid.text(labels, x, y, default.units="native", just=just, rot = rot,
                  gp=gpar(cex=heatmap.labels.cex(labels, type=type)))
        popViewport()
    }
    else
    {
        pushViewport(viewport(xscale=c(0, 1), yscale=c(0, 1)))
        x <- ((1:length(labels)) - 0.5) /length(labels)
        if (just == "left")
            y <- unit(1, "strwidth", "X")
        else if (just == "right")
            y <- unit(1, "npc") - unit(1, "strwidth", "X")
        else
            y <- 0.5
        grid.text(labels, x, y, default.units="native", just=just, rot=90 + rot,
                  gp=gpar(cex=heatmap.labels.cex(labels, type=type)))
        popViewport()
    }
    invisible(NULL)
}


#' @title Draw full heatmap figures
#' @name heatmap.simple
#'
#' @description Provides an easy way to produce full heatmap figures with the most
#' common features such as dendrograms from clustering, row and column
#' labels and a matrix of clinical variables. Many of the arguments
#' will take on sensible defaults if omitted (as described below).
#'
#' This function makes extensive use of the grid package and mixes it
#' with standard R graphics plots using gridBase. The figure is divided
#' into regions (with grid.layout which is similar to base R's layout
#' function) as directed by the layout.mat character matrix. The
#' strings in layout.mat will determine where each component of the
#' figure will be drawn. Other arguments are as described for the
#' heatmap functions above.
#'
#' @param x matrix of values for the main heatmap.
#' @param layout.mat a character matrix describing the layout of
#' the figure and its parts. The figure will be divided into rows and columns
#' corresponding to the rows and columns of layout.mat.
#' Each entry should be a string describing what should be plotted in that region.
#' The list of recognized entries are:
#' * heatmap
#' * key
#' * title
#' * row.labels.ljust
#' * row.labels.rjust
#' * col.labels.ljust
#' * col.labels.rjust
#' * row.ddr.left
#' * row.ddr.right
#' * col.ddr.up
#' * col.ddr.down
#' * clinical
#' * clinical.labels.ljust
#' * clinical.labels.rjust
#' If layout.mat is omitted, the function will try to provide a sensible
#' one depending on what other options have been supplied.
#' @param widths relative widths for the layout
#' @param heights relative heights for the layout.
#' If 'layout.mat' is omitted, care should be taken when supplying these
#' since they have to match the size of the default layout.mat, which can vary
#' depending on other arguments.
#' @param color.scheme a list defining the color scheme to be used for the
#' heatmap. See the heatmap.color.scheme() function for details
#' on the requirements for this variable
#' @param scale for now only supports "row", "none", or a user-supplied
#' function that will be applied to each row.
#' @param row.labels If omitted, will default to rownames of x.
#' @param col.labels If omitted, will default to colnames of x.
#' @param row.clust determines whether row should be clustered. If set to TRUE,
#' method is determined by row.clust.dist.fun and row.clust.fun
#' @param row.clust.dist.fun specifies how to compute the distance matrix for row clustering.
#' The acceptable character strings are "cor" and any method specification that can be
#' passed on to the dist() function. If 'dist.elem' is "cor", then 1-correlation
#' (between the columns of 'exprs') is used as distance.
#' @param row.clust.fun specifies how to compute the row dendrogram. Possible
#' specifications are: an object that can be coerced to "dendrogram" (using as.dendrogram()),
#' a function that takes a distance matrix and produces an object that can be coerced to a dendrogram,
#' a character string specifying how to compute, from a distance
#' matrix, an object that can be coerced to a dendrogram, or a character
#' string specifying how to compute a dendrogram from a distance matrix. (see dist() for more details)
#' @param col.clust determines whether columns should be clustered. If set to TRUE,
#' clustering method is determined by col.clust.dist.fun and col.clust.fun
#' @param col.clust.dist.fun specifies how to compute the distance matrix for column clustering.
#' The acceptable character strings are "cor" and any method specification that can be
#' passed on to the dist() function. If 'dist.elem' is "cor", then 1-correlation
#' (between the columns of 'exprs') is used as distance.
#' @param col.clust.fun specifies how to compute the row dendrogram. Possible
#' specifications are: an object that can be coerced to "dendrogram" (using as.dendrogram()),
#' a function that takes a distance matrix and produces an object that can be coerced to a dendrogram,
#' a character string specifying how to compute, from a distance
#' matrix, an object that can be coerced to a dendrogram, or a character
#' string specifying how to compute a dendrogram from a distance matrix. (see dist() for more details)
#' @param clinical a matrix of colors or a data frame to be passed to heatmap.clinical.
#' @param clinical.labels If omitted, will default to rownames of clinical'.
#' @param title text string to be printed as a title
#' @param lwd line width relative to the default (default=0.01)
#'
#' @keywords internal
#' @import grid stats gridBase
#'
heatmap.simple <- function(exprs,
                           layout.mat=NULL,
                           widths=NULL,
                           heights=NULL,
                           color.scheme=heatmap.color.scheme(),
                           key.min=utils::head(color.scheme$breaks, 1) - 1,
                           key.max=utils::tail(color.scheme$breaks, 1) + 1,
                           na.color="grey",
                           scale="row",
                           key.line.color="cyan",
                           row.labels=NULL,
                           col.labels=NULL,
                           row.clust=TRUE,
                           row.clust.dist.fun="cor",
                           row.clust.fun="average",
                           col.clust=TRUE,
                           col.clust.dist.fun="cor",
                           col.clust.fun="average",
                           col.grid=NA,
                           clinical=NULL,
                           clinical.labels=NULL,
                           clinical.na.color=na.color,
                           title=NULL,
                           use.raster=FALSE)
{

    ## save all standard R graphical parameters
    old.par <- graphics::par(no.readonly=TRUE)

    ## if clinical is a vector or data.frame, then convert it to a matrix
    if (!is.null(clinical))
    {
        if (is.vector(clinical))
            clinical <- matrix(clinical, nrow=1)
        else if (is.data.frame(clinical))
            clinical <- t(as.matrix(clinical))
    }

    ## get default values for arguments
    if (is.null(row.labels))
    {
        row.labels <- rownames(exprs)
        if (is.null(row.labels))
            row.labels <- rep("", nrow(exprs))
    }
    if (is.null(col.labels))
    {
        col.labels <- colnames(exprs)
        if (is.null(col.labels))
            col.labels <- rep("", ncol(exprs))
    }
    if (is.null(clinical.labels) && !is.null(clinical))
        clinical.labels <- rownames(clinical)

    ## if layout.mat is omitted, try to come up with a good default layout
    if (is.null(layout.mat))
    {
        layout.mat <- matrix(c("key", "col.ddr.up", "",
                               "row.ddr.left", "heatmap", "row.labels.ljust",
                               "", "col.labels.rjust", ""), ncol=3, byrow=TRUE)
        if (!is.null(clinical))
        {
            layout.mat <- rbind(layout.mat, c("clinical.labels.rjust", "clinical", ""))
            layout.mat <- rbind(layout.mat, c("", "", ""))
        }
        if (!is.null(title))
        {
            layout.mat <- rbind(c("", "title", ""), layout.mat)
        }
        if (is.null(widths))
        {
            widths <- c(0.3, 1, 0.3)
        }
        if (is.null(heights))
        {
            heights <- rep(0.3, nrow(layout.mat))
            heights[which(layout.mat == "heatmap", arr.ind=TRUE)[1]] <- 1
            if (!is.null(clinical))
            {
                heights[length(heights) - 1] <- 0.03*nrow(clinical)
                heights[length(heights)] <- 0.03
            }
            if (!is.null(title))
            {
                heights[1] <- 0.15
            }
        }
        if (!isTRUE(row.clust))
        {
            layout.mat[grep("row.ddr", layout.mat)] <- ""
        }
        if (!isTRUE(col.clust))
        {
            layout.mat[grep("col.ddr", layout.mat)] <- ""
        }

        stopifnot(length(widths) == ncol(layout.mat))
        stopifnot(length(heights) == nrow(layout.mat))
    }

    ## make sure layout is character matrix
    mode(layout.mat) <- "character"
    layout.mat <- as.matrix(layout.mat)

    ## scale the rows (if indicated by caller). When using correlation
    ## as distance, we definitely need to scale the rows *before*
    ## clustering since the mean of the columns don't really mean
    ## anything (expressions of different genes are not comparable in
    ## microarray experiments). Are there situations where we might want
    ## to cluster before scaling?
    if (is.function(scale))
    {
        exprs <- t(apply(exprs, 1, scale))
    }
    else
    {
        stopifnot(length(scale) == 1 && scale %in% c("row", "none"))
        ## scale the matrix if given
        if (scale == "row")
            exprs <- t(scale(t(exprs)))
    }

    ## get dendrograms and reorder the data matrix
    row.ddr <- NULL
    if (isTRUE(row.clust))
    {
        row.ddr <- rev(heatmap.ddr(t(exprs), row.clust.dist.fun, row.clust.fun))
        exprs <- exprs[order.dendrogram(rev(row.ddr)), ]
        row.labels <- row.labels[order.dendrogram(rev(row.ddr))]
    }

    col.ddr <- NULL
    if (isTRUE(col.clust))
    {
        col.ddr <- heatmap.ddr(exprs, col.clust.dist.fun, col.clust.fun)
        exprs <- exprs[ , order.dendrogram(col.ddr)]
        col.labels <- col.labels[order.dendrogram(col.ddr)]
        if (!is.null(clinical))
            clinical <- clinical[, order.dendrogram(col.ddr), drop=FALSE]
    }

    ## create the top viewport
    the.layout <- grid.layout(nrow(layout.mat), ncol(layout.mat), widths=widths, heights=heights)
    top.vp <- viewport(layout=the.layout, name="heatmap.top.vp")
    pushViewport(top.vp)

    ## get all unique strings in 'layout.mat'
    elements <- unique(as.vector(layout.mat))
    elements <- elements[elements != "" & !is.na(elements)]

    ## for the known strings in layout.mat, seek the named viewport and
    ## draw accordingly
    elem <- "heatmap"
    if (elem %in% elements)
    {
        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))
        heatmap.map(exprs, color.scheme=color.scheme, na.color=na.color, use.raster=use.raster, col.grid=col.grid)
        upViewport()
    }
    elem <- "key"
    if (elem %in% elements)
    {
        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))
        heatmap.key(exprs, color.scheme=color.scheme, key.min=key.min, key.max=key.max, key.line.color=key.line.color)
        upViewport()
    }
    elem <- "row.labels.ljust"
    if (elem %in% elements)
    {
        if (length(row.labels) != nrow(exprs))
            warning("the number of row labels is not equal to the number of rows of exprs")

        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))

        row.labels[is.na(row.labels)] <- ""
        heatmap.labels(row.labels, type="row.labels", just="left")
        upViewport()
    }
    elem <- "row.labels.rjust"
    if (elem %in% elements)
    {
        if (length(row.labels) != nrow(exprs))
            warning("the number of row labels is not equal to the number of rows of exprs")

        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))

        row.labels[is.na(row.labels)] <- ""
        heatmap.labels(row.labels, type="row.labels", just="right")
        upViewport()
    }
    elem <- "col.labels.ljust"
    if (elem %in% elements)
    {
        if (length(col.labels) != ncol(exprs))
            warning("the number of column labels is not equal to the number of columns of exprs")

        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))

        col.labels[is.na(col.labels)] <- ""
        heatmap.labels(col.labels, type="col.labels", just="left")
        upViewport()
    }
    elem <- "col.labels.rjust"
    if (elem %in% elements)
    {
        if (length(col.labels) != ncol(exprs))
            warning("the number of column labels is not equal to the number of columns of exprs")

        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))

        col.labels[is.na(col.labels)] <- ""
        heatmap.labels(col.labels, type="col.labels", just="right")
        upViewport()
    }
    elem <- "row.ddr.left"
    if (elem %in% elements && !is.null(row.ddr))
    {
        if (attr(row.ddr, "members") != nrow(exprs))
            warning("the number of row dendrogram leaves is not equal to the number of rows of exprs")

        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2]), clip="off"))
        graphics::par(new=TRUE, fig=gridFIG(), mar=c(0,0,0,0))
        graphics::plot(row.ddr, horiz=TRUE, leaflab="none", xaxt="n", yaxt="n", xaxs="r", yaxs="i")
        upViewport()
    }
    elem <- "row.ddr.right"
    if (elem %in% elements && !is.null(row.ddr))
    {
        if (attr(row.ddr, "members") != nrow(exprs))
            warning("the number of row dendrogram leaves is not equal to the number of rows of exprs")

        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))
        graphics::par(new=TRUE, fig=gridFIG(), mar=c(0,0,0,0))
        graphics::plot(row.ddr, horiz=TRUE, leaflab="none", xaxt="n", yaxt="n", xaxs="r", yaxs="i",
             xlim=c(0, attr(row.ddr, "height")))
        upViewport()
    }
    elem <- "col.ddr.up"
    if (elem %in% elements && !is.null(col.ddr))
    {
        if (attr(col.ddr, "members") != ncol(exprs))
            warning("the number of column dendrogram leaves is not equal to the number of columns of exprs")

        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))
        graphics::par(new=TRUE, fig=gridFIG(), mar=c(0,0,0,0))
        graphics::plot(col.ddr, horiz=FALSE, leaflab="none", xaxt="n", yaxt="n", xaxs="i", yaxs="r")
        upViewport()
    }
    elem <- "col.ddr.down"
    if (elem %in% elements && !is.null(col.ddr))
    {
        if (attr(col.ddr, "members") != ncol(exprs))
            warning("the number of column dendrogram leaves is not equal to the number of columns of exprs")

        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))
        graphics::par(new=TRUE, fig=gridFIG(), mar=c(0,0,0,0))
        graphics::plot(col.ddr, horiz=FALSE, leaflab="none", xaxt="n", yaxt="n", xaxs="i", yaxs="r",
             ylim=c(attr(col.ddr, "height"), 0))
        upViewport()
    }
    elem <- "clinical"
    if (elem %in% elements && !is.null(clinical))
    {
        if (ncol(clinical) != ncol(exprs))
            warning("the size of clinical does not correspond to the number of columns of exprs")

        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))
        heatmap.clinical(clinical, clinical.na.color)
        upViewport()
    }
    elem <- "clinical.labels.ljust"
    if (elem %in% elements && !is.null(clinical.labels))
    {
        if (length(clinical.labels) != nrow(clinical))
            warning("the number of clinical labels is not equal to the number of variables in clinical")

        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))

        clinical.labels[is.na(clinical.labels)] <- ""
        heatmap.labels(labels=clinical.labels, type="row.labels", just="left")
        upViewport()
    }
    elem <- "clinical.labels.rjust"
    if (elem %in% elements && !is.null(clinical.labels))
    {
        if (length(clinical.labels) != nrow(clinical))
            warning("the number of clinical labels is not equal to the number of variables in clinical")

        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))

        clinical.labels[is.na(clinical.labels)] <- ""
        heatmap.labels(labels=clinical.labels, type="row.labels", just="right")
        upViewport()
    }
    elem <- "title"
    if (elem %in% elements && !is.null(title))
    {
        idx <- which(layout.mat == elem, arr.ind=TRUE)
        pushViewport(viewport(name=elem,
                              layout.pos.row=unique(idx[,1]),
                              layout.pos.col=unique(idx[,2])))

        grid.text(title, gp=gpar(cex=0.9 * heatmap.labels.cex(title)))
        upViewport()
    }

    ## restore all standard R graphical parameters
    graphics::par(old.par)
    upViewport()


    row.ind <- 1:nrow(exprs)
    col.ind <- 1:ncol(exprs)
    if (!is.null(row.ddr)) row.ind <- order.dendrogram(rev(row.ddr))
    if (!is.null(col.ddr)) col.ind <- order.dendrogram(col.ddr)

    return(invisible(list(layout.mat=layout.mat,
                          widths=widths,
                          heights=heights,
                          color.scheme=color.scheme,
                          key.min=key.min,
                          key.max=key.max,
                          na.color=na.color,
                          scale=scale,
                          key.line.color=key.line.color,
                          row.labels=row.labels,
                          col.labels=col.labels,
                          row.clust=row.clust,
                          row.clust.dist.fun=row.clust.dist.fun,
                          row.clust.fun=row.clust.fun,
                          col.clust=col.clust,
                          col.clust.dist.fun=col.clust.dist.fun,
                          col.clust.fun=col.clust.fun,
                          clinical=clinical,
                          clinical.labels=clinical.labels,
                          clinical.na.color=clinical.na.color,
                          title=title,
                          row.ind=row.ind, col.ind=col.ind,
                          row.ddr=row.ddr, col.ddr=col.ddr)))
}
