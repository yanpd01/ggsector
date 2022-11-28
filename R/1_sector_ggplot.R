### legend -----------------------------------------------

#' draw_key_sector
#'
#' @rdname draw_key_sector
#'
#' @inheritParams ggplot2::draw_key_polygon
#' @return ggplot legend
#' @export
draw_key_sector <- function(data, params, size) {
    if (getOption("debug_sector_legend", FALSE)) {
        # print(head(data))
        # print(head(size))
    }
    # grob
    sectorGrob(
        r = 0.5,
        vp = viewport(0.25, 0.25),
        gp = gpar(
            col = data$colour %||% NA,
            fill = alpha(data$fill %||% "grey20", data$alpha)
        )
    )
}


### geom  -----------------------------------------------
# Draw in vector form

#' @rdname geom_sector
#' @format NULL
#' @usage NULL
#' @export
GeomSectorPanel <- ggproto(
    "GeomSectorPanel",
    Geom,
    required_aes = c("x", "y", "theta"),
    non_missing_aes = c("r", "start", "r_start", "type"),
    default_aes = aes(
        r = 0.45,
        start = 0,
        r_start = 0,
        type = "percent",
        colour = "black",
        fill = "transparent",
        alpha = NA
    ),
    setup_data = function(data, params) {
        data$width_raw <- resolution(data$x, FALSE)
        data$height_raw <- resolution(data$y, FALSE)
        return(data)
    },
    draw_panel = function(data, panel_params, coord, na.rm = FALSE) {
        ## coord transform
        coords <- coord$transform(data, panel_params)
        ## Calculate the range of the radius
        tmp_width <- diff(scales::rescale(
            c(0, min(coords$width_raw)),
            from = panel_params$x.range
        ))
        tmp_height <- diff(scales::rescale(
            c(0, min(coords$height_raw)),
            from = panel_params$y.range
        ))
        tmp_hw <- min(tmp_width, tmp_height)
        ## grob
        sectorGrob(
            x = coords$x,
            y = coords$y,
            theta = coords$theta,
            r = coords$r * tmp_hw,
            start = coords$start,
            r_start = coords$r_start * tmp_hw,
            type = coords$type,
            vp = NULL,
            gp = gpar(
                col = alpha(coords$colour, coords$alpha),
                fill = alpha(coords$fill, coords$alpha)
            )
        )
    },
    draw_key = draw_key_sector
)


# Draw in a single form

#' @rdname geom_sector
#' @format NULL
#' @usage NULL
#' @export
GeomSectorIndividual <- ggproto(
    "GeomSectorIndividual",
    Geom,
    required_aes = c("x", "y", "theta"),
    non_missing_aes = c("r", "start", "r_start", "type"),
    default_aes = aes(
        # size = 0.9,
        r = 0.45,
        start = 0,
        r_start = 0,
        type = "percent",
        colour = "black",
        fill = "transparent",
        alpha = NA
    ),
    setup_data = function(data, params) {
        data$width_raw <- resolution(data$x, FALSE)
        data$height_raw <- resolution(data$y, FALSE)
        return(data)
    },
    draw_group = function(data, panel_params, coord, na.rm = FALSE) {
        ## coord trans
        coords <- coord$transform(data, panel_params)
        ## Calculate the range of the radius
        tmp_width <- diff(scales::rescale(
            c(0, min(coords$width_raw)),
            from = panel_params$x.range
        ))
        tmp_height <- diff(scales::rescale(
            c(0, min(coords$height_raw)),
            from = panel_params$y.range
        ))
        tmp_hw <- min(tmp_width, tmp_height)
        # grob
        sectorGrob(
            x = 0.5,
            y = 0.5,
            theta = coords$theta,
            r = coords$r,
            start = coords$start,
            r_start = coords$r_start,
            type = coords$type,
            vp = viewport(
                x = coords$x,
                y = coords$y,
                # default.units = "native",
                width = unit(tmp_hw, "snpc"),
                height = unit(tmp_hw, "snpc")
            ),
            gp = gpar(
                col = alpha(coords$colour, coords$alpha),
                fill = alpha(coords$fill, coords$alpha)
            )
        )
    },
    draw_key = draw_key_sector
)

### layer -----------------

#' ggplot sector
#'
#' Draw sector with ggplot2.
#'
#' When "individual=FALSE", draw very quickly with a vector form,
#' when "individual=TRUE", draw individually at a slower speed.
#'
#' The required parameters in mapping are "x", "y", "theta", and
#' the additional modifiable parameters are "r", "start", "r_start", "type", "colour", "fill".
#'
#' When r = 0.5, the sector just fills the square.
#'
#' For details, please check the [grid.sector()].
#'
#' For more details, please type `vignette("ggsector")`.
#'
#' @rdname geom_sector
#'
#' @inheritParams ggplot2::geom_point
#' @param individual Logical, default is TRUE.
#' When "individual=FALSE", draw very quickly with a vector form,
#' when "individual=TRUE", draw individually at a slower speed.
#' Vector form cannot control the deformation of the sector, you need to add coord_fixed()
#' @param verbose Logical, default is TRUE. Whether to display reminder information.
#' @return ggplot object
#'
#' @examples
#' \donttest{
#' # prepare
#' library(ggsector)
#' library(reshape2)
#' df <- cor(mtcars)[1:3, 1:5] %>%
#'     abs() %>%
#'     melt(varnames = c("x", "y"))
#'
#' # individual
#' ggplot(df) +
#'     geom_sector(aes(x, y), theta = 75, fill = 2, individual = FALSE) +
#'     theme_bw() +
#'     theme(axis.title = element_blank())
#' ggplot(df) +
#'     geom_sector(aes(x, y), theta = 75, fill = 2, individual = FALSE) +
#'     coord_fixed(ratio = 3 / 5) +
#'     theme_bw() +
#'     theme(axis.title = element_blank())
#' ggplot(df) +
#'     geom_sector(aes(x, y), theta = 75, fill = 2, individual = TRUE) +
#'     theme_bw() +
#'     theme(axis.title = element_blank())
#'
#' # theta
#' ggplot(df) +
#'     ## type = "percent", theta = 0-100
#'     geom_sector(
#'         aes(y, x, theta = value * 100),
#'         type = "percent",
#'         color = "blue",
#'         individual = TRUE
#'     ) +
#'     ## type = "degree", theta = 0-360
#'     geom_sector(
#'         aes(y, x, theta = value * 360),
#'         type = "degree",
#'         color = "red",
#'         alpha = 0.5,
#'         individual = TRUE
#'     ) +
#'     theme_bw() +
#'     theme(axis.title = element_blank())
#'
#' # r
#' ggplot(df) +
#'     geom_sector(
#'         aes(y, x, theta = value * 100),
#'         r = rep(c(0.15, 0.3, 0.45), 5),
#'         fill = 2,
#'         individual = TRUE
#'     ) +
#'     theme_bw() +
#'     theme(axis.title = element_blank())
#'
#' # start
#' ggplot(df) +
#'     geom_sector(
#'         aes(y, x, theta = value * 100),
#'         start = rep(c(60, 40, 20), 5),
#'         fill = 2,
#'         individual = TRUE
#'     ) +
#'     theme_bw() +
#'     theme(axis.title = element_blank())
#'
#' # r_start
#' ggplot(df) +
#'     geom_sector(
#'         aes(y, x, theta = value * 100),
#'         r_start = rep(c(0.15, 0.25, 0.35), 5),
#'         fill = 2,
#'         individual = TRUE
#'     ) +
#'     theme_bw() +
#'     theme(axis.title = element_blank())
#' }
#'
#' @export
geom_sector <- function(mapping = NULL, data = NULL,
                        stat = "identity", position = "identity",
                        ...,
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE,
                        individual = TRUE,
                        verbose = TRUE
                        #
) {
    if (individual) {
        layer(
            data = data,
            mapping = mapping,
            stat = stat,
            geom = GeomSectorIndividual,
            position = position,
            show.legend = show.legend,
            inherit.aes = inherit.aes,
            params = list(
                na.rm = na.rm,
                ...
            )
        )
    } else {
        if (verbose) message("For better display effect, please add
`coord_fixed(ratio = length(unique(x)) / length(unique(y)))`")
        layer(
            data = data,
            mapping = mapping,
            stat = stat,
            geom = GeomSectorPanel,
            position = position,
            show.legend = show.legend,
            inherit.aes = inherit.aes,
            params = list(
                na.rm = na.rm,
                ...
            )
        )
    }
}
