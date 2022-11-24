## sector df -------------------------
sector_df_100 <- data.frame(
    x = round(sin(seq(0, 4 * pi, length.out = 201)), 10),
    y = round(cos(seq(0, 4 * pi, length.out = 201)), 10)
)
sector_df_360 <- data.frame(
    x = round(sin(seq(0, 4 * pi, length.out = 721)), 10),
    y = round(cos(seq(0, 4 * pi, length.out = 721)), 10)
)


#' sector coordinates
#'
#' According to the input center position, radius and angle,
#' get the polygon coordinates of a sector.
#'
#' [sector_df()] Only one value can be passed in for each parameter,
#'  and a sector coordinate is returned.
#'
#' [sector_df_multiple()] Each parameter can pass in multiple values,
#'  and return multiple sector coordinates
#'
#' For more details, please type `vignette("ggsector")`.
#'
#' @rdname sector_df
#' @param x Numeric, the x-axis coordinate of the sector center.
#' @param y Numeric, the y-axis coordinate of the sector center.
#' @param theta Numeric, the angle of the sector,
#' if 'type = "percent"', the input is a percentage(0-100),
#' if 'type = "degree"', the input is an angle(0-360).
#' @param r Numeric, radius of the outer circle of the sector(0-0.5).
#' @param r_start Numeric, radius of the inner circle of the sector(0-r).
#' @param start Numeric, starting angle of sector.
#' @param type "percent" or "degree".
#'
#' @return coordinates of sector.
#'
#' @examples
#' \donttest{
#' ## coordinates of single sector
#' # type of percent, start = 0, r_start = 0
#' tmp_df <- sector_df(x = 0.5, y = 0.5, theta = 25, r = 0.4, start = 0, r_start = 0)
#' tmp_df
#' grid.newpage()
#' grid.polygon(
#'     tmp_df$x, tmp_df$y,
#'     vp = viewport(height = unit(1, "snpc"), width = unit(1, "snpc"))
#' )
#' # type of percent, start = 50, r_start = 0.2
#' tmp_df <- sector_df(x = 0.5, y = 0.5, theta = 25, r = 0.4, start = 50, r_start = 0.2)
#' tmp_df
#' grid.newpage()
#' grid.polygon(
#'     tmp_df$x, tmp_df$y,
#'     vp = viewport(height = unit(1, "snpc"), width = unit(1, "snpc"))
#' )
#'
#' # type of degree, start = 90, r_start = 0
#' tmp_df <- sector_df(
#'     x = 0.5, y = 0.5, theta = 180, r = 0.4,
#'     start = 90, r_start = 0, type = "degree"
#' )
#' tmp_df
#' grid.newpage()
#' grid.polygon(
#'     tmp_df$x, tmp_df$y,
#'     vp = viewport(height = unit(1, "snpc"), width = unit(1, "snpc"))
#' )
#' # type of degree, start = 180, r_start = 0.2
#' tmp_df <- sector_df(
#'     x = 0.5, y = 0.5, theta = 180, r = 0.4,
#'     start = 270, r_start = 0.2, type = "degree"
#' )
#' tmp_df
#' grid.newpage()
#' grid.polygon(
#'     tmp_df$x, tmp_df$y,
#'     vp = viewport(height = unit(1, "snpc"), width = unit(1, "snpc"))
#' )
#'
#' # Coordinates of Multiple Sectors
#' tmp_df <- sector_df_multiple(
#'     x = c(0.2, 0.5, 0.8),
#'     theta = c(25, 50, 75),
#'     r = 0.15,
#'     start = c(75, 50, 100),
#'     r_start = c(0, 0.05, 0.1),
#'     type = "percent"
#' )
#' tmp_df
#' grid.newpage()
#' grid.polygon(
#'     tmp_df$x,
#'     tmp_df$y,
#'     id = tmp_df$group,
#'     vp = viewport(height = unit(1, "snpc"), width = unit(1, "snpc")),
#'     gp = gpar(
#'         fill = 3:1, col = 1:3
#'     )
#' )
#' }
#' @export
sector_df <- function(x = 0.5,
                      y = 0.5,
                      theta = 25,
                      r = 0.5,
                      start = 0,
                      r_start = 0,
                      type = c("percent", "degree")
                      #
) {
    type <- match.arg(type)
    if (type == "percent") {
        if (theta < 0 || theta > 100) stop('The "percent" should be between [0-100]')
        if (start < 0 || start > 100) stop('The "start of percent" should be between [0-100]')
        if (r_start < 0 || r_start >= r) stop('The "r_start" should be between [0-r)')
        theta <- round(theta)
        df_in <- sector_df_100[(start + 1):(start + theta + 1), ] * r
        if (r_start != 0) {
            df_start <- sector_df_100[(start + 1):(start + theta + 1), ] * r_start
            tmp_x <- c(df_in$x + x, rev(df_start$x + x))
            tmp_y <- c(df_in$y + y, rev(df_start$y + y))
        } else {
            tmp_x <- c(x, df_in$x + x, x)
            tmp_y <- c(y, df_in$y + y, y)
        }
    }
    if (type == "degree") {
        if (theta < 0 || theta > 360) stop('The "degree" should be between [0-360]')
        if (start < 0 || start > 360) stop('The "start of degree" should be between [0-360]')
        if (r_start < 0 || r_start >= r) stop('The "r_start" should be between [0-r)')
        theta <- round(theta)
        df_in <- sector_df_360[(start + 1):(start + theta + 1), ] * r
        if (r_start != 0) {
            df_start <- sector_df_360[(start + 1):(start + theta + 1), ] * r_start
            tmp_x <- c(df_in$x + x, rev(df_start$x + x))
            tmp_y <- c(df_in$y + y, rev(df_start$y + y))
        } else {
            tmp_x <- c(x, df_in$x + x, x)
            tmp_y <- c(y, df_in$y + y, y)
        }
    }
    return(data.frame(x = tmp_x, y = tmp_y))
}


#' @rdname sector_df
#' @param group The id provided when drawing a sector with [grid::polygonGrob()].
#' If missing, it will be automatically added as a number.
#'
#' @return  coordinates of sectors.
#' @export
sector_df_multiple <- function(x = 0.5,
                               y = 0.5,
                               theta = 25,
                               r = 0.5,
                               start = 0,
                               r_start = 0,
                               type = "percent",
                               group
                               #
) {
    df_in <- data.frame(
        x = x,
        y = y,
        theta = theta,
        r = r,
        start = start,
        r_start = r_start,
        type = type
    )

    if (missing(group)) group <- seq_len(nrow(df_in))
    if (length(group) != nrow(df_in)) stop("Variables should be of uniform length")

    out_list <- lapply(seq_len(nrow(df_in)), function(x) {
        df_tmp <- do.call(sector_df, df_in[x, ])
        df_tmp$group <- group[x]
        return(df_tmp)
    })
    out_df <- do.call(rbind, out_list)
    return(out_df)
}
