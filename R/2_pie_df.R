


# sector_df_multiple <- function(x = 0.5,
#                                y = 0.5,
#                                theta = 25,
#                                r = 0.5,
#                                start = 0,
#                                r_start = 0,
#                                type = "percent",
#                                group
#                                #
# ) {
#     df_in <- data.frame(
#         x = x,
#         y = y,
#         theta = theta,
#         r = r,
#         start = start,
#         r_start = r_start,
#         type = type
#     )

#     if (missing(group)) group <- seq_len(nrow(df_in))
#     if (length(group) != nrow(df_in)) stop("Variables should be of uniform length")

#     out_list <- lapply(seq_len(nrow(df_in)), function(x) {
#         df_tmp <- do.call(sector_df, df_in[x, ])
#         df_tmp$group <- group[x]
#         return(df_tmp)
#     })
#     out_df <- do.call(rbind, out_list)
#     return(out_df)
# }

# pie_df <- function(x = 0.5,
#                    y = 0.5,
#                    theta,
#                    r = 0.5,
#                    start = 0,
#                    r_start = 0,
#                    type = "percent",
#                    group = 1:4
#                    #
# ) {

# }


# sector_df_multiple()

# # grid.newpage(); grid.sector(start = c(0, 25, 50, 75), theta = 25, gp = gpar(fill = 11:14))
