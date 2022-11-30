# # (x = 0.5,
# #                                y = 0.5,
# #                                theta = 25,
# #                                r = 0.5,
# #                                start = 0,
# #                                r_start = 0,
# #                                type = "percent",
# #                                ratio = 1,
# #                                group
# #                                #
# # )


# pie_df <- function(x = 0.5,
#                    y = 0.5,
#                    terms,
#                    values = 1,
#                    highlight = FALSE,
#                    highlight_ratio = 0.1,
#                    r = 0.5,
#                    start = 0,
#                    r_start = 0,
#                    type = "percent",
#                    ratio = 1
#                    #
# ) {
#     df_in <- data.frame(
#         x = x,
#         y = y,
#         values = values,
#         highlight = highlight,
#         highlight_ratio = highlight_ratio,
#         r = r,
#         start = start,
#         r_start = r_start,
#         type = type,
#         ratio = ratio
#     )


#     v_sum <- sum(values)
#     pct <- values / v_sum
#     theta <- round(pct * n)
#     start_in <- cumsum(c(start, theta[-length(theta)]))
#     df_out <- sector_df_multiple(
#         x = x,
#         y = y,
#         theta = theta,
#         r = r,
#         start = start_in,
#         r_start = r_start,
#         type = type,
#         ratio = ratio
#     )
#     len <- nrow(df_out) + 1
#     rank_x <- as.numeric(factor(df_out$x)) * len^2
#     rank_y <- as.numeric(factor(df_out$y)) * len
#     df_out$group <- rank_x + rank_y + df_out$group
#     return(df_out)
# }




# # # 突出某个sector
# # # 从xy的坐标点上动手
# # # group_new=paste0（x，y，group）

# # # 颜色不用考虑，ggplot会帮我们搞定的
# # # 在grob中研究level，numic决定颜色，color就设置透明色或者黑色

# grid.newpage()
# grid.sector(start = c(0, 25, 50, 75), theta = 25, gp = gpar(fill = 11:14))

# sector_df_multiple(start = c(0, 25, 50, 75))$group
