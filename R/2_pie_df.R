


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
#     1
# }

# pie_df <- function(x = 0.5,
#                    y = 0.5,
#                    terms,
#                    values = 1,
#                    highlight = FALSE,
#                    highlight_ratio = 0.1,
#                    r = 0.5,
#                    start = 0,
#                    r_start = 0,
#                    type = "percent"
#                    #
# ) {
#     if (missing(terms)) terms <- seq_along(values)
#     if (type == "percent") {
#         n <- 100
#     } else if (type == "degree") {
#         n <- 360
#     } else {
#         n <- round(as.numeric(type))

#     }
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
#         group = paste("a", x, y, seq_along(values), sep = ".")
#     )
#     return(df_out)
# }
# pie_df()


# # 突出某个sector
# # 从xy的坐标点上动手
# # group_new=paste0（x，y，group）

# # 颜色不用考虑，ggplot会帮我们搞定的
# # 在grob中研究level，numic决定颜色，color就设置透明色或者黑色
