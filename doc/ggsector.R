## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
    collapse = TRUE,
    comment = "#>",
    fig.width = 8,
    fig.height = 8
)

## ----setup, message = FALSE---------------------------------------------------
library(magrittr)
library(grid)
library(ComplexHeatmap)
library(Seurat)
library(ggplot2)
library(ggsector)

## ---- fig.width=3, fig.height=3-----------------------------------------------
tmp_df <- sector_df(x = 0.5, y = 0.5, theta = 25, r = 0.4, start = 0, r_start = 0)
tmp_df
grid.newpage()
grid.polygon(
    tmp_df$x, tmp_df$y,
    vp = viewport(height = unit(1, "snpc"), width = unit(1, "snpc"))
)

## ---- fig.width=3, fig.height=3-----------------------------------------------
tmp_df <- sector_df(x = 0.5, y = 0.5, theta = 25, r = 0.4, start = 50, r_start = 0.2)
tmp_df
grid.newpage()
grid.polygon(
    tmp_df$x, tmp_df$y,
    vp = viewport(height = unit(1, "snpc"), width = unit(1, "snpc"))
)

## ---- fig.width=3, fig.height=3-----------------------------------------------
tmp_df <- sector_df(
    x = 0.5, y = 0.5, theta = 180, r = 0.4,
    start = 90, r_start = 0, type = "degree"
)
tmp_df
grid.newpage()
grid.polygon(
    tmp_df$x, tmp_df$y,
    vp = viewport(height = unit(1, "snpc"), width = unit(1, "snpc"))
)

## ---- fig.width=3, fig.height=3-----------------------------------------------
tmp_df <- sector_df(
    x = 0.5, y = 0.5, theta = 180, r = 0.4,
    start = 270, r_start = 0.2, type = "degree"
)
tmp_df
grid.newpage()
grid.polygon(
    tmp_df$x, tmp_df$y,
    vp = viewport(height = unit(1, "snpc"), width = unit(1, "snpc"))
)

## ---- fig.width=3, fig.height=3-----------------------------------------------
tmp_df <- sector_df_multiple(
    x = c(0.2, 0.5, 0.8),
    theta = c(25, 50, 75),
    r = 0.15,
    start = c(75, 50, 100),
    r_start = c(0, 0.05, 0.1),
    type = "percent"
)
tmp_df
grid.newpage()
grid.polygon(
    tmp_df$x,
    tmp_df$y,
    id = tmp_df$group,
    vp = viewport(height = unit(1, "snpc"), width = unit(1, "snpc")),
    gp = gpar(
        fill = 3:1, col = 1:3
    )
)

## ---- fig.width=5, fig.height=4-----------------------------------------------
grid.newpage()
gp <- sectorGrob(
    x = unit(c(3, 9, 15), "cm"),
    y = unit(c(5, 9, 15), "cm"),
    theta = c(90, 180, 270),
    r = 1,
    start = c(180, 180, 270),
    r_start = c(0.6, 0.3, 0),
    type = "degree",
    group = factor(1:3, levels = c(2, 3, 1)),
    gp = gpar(fill = c("green", "red", "blue"))
)
grid.draw(gp)

## ---- fig.width=5, fig.height=4-----------------------------------------------
grid.newpage()
grid.sector(
    x = c(0.1, 0.5, 0.9),
    y = c(0.9, 0.6, 0.1),
    theta = c(25, 50, 90),
    r = .1,
    start = c(25, 50, 100),
    r_start = c(0.06, 0.03, 0),
    type = "percent",
    group = factor(1:3, levels = c(2, 3, 1)),
    gp = gpar(col = c("green", "red", "blue"), fill = 2:4),
    default.units = "npc"
)

## -----------------------------------------------------------------------------
library(magrittr)
library(ComplexHeatmap)
t0 <- cor(mtcars) %>%
    set_colnames(paste("y_", colnames(.))) %>%
    set_rownames(paste("x_", rownames(.)))
mat <- abs(t0)
mat[1:5, 1:5]

## -----------------------------------------------------------------------------
set.seed(1)
Heatmap(
    mat,
    name = "vp",
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(
            x = x, y = y, width = width, height = height,
            gp = gpar(col = "grey", fill = NA)
        )
        grid.sector(
            theta = mat[i, j] * 100,
            r = 0.5,
            start = mat[i, j] * 100 * runif(1),
            r_start = mat[i, j] * 0.49 * runif(1),
            vp = viewport(x, y, width, height),
            gp = gpar(fill = fill, col = "transparent")
        )
    },
    width = unit(.7, "snpc"),
    height = unit(.7, "snpc")
)

## -----------------------------------------------------------------------------
# The default viewport locks the horizontal and vertical axes
# so that the sector does not deform, which needs to be removed here.
# The radius 'r' is half the min(length, width).
set.seed(2)
Heatmap(
    mat,
    name = "xy + r",
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(
            x = x, y = y, width = width, height = height,
            gp = gpar(col = "grey", fill = NA)
        )
        r <- as.numeric(min(width, height)) / 2
        grid.sector(
            x,
            y,
            theta = mat[i, j] * 100,
            r = r,
            start = mat[i, j] * 100 * runif(1),
            r_start = mat[i, j] * r * 0.9 * runif(1),
            vp = NULL,
            gp = gpar(fill = fill, col = "transparent")
        )
    },
    width = unit(.7, "snpc"),
    height = unit(.7, "snpc")
)

## -----------------------------------------------------------------------------
# The input matrix needs to be extracted with pindex(mat, i, j)
set.seed(3)
Heatmap(
    mat,
    name = "layer",
    rect_gp = gpar(type = "none"),
    layer_fun = function(j, i, x, y, width, height, fill) {
        grid.rect(
            x = x, y = y, width = width, height = height,
            gp = gpar(col = "grey", fill = NA)
        )
        r <- as.numeric(min(width, height)) / 2
        grid.sector(
            x,
            y,
            theta = pindex(mat, i, j) * 100,
            r = r,
            start = pindex(mat, i, j) * 100 * runif(nrow(mat) * ncol(mat)),
            r_start = pindex(mat, i, j) * r * 0.9 * runif(nrow(mat) * ncol(mat)),
            vp = NULL,
            gp = gpar(fill = fill, col = "transparent")
        )
    },
    width = unit(.7, "snpc"),
    height = unit(.7, "snpc")
)

## -----------------------------------------------------------------------------
set.seed(1)
t0 <- cor(mtcars) %>%
    set_colnames(paste("y_", colnames(.))) %>%
    set_rownames(paste("x_", rownames(.)))
t1 <- t0 %>%
    reshape2::melt()
t1$group <- sample(letters[1:5], nrow(t1), TRUE)
t1$r <- sample(c(0.24, 0.48), nrow(t1), TRUE)
t1$start <- runif(nrow(t1), max = 50)
t1$r_start <- sample(c(0, 0.08, 0.16), nrow(t1), TRUE)
t1$v <- abs(t1$value)
head(t1)

## -----------------------------------------------------------------------------
ggplot(t1) +
    aes(x = Var1, y = Var2, theta = v * 100) +
    geom_sector(
        aes(fill = group, r = r, start = start, r_start = r_start),
        color = "transparent"
    ) +
    theme_bw()

## -----------------------------------------------------------------------------
ggplot(t1) +
    aes(x = Var1, y = Var2, theta = v * 360) +
    geom_sector(
        aes(fill = group, r = r, start = start * 3.6, r_start = r_start),
        color = "transparent",
        type = "degree"
    ) +
    coord_fixed() +
    theme_bw()

## -----------------------------------------------------------------------------
ggplot(t1) +
    aes(x = Var1, y = Var2, theta = v * 100, r_start = r_start) +
    geom_sector(
        aes(fill = group, r = r, start = start),
        color = "transparent",
        individual = TRUE
    ) +
    theme_bw()

## ---- eval = FALSE------------------------------------------------------------
#  ## Download pbmc data from
#  # https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
#  pbmc.data <- Read10X(data.dir = "../filtered_gene_bc_matrices/hg19")
#  pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
#  pbmc <- NormalizeData(pbmc)
#  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#  pbmc <- ScaleData(pbmc, features = rownames(pbmc))
#  pbmc <- RunPCA(pbmc)
#  pbmc <- RunUMAP(pbmc, dim = 1:10)
#  pbmc <- FindNeighbors(pbmc, dims = 1:10)
#  pbmc <- FindClusters(pbmc, resolution = 0.5)
#  mks <- tibble::tribble(
#      ~type, ~marker,
#      "Naive CD4+ T", "IL7R,CCR7",
#      "CD14+ Mono", "CD14,LYZ",
#      "Memory CD4+", "IL7R,S100A4",
#      "B", "MS4A1",
#      "CD8+ T", "CD8A",
#      "FCGR3A+ Mono", "FCGR3A,MS4A7",
#      "NK", "GNLY,NKG7",
#      "DC", "FCER1A,CST3",
#      "Platelet", "PPBP",
#  ) %>%
#      tidyr::separate_rows(marker, sep = ", *") %>%
#      dplyr::distinct()
#  # Dotplot
#  DotPlot(pbmc, features = unique(mks$marker)) + coord_flip()
#  # contrast with DotPlot
#  SectorPlot(pbmc, c(mks$marker, "fsdd"), features_level = unique(rev(mks$marker)))
