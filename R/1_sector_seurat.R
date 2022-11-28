
#' Draw sector for seurat object
#'
#' A better alternative to [Seurat::DotPlot()].
#' For more details, please type `vignette("ggsector")`.
#'
#' @rdname seurat
#' @param object Seurat object
#' @param features Input vector of genes list.
#' @param features_level Levels of genes list.
#' @param assay Specific assay to get data from or set data for; defaults to the default assay.
#' @param slot Specific assay data to get or set
#' @param group.by Column of metadata to group the cells by, default is Idents().
#' @param group_level Levels of group.
#' @param col_low Colours for low ends of the gradient.
#' @param col_mid Colour for mid point.
#' @param col_high Colours for high ends of the gradient.
#' @param col_midpoint The midpoint (in data value) of the diverging scale.
#' Defaults to quantile(exp, 0.5)
#'
#' @return ggplot
#'
#' @examples
#' \donttest{
#' ## Download pbmc data from
#' # https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
#' library(Seurat)
#' path <- paste0(tempdir(), "/pbmc3k.tar.gz")
#' file <- paste0(tempdir(), "/filtered_gene_bc_matrices/hg19")
#' download.file(
#'     "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
#'     path
#' )
#' untar(path, exdir = tempdir())
#' pbmc.data <- Read10X(data.dir = file)
#' pbmc <- CreateSeuratObject(
#'     counts = pbmc.data,
#'     project = "pbmc3k",
#'     min.cells = 3,
#'     min.features = 200
#' )
#' pbmc <- NormalizeData(pbmc)
#' pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
#' pbmc <- ScaleData(pbmc, features = rownames(pbmc))
#' pbmc <- RunPCA(pbmc)
#' pbmc <- RunUMAP(pbmc, dim = 1:10)
#' pbmc <- FindNeighbors(pbmc, dims = 1:10)
#' pbmc <- FindClusters(pbmc, resolution = 1)
#' pbmc <- FindClusters(pbmc, resolution = 0.5)
#' markers <- tibble::tribble(
#'     ~type, ~marker,
#'     "Naive CD4+ T", "IL7R,CCR7",
#'     "CD14+ Mono", "CD14,LYZ",
#'     "Memory CD4+", "IL7R,S100A4",
#'     "B", "MS4A1",
#'     "CD8+ T", "CD8A",
#'     "FCGR3A+ Mono", "FCGR3A,MS4A7",
#'     "NK", "GNLY,NKG7",
#'     "DC", "FCER1A,CST3",
#'     "Platelet", "PPBP",
#' ) %>%
#'     tidyr::separate_rows(marker, sep = ", *") %>%
#'     dplyr::distinct()
#' # Dotplot
#' DotPlot(pbmc, features = unique(markers$marker)) + coord_flip()
#' # contrast with DotPlot
#' SectorPlot(pbmc, markers$marker, features_level = unique(rev(markers$marker)))
#' SectorPlot(pbmc, markers$marker, group.by = "RNA_snn_res.1")
#' }
#'
#' @export
SectorPlot <- function(object,
                       features,
                       features_level,
                       assay,
                       slot = c("data", "scale.data", "counts"),
                       group.by,
                       group_level,
                       col_low = "blue",
                       col_mid = "white",
                       col_high = "red",
                       col_midpoint
                       #
) {
    ## define var
    . <- NULL
    ## para
    sob <- object
    if (missing(assay)) assay <- Seurat::DefaultAssay(sob)
    slot <- match.arg(slot)
    if (missing(group.by)) {
        group <- Seurat::Idents(sob)
    } else {
        group <- sob[[group.by]][[1]]
    }
    # select gene
    mk_gene <- intersect(features, rownames(sob))
    diff_gene <- setdiff(features, mk_gene)
    if (length(diff_gene)) warning(paste("These genes were not found in the Seurat object:\n", diff_gene))

    ## treat data
    df_sob <- Seurat::GetAssayData(sob, slot = slot, assay = assay)[mk_gene, ]
    df_sum <-
        Matrix::t(df_sob) %>%
        tibble::as_tibble() %>%
        dplyr::bind_cols(group = group, .) %>%
        dplyr::group_by(group) %>%
        dplyr::summarise_all(list(exp = mean, pct = function(x) sum(x > 0) / length(x))) %>%
        tidyr::pivot_longer(
            -1,
            names_to = c("gene", ".value"),
            names_pattern = "(^.*)_(.{3}$)"
        )
    ## col and level
    if (missing(col_midpoint)) col_midpoint <- stats::quantile(df_sum$exp, 0.5)
    if (missing(features_level)) features_level <- mk_gene
    if (!missing(group_level)) {
        g_level <- unique(group_level)
        if (length(g_level) != length(group_level)) warning("Duplicate values in group_level")
        df_sum$group <- factor(df_sum$group, levels = g_level)
    }
    f_level <- unique(features_level)
    if (length(f_level) != length(features_level)) warning("Duplicate values in features_level")
    df_sum$gene <- factor(df_sum$gene, levels = rev(f_level))
    # ggplot
    ggplot(df_sum) +
        aes_string("group", "gene", theta = "pct * 100", fill = "exp") +
        geom_sector(verbose = FALSE) +
        coord_fixed(length(unique(group)) / length(unique(mk_gene))) +
        scale_fill_gradient2(low = col_low, mid = col_mid, high = col_high, midpoint = col_midpoint) +
        theme_classic()
}
