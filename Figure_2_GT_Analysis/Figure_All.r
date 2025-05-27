suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(gridExtra)
  library(mclust)
  library(pheatmap)
  library(SpatialExperiment)
  library(scran)
  library(ggspavis)
  library(ggrepel)
  library(dbscan)
  library(scater)
  library(tibble)
  library(limma)
  library(bluster)
  library(reshape2)
  library(spatialLIBD)
})

######################### Figure 2A #########################
# Load data
spe <- fetch_data("spe")
spe <- spe[, spe$sample_id == "151673"]

# Get spe color
sce_layer <- fetch_data("sce_layer")
col_pal <- get_colors(libd_layer_colors, sce_layer$layer_guess_reordered_short)

# Using existing color table
p_2A <- vis_clus(spe = spe,
                 sampleid = "151673",
                 clustervar = "layer_guess_reordered_short",
                 colors = col_pal,
                 image_id = "lowres") + theme(legend.position = "none")


######################### Figure 2B #########################
# Load data

# QUALITY CONTROL (QC)
# subset to keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
# calculate per-spot QC metrics
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
# select QC thresholds
qc_lib_size <- colData(spe)$sum < 600
qc_detected <- colData(spe)$detected < 400
qc_mito <- colData(spe)$subsets_mito_percent > 28
qc_cell_count <- colData(spe)$cell_count > 10
# combined set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
colData(spe)$discard <- discard
# filter low-quality spots
spe <- spe[, !colData(spe)$discard]

# Normalization
spe <- logNormCounts(spe)

spe <- spe[!is_mito, ]

spe$ground_truth <- spe$layer_guess_reordered_short
keep <- !is.na(spe$ground_truth)
spe <- spe[,keep]

# fit mean-variance relationship
dec <- modelGeneVar(spe)

# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)

# Run PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top_hvgs)

# Get the border spots
knn <- kNN(colData(spe)[,c("array_row","array_col")], 4)
knn$id[knn$dist > 1.5] <- NA
knn$dist[knn$dist > 1.5] <- NA

nei_classes <- apply(knn$id, 1, function(u) {
  id <- na.omit(u)
  class <- spe$ground_truth[id]
  length(table(as.character(class)))
})

table(nei_classes)
spe$border <- nei_classes

############### Plot Figure 2A (right) ###############

d <- as.data.frame(cbind(colData(spe_sub), SpatialExperiment::spatialCoords(spe_sub)), optional = TRUE)
clustervar <- "ground_truth"
sampleid = unique(spe$sample_id)[1]
colors <- get_colors(libd_layer_colors, sce_layer$layer_guess_reordered_short)
image_id <- "lowres"
point_size = 2
na_color = "#CCCCCC00"

d$border <- ifelse(d$border==1, "internal","border")

pxl_row_in_fullres <- pxl_col_in_fullres <- key <- NULL
# stopifnot(all(c("pxl_col_in_fullres", "pxl_row_in_fullres", "key") %in% colnames(d)))

img <- SpatialExperiment::imgRaster(spe, sample_id = sampleid, image_id = image_id)

frame_lims <- spatialLIBD::frame_limits(spe, sampleid = sampleid, image_id = image_id)
img <- img[frame_lims$y_min:frame_lims$y_max, frame_lims$x_min:frame_lims$x_max]
adjust <- list(x = frame_lims$x_min, y = frame_lims$y_min)

p <- ggplot(
    d,
    aes(
        x = pxl_col_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sampleid, image_id = image_id) - adjust$x,
        y = pxl_row_in_fullres * SpatialExperiment::scaleFactors(spe, sample_id = sampleid, image_id = image_id) - adjust$y,
        fill = factor(!!sym(clustervar)),
        key = key,
        color = border
    )
)

grob <- grid::rasterGrob(img,
        width = grid::unit(1, "npc"),
        height = grid::unit(1, "npc")
    )
p <- p + geom_spatial(
        data = tibble::tibble(grob = list(grob)),
        aes(grob = grob),
        x = 0.5,
        y = 0.5
    )

p <- p +
    geom_point(
        pch = 21,
        size = point_size,
        stroke = 0.5
    ) +
    coord_fixed(expand = FALSE) +
    scale_fill_manual(values = colors, na.value = na_color) +
    scale_colour_manual(values = c("border" = "black", "internal"="#CCCCCC00")) + 
    xlim(0, ncol(img)) +
    ylim(nrow(img), 0) +
    xlab("") + ylab("") +
    labs(fill = NULL) +
    theme_set(theme_bw(base_size = 20)) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.box.spacing = unit(0, "pt")
    )

p_2A <- p

# Formulate plotting dataframe
cd <- colData(spe)
df <- data.frame(cd, 
                 PCA = reducedDim(spe, "PCA")[,1:3])

df$border <- ifelse(df$border==1, "internal","border")
df$`ground truth` <- df$ground_truth

# write.table(df, "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Cluster_Benchmark/Draft_Fig2/Figure_Dataframe/Figure2B.tsv", sep="\t")

############### Plot Figure 2B (right) ###############
p_2B <- ggplot(df %>% filter(ground_truth != "NA") %>%
         arrange(border),
       aes(PCA.PC1, PCA.PC2, fill = `ground truth`)) +
  geom_point(size = 4, pch = 21, stroke = 0) +
  scale_fill_manual(values = get_colors(libd_layer_colors, sce_layer$layer_guess_reordered_short)) +
  geom_point(data = df %>% filter(border=="border"),
             mapping = aes(colour = border), size = 4, pch = 21) +
  scale_colour_manual(values = c("border" = "black")) + 
  theme_classic() +
  labs(y="PC2", x="PC1") + 
  theme(axis.ticks = element_blank(), 
        axis.text = element_blank(),
        legend.position="inside",
        legend.position.inside = c(0.92, 0.75),
        legend.box.background = element_rect(colour = "black"))

######################### Figure 2C #########################
# Common function: plot_clusters
plot_clusters <- function(label, data, color_palette, ps=0.02, legend_position="none"){
  color_palette <- color_palette[1:length(unique(data[[label]]))] %>% setNames(unique(data[[label]]))
  con_p <- ggplot(data, aes(x, y, colour=as.factor(.data[[label]]))) + 
    geom_point(size = ps) + #%>% ggrastr::rasterise() +
    theme_classic(base_size=12) +
    theme(legend.position = legend_position, 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line=element_blank()) +
    # xlab("") + ylab("") +
    labs(colour="label") + 
    scale_colour_manual(values = color_palette) +
    scale_x_continuous(expand = expansion(0,1)) +
    scale_y_continuous(expand = expansion(0,1)) +
    ggtitle(label) + coord_fixed()
  return(con_p)
}

# Figure 2C: the cellcharter plot
XMBdata <- read.table("data/XMB/combined_results.tsv")

colP <- c("CellCharter_5", "CellCharter_7", "CellCharter_9", "CellCharter_11")

palette <- c("#75E673", "#8ABDDB", "#DF5782", "#D96AD0", "#D0D4BB", "#9A40E6", "#DC9560",
             "#D3DB6A", "#807DD4", "#D9A9CB", "#7CDEC2")

p_list <- lapply(colP, \(cp) {
  plt <- plot_clusters(label=cp, data=XMBdata, color_palette=palette, ps=0.001) + 
          geom_rect(aes(xmin = 4150, xmax = 6150, ymin = 5200, ymax = 6800),
                    fill = "transparent", linetype = "dashed", 
                    color = "blue", linewidth = 1) + ggtitle("")
  return(plt)
})

# Refine the first 5 cluster one
palette_5 <- c("#75E673", "#DF5782", "#8ABDDB", "#D0D4BB", "#9A40E6")
p_list[[1]] <- plot_clusters(label="CellCharter_5", data=XMBdata, color_palette=palette_5, ps=0.001) + 
          geom_rect(aes(xmin = 4150, xmax = 6150, ymin = 5200, ymax = 6800),
                    fill = "transparent", linetype = "dashed", 
                    color = "blue", linewidth = 1) + ggtitle("")

p_2C <- arrangeGrob(grobs=p_list, nrow=1)

######################### Figure 2D #########################
palette <- c("#DA8CBC", "#82E06B", "#CDBFD3", "#9486DB", "#D1D5B4", "#73BAD5", "#DCD263",
            "#DF61C9", "#D87B65", "#81E3C4", "#9F42DD")
p_2D <- plot_clusters("label", XMBdata, color_palette = palette, ps=0.2, legend_position = "bottom") + 
                    geom_rect(aes(xmin = 4150, xmax = 6150, ymin = 5200, ymax = 6800),
                    fill = "transparent", linetype = "dashed", 
                    color = "blue", linewidth = 1) +
                guides(color = guide_legend(override.aes = list(size = 3))) 

######################### Figure 2F #########################
Sterodata <- read.table("data/StereoSeq_liver/combined_results.tsv")

p_2F <- plot_clusters("label", Sterodata, color_palette = palette, ps=0.1, legend_position = "bottom") + 
                guides(color = guide_legend(override.aes = list(size = 3))) 

######################### Figure 2H #########################
VisiumData <- read.table("data/Visium_BC/combined_results.tsv")

p_2H <- plot_clusters("label", VisiumData, color_palette = palette, ps=1, legend_position = "bottom") + 
                guides(color = guide_legend(override.aes = list(size = 3))) 

######################### Figure 2EGI #########################
label_ARI_XMB <- read.table("data/XMB/label_ARI.tsv")
label_ARI_stereoseq <- read.table("data/StereoSeq_liver/label_ARI.tsv")
label_ARI_newlab <- read.table("data/StereoSeq_NewLab/label_ARI.tsv")
label_ARI_BC <- read.table("data/Visium_BC/label_ARI.tsv")
label_ARI_BC <- label_ARI_BC[label_ARI_BC$cluster_number < 8, ]
label_ARI_stereoseq <- label_ARI_stereoseq[label_ARI_stereoseq$cluster_number <= 9, ]
label_ARI_newlab <- label_ARI_newlab[label_ARI_newlab$cluster_number <= 9, ]

#label_ARI_selected$Dataset <- "Selected_StereoSeq"
label_ARI_BC$Dataset <- "Visium"
label_ARI_XMB$Dataset <- "Xenium"
label_ARI_stereoseq$Dataset <- "StereoSeq"
label_ARI_newlab$Dataset <- "StereoSeq"

label_ARI_BC$Label <- "Original"
label_ARI_XMB$Label <- "Original"
label_ARI_stereoseq$Label <- "Original"
label_ARI_newlab$Label <- "Aggregated"

ARI_df <- rbind(label_ARI_XMB, label_ARI_BC, label_ARI_stereoseq, label_ARI_newlab)
ARI_df <- ARI_df[!is.na(ARI_df$cluster_number), ]
#ARI_df <- ARI_df[!(ARI_df$method == "SpaceFlow" & ARI_df$Dataset == "StereoSeq"), ]
# Remove some results

vline_df <- data.frame(Dataset = c("Xenium", "Visium", "StereoSeq", "StereoSeq"), GT = c(11, 4, 3, 9))

ari_palette <- c("#00A08A", "#FF0000", "#F2AD00",
                 "#F98400", "#5BBCD6") %>% 
                 setNames(c("CellCharter", "SCANIT", "STAGATE", "SEDR", "SpaceFlow"))
linetype_palette <- c("solid", "dashed") %>% setNames(c("Original", "Aggregated"))

p_2EGI <- ggplot(ARI_df, aes(x = cluster_number, y = ARI, colour = method)) + 
            scale_x_continuous(breaks = unique(ARI_df$cluster_number)) + 
            geom_line(aes(linetype=Label), linewidth=1) +
            geom_point(shape = 15, size = 3) + 
            scale_color_manual(values = ari_palette) + 
            scale_linetype_manual(values = linetype_palette) + 
            theme_classic(base_size=12) + 
            facet_grid(.~factor(Dataset, c("Xenium", "Visium", "StereoSeq")), scales = "free_x") + 
            geom_vline(data = vline_df, aes(xintercept = GT),
                       linetype = "dashed", color = "blue",
                       linewidth = 0.5) +
            xlab("Number of clusters") + 
            theme(legend.position = "bottom", 
                  strip.background = element_blank()) +
            geom_text(data = vline_df, aes(x=GT, y = 0.05, label="\n GT"), colour="blue", angle=90, size=4)

######################### Assemble everything together #########################

ggsave("All_Subgraphs/2A.pdf", p_2A, width = 8, height = 8)
ggsave("All_Subgraphs/2B.pdf", p_2B, width = 10, height = 6)
ggsave("All_Subgraphs/2C.pdf", p_2C, width = 24, height = 6)
ggsave("All_Subgraphs/2D.pdf", p_2D, width = 6, height = 6)
ggsave("All_Subgraphs/2F.pdf", p_2F, width = 6, height = 6)
ggsave("All_Subgraphs/2H.pdf", p_2E, width = 6, height = 6)
ggsave("All_Subgraphs/2EGI.pdf", p_2EGI, width = 13, height = 6)

