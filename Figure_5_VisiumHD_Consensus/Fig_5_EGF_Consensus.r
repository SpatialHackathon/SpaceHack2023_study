# Change this to your github repo dir
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(cowplot)
  library(clue)
  library(khroma)
  library(scran)
  library(limma)
  library(tibble)
  library(readr)
  library(ggrepel)
  library(mclust)
  library(pheatmap)
  library(fastDummies)
  library(randomcoloR)
  library(reshape2)
  library(poem)
  library(gridExtra)
  library(optparse)
  library(yaml)
})
source("utils.R")

################## Load dataset 
df <- read.delim("data/combiend_results.tsv", header=TRUE, row.names=1)

cluster_instance <- grep("_label$", colnames(df), value=TRUE)
all_clusters <- df[, cluster_instance]

################## Define plotting function and parameters
set.seed(2025)
palette <- distinctColorPalette(k=10)

plot_clusters <- function(label, data, color_palette, ps=0.01){
  color_palette <- color_palette[1:length(unique(data[[label]]))] %>% setNames(unique(data[[label]]))
  con_p <- ggplot(data, aes(col, row, colour=as.factor(.data[[label]]))) + 
    geom_point(size = ps) +
    theme_classic() +
    theme(legend.position = "none", 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
    # xlab("") + ylab("") +
    labs(colour="label") + 
    scale_colour_manual(values = color_palette) +
    scale_x_continuous(expand = expansion(0,1)) +
    scale_y_continuous(expand = expansion(0,1)) +
    ggtitle(label) + coord_fixed()
  return(con_p)
}

###### 4E Pairwise ARI results
aris_all <- calc_aris(all_clusters)
greens <- colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Greens"))(100)

ph_all <- pheatmap(aris_all, 
                   fontsize = 8, 
                   treeheight_row = 0,
                   treeheight_col = 20,
                   main='Cross method ARI', 
                   labels_row = sub("_.*$", "", rownames(aris_all)),
                   labels_col = sub("_.*$", "", colnames(aris_all)),
                   color = greens)

ggsave("VisiumHD_ARI.pdf", ph_all[[4]], width = 4, height = 4)

###### Consensus Calling

# Select base clusterings based on tree cuts
hl <- cutree(ph_all$tree_row, 3)
meths <- names(hl)[hl == as.numeric(names(which.max(table(hl))))]
nclust <- lapply(meths, \(x) length(unique(df[[x]]))) %>% unlist()
if (!all(nclust = max(nclust))){warning("Result instances have different nclust, use largest one for consensus.")}

# Subset
base_clusterings <- all_clusters %>%
                    dplyr::select(meths)

# Find reference for alignment
ses <- apply(base_clusterings, 2, function(u) spot_entropy(df[,c("row","col")], u))
mses <- colMeans(ses)
base_clusterings <- align_classes(base_clusterings, names(which.min(mses)))

# consensus calling
df[[sprintf("consensus_lca_%s", nclust)]] <- diceR:::LCA(base_clusterings,
                                            is.relabelled = FALSE, seed = 100)

############### Figure 5G Consensus plot
p_con <- plot_clusters(data=df, 
                       label=sprintf("consensus_lca_%s", nclust), 
                       color_palette=palette)

ggsave("5G_consensus.pdf", plot = p_con, width=8, height=8)

############### Figure 5F CME plot
lv <- levels(base_clusterings[[meths[1]]])
ps <- apply(base_clusterings, 1, function(u) table(factor(u, levels=lv))) %>% t
df$ent_bc <- apply(ps, 1, calc_entropy)

write.table(df[,c(meths, "ent_bc")], "data/con_ent_results.tsv", sep="\t")

cme_p <- ggplot(df, aes(col, row, colour=ent_bc)) +
  geom_point(size = 0.01) +
  scale_colour_gradient(low="gray90", high="deeppink4") +
  theme_classic() +
  theme(legend.position = "right", 
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_blank(),
      axis.text.y=element_blank(),
      axis.ticks.y=element_blank()) +
  scale_x_continuous(expand = expansion(0,1)) +
  scale_y_continuous(expand = expansion(0,1)) +
  coord_fixed() +
  labs(subtitle = sprintf("mean entropy: %.4f", mean(df[["ent_bc"]])), 
       colour = "CME")

ggsave("5F_CME.pdf", plot = cme_p, width=8, height=8)

###### Supplementary figure: Base clustering results
p_vhd <- lapply(names(n_clust), 
                plot_clusters, data=df, color_palette=palette, ps=0.01)

ggsave("supp_VisiumHD_BCs.pdf", 
       plot = marrangeGrob(p_vhd, ncol=4, nrow=3),
       width=8, height=8)
