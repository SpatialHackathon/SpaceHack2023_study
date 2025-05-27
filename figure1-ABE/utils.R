color_palette_fig1A <- c(
  "original_paper" = "black",
  "benchmark_Hu2024" = "#E41A1C", 
  "scMEB" = "#377EB8", 
  "conST"  = "#4DAF4A", 
  "BANKSY" = "#984EA3", 
  "SpaceFlow"= "#FF7F00", 
  "DeepST"  = "#FFFF33", 
  "CellCharter" = "#A65628", 
  "STAGATE" = "#F781BF", 
  "SPICEMIX" = "#999999", 
  "DR-SC" = "#66C2A5", 
  "BASS" = "#FC8D62", 
  "SpaGCN" = "#8DA0CB", 
  "GraphST" = "#E78AC3",
  "BayesSpace" = "#360EC9"
)

# Keep method names consistent
recoding_map <- c(
  "precast" = "PRECAST",
  "SpiceMix" = "SPICEMIX",
  "seurat" = "Seurat",
  "DRSC" = "DR-SC",
  "SCANIT" = "SCAN-IT",
  "bass" = "BASS",
  "stardust" = "Stardust",
  "SC_MEB" = "scMEB",
  "meringue" = "MERINGUE",
  "maple" = "MAPLE",
  "deepST" = "DeepST",
  "spaGCN" = "SpaGCN"
)

## violin plots
theme_pub <- function(
    font_size = 14,
    font_family = "",
    line_size = .5,
    rel_small = 12/14,
    rel_tiny = 11/14,
    rel_large = 16/14
) {
  half_line <- font_size / 2
  small_size <- rel_small * font_size
  
  # work off of theme_grey just in case some new theme element comes along
  theme_grey(base_size = font_size, base_family = font_family) %+replace%
    theme(
      line              = element_line(colour = "black", size = line_size, linetype = 1, lineend = "butt"),
      rect              = element_rect(fill = NA, colour = NA, size = line_size, linetype = 1),
      text              = element_text(family = font_family, face = "plain", colour = "black",
                                       size = font_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = .9,
                                       margin = margin(), debug = FALSE),
      
      axis.line         = element_line(colour = "black", size = line_size, lineend = "square"),
      axis.line.x       = NULL,
      axis.line.y       = NULL,
      axis.text         = element_text(colour = "black", size = small_size),
      axis.text.x       = element_text(margin = margin(t = small_size / 4), vjust = 1),
      axis.text.x.top   = element_text(margin = margin(b = small_size / 4), vjust = 0),
      axis.text.y       = element_text(margin = margin(r = small_size / 4), hjust = 1),
      axis.text.y.right = element_text(margin = margin(l = small_size / 4), hjust = 0),
      axis.ticks        = element_line(colour = "black", size = line_size),
      axis.ticks.length = unit(half_line / 2, "pt"),
      axis.title.x      = element_text(
        margin = margin(t = half_line / 2),
        vjust = 1
      ),
      axis.title.x.top  = element_text(
        margin = margin(b = half_line / 2),
        vjust = 0
      ),
      axis.title.y      = element_text(
        angle = 90,
        margin = margin(r = half_line / 2),
        vjust = 1
      ),
      axis.title.y.right = element_text(
        angle = -90,
        margin = margin(l = half_line / 2),
        vjust = 0
      ),
      legend.background = element_blank(),
      legend.spacing    = unit(font_size, "pt"),
      legend.spacing.x  = NULL,
      legend.spacing.y  = NULL,
      legend.margin     = margin(0, 0, 0, 0),
      legend.key        = element_blank(),
      legend.key.size   = unit(1.1 * font_size, "pt"),
      legend.key.height = NULL,
      legend.key.width  = NULL,
      legend.text       = element_text(size = rel(rel_small)),
      legend.text.align = NULL,
      legend.title      = element_text(hjust = 0),
      legend.title.align = NULL,
      legend.position   = "right",
      legend.direction  = NULL,
      legend.justification = c("left", "center"),
      legend.box        = NULL,
      legend.box.margin =  margin(0, 0, 0, 0),
      legend.box.background = element_blank(),
      legend.box.spacing = unit(font_size, "pt"),
      
      panel.background  = element_blank(),
      panel.border      = element_blank(),
      panel.grid        = element_blank(),
      panel.grid.major  = NULL,
      panel.grid.minor  = NULL,
      panel.grid.major.x = NULL,
      panel.grid.major.y = NULL,
      panel.grid.minor.x = NULL,
      panel.grid.minor.y = NULL,
      panel.spacing     = unit(half_line, "pt"),
      panel.spacing.x   = NULL,
      panel.spacing.y   = NULL,
      panel.ontop       = FALSE,
      
      strip.background  = element_rect(fill = "grey80"),
      strip.text        = element_text(
        size = rel(rel_small),
        margin = margin(half_line / 2, half_line / 2,
                        half_line / 2, half_line / 2)
      ),
      strip.text.x      = NULL,
      strip.text.y      = element_text(angle = -90),
      strip.placement   = "inside",
      strip.placement.x =  NULL,
      strip.placement.y =  NULL,
      strip.switch.pad.grid = unit(half_line / 2, "pt"),
      strip.switch.pad.wrap = unit(half_line / 2, "pt"),
      
      plot.background   = element_blank(),
      plot.title        = element_text(
        face = "bold",
        size = rel(rel_large),
        hjust = 0.5, vjust = 1,
        margin = margin(b = half_line)
      ),
      plot.subtitle     = element_text(
        size = rel(rel_small),
        hjust = 0, vjust = 1,
        margin = margin(b = half_line)
      ),
      plot.caption      = element_text(
        size = rel(rel_tiny),
        hjust = 1, vjust = 1,
        margin = margin(t = half_line)
      ),
      plot.tag           = element_text(
        face = "bold",
        hjust = 0, vjust = 0.7
      ),
      plot.tag.position = c(0, 1),
      plot.margin       = margin(half_line, half_line, half_line, half_line),
      
      complete = TRUE
    )
}

wes_palette = c(
  "#FD6467","#0B775E","#2F343B",
  "#87B2C7", "#9B5094", "#F1BB7B", "#3E442B",
  "#F98400","#09ACEC","#FFDE05", "#DF9D5A", "#F2D3BC", "#DAD6D6"
)

# https://stackoverflow.com/questions/35717353/split-violin-plot-with-ggplot2
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
  data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
  grp <- data[1, "group"]
  newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
  newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
  newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])

  if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
    stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
      1))
    quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
    aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
    aesthetics$alpha <- rep(1, nrow(quantiles))
    both <- cbind(quantiles, aesthetics)
    quantile_grob <- GeomPath$draw_panel(both, ...)
    ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
  }
  else {
    ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
  }
})

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}