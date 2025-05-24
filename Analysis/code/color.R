
library(viridis)
library(ggsci)
library(RColorBrewer)

# label 2 color
label2Col <- function(labels) {
  if (!is.factor(labels)) labels <- factor(labels)
  num <- length(levels(labels))
  colors <-
    if (num <= 9) {
      RColorBrewer::brewer.pal(num, "Set1")
    } else {
      gg_color_hue(num)
    }
  colors[labels]
}

# more colors (ggplot)
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  grDevices::hcl(h=hues, l=65, c=100)[1:n]
}

# Generate colors from a customed color palette
scPalette <- function(n) {
  colorSpace <- c('#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC','#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A',
                  '#E3BE00','#FB9A99','#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B','#00441B','#DEDC00','#B3DE69','#8DD3C7','#999999')
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  } else {
    colors <- grDevices::colorRampPalette(colorSpace)(n)
  }
  return(colors)
}

# three colors for heatmap
c("#476fa9", "#ffffff", "#ca3226")
# or
c("#3658d3", "#ffffff", "#900504")

# negative and positive
colors_2 = c("#476fa9", "#ca3226")

# Un, neg, pos
colors_3 = c("#d6d6d6", "#476fa9", "#ca3226")


# 2 color for 2 cats (e.g. postive, negtive)
colors_2_d = c("darkred", "navyblue")
# alpha .3
c('#8B00004D', '#0000804D')


# specify color for cluster names
# 74 colors, from: https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
# then manually remove similar ones
# 54 colors
cluster_color_maps <- stats::setNames(c("#7FC97F", "#BEAED4", "#FFD92F", "#FB8072", "#386CB0", "#F0027F",
                                       "#BF5B17", "#666666", "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
                                       "#66A61E", "#E6AB02", "#E31A1C", "#A6CEE3", "#1F78B4", "#B2DF8A",
                                       "#33A02C", "#FB9A99", "#FDBF6F", "#A6761D", "#FF7F00", "#6A3D9A",
                                       "#B15928", "#FBB4AE", "#B3CDE3", "#CCEBC5", "#FED9A6", "#E5D8BD",
                                       "#F2F2F2", "#FDCDAC", "#FFED6F", "#F1E2CC", "#CCCCCC", "#377EB8",
                                       "#984EA3", "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62",
                                       "#8DA0CB", "#E78AC3", "#A6D854", "#E5C494", "#B3B3B3", "#8DD3C7",
                                       "#FFFF99", "#80B1D3", "#FCCDE5", "#D9D9D9", "#BC80BD", "#E6F5C9",
                                       "#FDC086"),
                                     as.character(0:54)
                                     )

# 7 colors
cell_color_maps <- c(
  'T cells' = '#197EC0FF',
  'B cells' = '#709AE1FF',
  'Myeloid' = '#8A9197FF',
  'Epithelia' = '#D2AF81FF',
  'Fibroblasts' = '#FD7446FF',
  'Endothelia' = '#71D0F5FF',
  'Platelets' = '#D5E4A2FF'
)

# 5 colors
tissue_color_maps <- c(
  "PBMC" = "#0073C2FF",
  "LN" = "#EFC000FF",
  "Normal" = "#7AA6DCFF",
  "Adjacent" = "#868686FF",
  "Tumor" = "#CD534CFF"
)

# 7 colors
origin2_color_maps <- c(
  'prePBMC' = '#023fa5', 'postPBMC' = '#7d87b9',
  'nLN' = '#ff7f03', 'pLN' = '#EFC000FF',
  'Normal' = '#7AA6DCFF', 'Adjacent' = '#868686FF',
  "Tumor" = "#CD534CFF"
)
# 7 colors
origin3_color_maps <- c(
  'PBMC1' = '#023fa5', 'PBMC2' = '#7d87b9',
  'LN_Down' = '#ff7f03', 'LN_Up' = '#EFC000FF',
  'Normal' = '#7AA6DCFF', 'Adjacent' = '#868686FF',
  "Tumor" = "#CD534CFF"
)

# 7 colors
origin4_color_maps <- c(
  'prePBMC' = '#023fa5', 'postPBMC' = '#7d87b9',
  'nLN' = '#bb7784', 'm_nLN' = '#ff7f03', 'm_pLN' = '#EFC000FF',
  'Normal' = '#7AA6DCFF', 'Adjacent' = '#868686FF',
  "Tumor" = "#CD534CFF"
)

# 6 colors
patient_color_maps <- c(
  "S0619" = "#4DBBD5FF",
  "S0730" = "#00A087FF",
  "S0819" = "#3C5488FF",
  "S0920" = "#F39B7FFF",
  "S1125" = "#8491B4FF",
  "S1204" = "#91D1C2FF"
)

# 17 colors
origin_color_maps <- c(
  'PBMC1' = '#023fa5', 'PBMC2' = '#7d87b9',
  'LN-LR' = '#336600', 'LN-RR' = '#d6bcc0', 'LN-UTP1' = '#bb7784',
  'LN-UTP2' = '#8e063b', 'LN-ITP' = '#4a6fe3', 'LN-LTP' = '#e07b91',
  'LN-LP' = '#d33f6a', 'LN-LP1' = '#11c638', 'LN-LP2' = '#8dd593',
  'LN-LG' = '#c6dec7',
  'Normal' = '#1CE6FF', 'Adjacent' = '#f0b98d',
  'Tumor_core' = '#ef9708', 'Tumor_invasion' = '#0fcfc0', 'Tumor_middle' = '#7f7f7f'
)




# 4 colors
vdj_color_maps <- c(
  'TCR' = '#36648b',
  'BCR' = '#7ccd7c',
  'Both' = '#ffd700',
  'None' = '#d6d6d6'
)

# 5 colors
ig_color_maps <- c(
  'IgM' = '#374E55FF',
  'IgD' = '#DF8F44FF',
  'IgG' = '#00A1D5FF',
  'IgA' = '#B24745FF',
  'None' = '#d6d6d6'
)

# IGHC colors
ighc_color_maps <- c(
  'IGHM' = '#17becf',
  'IGHD' = '#b5bd61',
  'IGHA1' = '#7f7f7f',
  'IGHA2' = '#e377c2',
  'IGHG1' = '#8c564b',
  'IGHG2' = '#aa40fc',
  'IGHG3' = '#d62728',
  'IGHG4' = '#279e68',
  'None' = '#d6d6d6'
)

# 3 colors
metastasis_color_maps <- c(
  'None' = '#d6d6d6',
  'N' = '#4DBBD5FF',
  'P' = '#E64B35FF'
)

# 3 colors
drainage_color_maps <- c(
  'None' = '#d6d6d6',
  'Down' = '#4DBBD5FF',
  'Up' = '#E64B35FF'
)

# 3 colors
epi_malig_colors <- c(
  "NonMalig." = "#7AA6DCFF",
  "PreMalig." = "#374E55FF",
  "Malig." = "#E64B35FF"
)


#####################################################
# T cells
#####################################################

t_level_1_color <- c(
  "CD4" = "#374E55FF",
  "CD8" = "#00A1D5FF",
  "Treg" = "#DF8F44FF",
  "γδT" = "#279e68",
  "NK/NKT" = "#B24745FF",
  "Unknown" = "#d6d6d6"
)

t_level_2_color = c(
  "Tn" = "#aec7e8",
  "Tcm" = "#1f77b4",
  "Tem" = "#ff7f0e",
  "Trm" = "#aa40fc",
  "Tfh" = "#e377c2",
  "Teff" = "#8c564b",
  "Tex" = "#d62728",
  "Treg" = "#DF8F44FF",
  "γδT" = "#279e68",
  "NK/NKT" = "#B24745FF",
  "Unknown" = "#d6d6d6"
)

t_level_3_color = c(
  # CD4
  "CD4-C1-Tcm" = "#1f77b4",
  "CD4-C2-Tn" = "#ff7f0e",
  "CD4-C3-Tn" = "#279e68",
  "CD4-C4-Tn" = "#d62728",
  "CD4-C5-Tcm" = "#aa40fc",
  "CD4-C6-Tn" = "#8c564b",
  "CD4-C7-Tn" = "#e377c2",
  "CD4-C8-Tfh" = "#b5bd61",
  "CD4-C9-Tn" = "#17becf",
  "CD4-C10-Tn" = "#aec7e8",
  "CD4-C11-Tem" = "#ffbb78",
  "CD4-C12-Tcm" = "#98df8a",
  # CD8
  "CD8-C1-Tcm" = "#1f77b4",
  "CD8-C2-Tex" = "#ff7f0e",
  "CD8-C3-Tem" = "#279e68",
  "CD8-C4-Teff" = "#d62728",
  "CD8-C5-Teff" = "#aa40fc",
  "CD8-C6-Teff" = "#8c564b",
  "CD8-C7-Tex" = "#e377c2",
  "CD8-C8-Tn" = "#b5bd61",
  "CD8-C9-Trm" = "#17becf",
  "CD8-C10-Teff" = "#aec7e8",
  "CD8-C11-Tex" = "#ffbb78",
  "CD8-C12-Tex" = "#98df8a",
  # Treg
  "Treg-C1" = "#1f77b4",
  "Treg-C2" = "#ff7f0e",
  "Treg-C3" = "#279e68",
  # γδT
  "γδT-C1" = "#1f77b4",
  "γδT-C2" = "#ff7f0e",
  "NK/NKT" = "#B24745FF",
  "Unknown" = "#d6d6d6"
)

#####################################################
# Myeloid
#####################################################

mye_level_1_color <- c(
  "Mono" = "#00A1D5FF",
  "MΦ" = "#374E55FF",
  "DC" = "#DF8F44FF",
  "Mast" = "#279e68",
  "CD34+ cells" = "#B24745FF"
)

# c("cMono", "ncMono", "MΦ", "cDC1", "cDC2", "tDC", "pDC", "Mast", "CD34+ cells")
mye_level_2_color = c(
  "cMono" = "#1f77b4",
  "ncMono" = "#aec7e8",
  "MΦ" = "#ff7f0e",
  "cDC1" = "#aa40fc",
  "cDC2" = "#e377c2",
  "mregDC" = "#8c564b",
  "pDC" = "#d62728",
  "Mast" = "#DF8F44FF",
  "γδT" = "#279e68",
  "CD34+ cells" = "#B24745FF"
)

# c("cMono-C1", "cMono-C2", "cMono-C3", "cMono-C4", "ncMono-C5", "cMono-C6",
# "MΦ-C1", "MΦ-C2", "MΦ-C3", "MΦ-C4", "MΦ-C5", "MΦ-C6", "MΦ-C7",
# "cDC2-C1", "cDC2-C2", "tDC-C3", "pDC-C4", "cDC1-C5", "cDC1-C6", "pDC-C7", "pDC-C8",
#"Mast-C1", "Mast-C2", "Mast-C3",
# "CD34+ cells")
mye_level_3_color = c(
  # Mono
  "cMono-C1" = "#1f77b4",
  "cMono-C2" = "#ff7f0e",
  "cMono-C3" = "#279e68",
  "cMono-C4" = "#d62728",
  "ncMono-C5" = "#aa40fc",
  "cMono-C6" = "#8c564b",
  "MΦ-C1" = "#e377c2",
  "MΦ-C2" = "#b5bd61",
  "MΦ-C3" = "#17becf",
  "MΦ-C4" = "#aec7e8",
  "MΦ-C5" = "#ffbb78",
  "MΦ-C6" = "#98df8a",
  "MΦ-C7" = "#ff9896",
  # DC
  "cDC2-C1" = "#1f77b4",
  "cDC2-C2" = "#ff7f0e",
  "mregDC-C3" = "#279e68",
  "pDC-C4" = "#d62728",
  "cDC1-C5" = "#aa40fc",
  "cDC1-C6" = "#8c564b",
  "pDC-C7" = "#e377c2",
  "pDC-C8" = "#b5bd61",
  # Mast
  "Mast-C1" = "#1f77b4",
  "Mast-C2" = "#ff7f0e",
  "Mast-C3" = "#279e68",
  # CD34+
  "CD34+ cells" = "#B24745FF"
)

#####################################################
# endothelia
#####################################################

endo_level_1_color <- c(
  "blood" = "#374E55FF",
  "lymphatic" = "#00A1D5FF"
)

# in case I will need
endo_level_2_color_old = c(
  "arteries" = "#aec7e8",
  "pcv-C1" = "#1f77b4",
  "pcv-C2" = "#ff7f0e",
  "pcv-C3" = "#aa40fc",
  "pcv-C4" = "#e377c2",
  "pcv-C5" = "#8c564b",
  "immature" = "#d62728",
  "tip cell" = "#DF8F44FF",
  "normal lymphatics" = "#279e68",
  "tumor lymphatics" = "#B24745FF"
)

endo_level_2_color = c(
  "arteries" = "#aec7e8",
  "pcv-C1" = "#1f77b4",
  "pcv-C2" = "#ff7f0e",
  "pcv-C3" = "#aa40fc",
  "pcv-C4" = "#e377c2",
  "pcv-C5" = "#8c564b",
  "immature" = "#d62728",
  "tip cell" = "#DF8F44FF",
  "LEC-C1" = "#279e68",
  "LEC-C2" = "#B24745FF"
)

endo_level_1.5_color = c(
  "arteries" = "#aec7e8",
  "pcv" = "#1f77b4",
  "immature" = "#d62728",
  "tip cell" = "#DF8F44FF",
  "LEC" = "#B24745FF"
)

#####################################################
# fibroblasts
#####################################################

fib_level_1_color <- c(
  "fibroblasts" = "#374E55FF",
  "pericyte" = "#00A1D5FF",
  "VSMC" = "#DF8F44FF"
)


fib_level_2_color = c(
  "NMF" = "#aec7e8",
  "myCAF" = "#1f77b4",
  "iCAF" = "#ff7f0e",
  "apCAF" = "#aa40fc",
  "pericyte" = "#279e68",
  "VSMC" = "#B24745FF"
)

fib_level_3_color = c(
  "NMF-C1" = "#aec7e8",
  "NMF-C2" = "#e377c2",
  "NMF-C3" = "#b5bd61",
  "NMF-C4" = "#17becf",
  "myCAF-C1" = "#1f77b4",
  "myCAF-C2" = "#ffbb78",
  "myCAF-C3" = "#98df8a",
  "iCAF-C1" = "#ff7f0e",
  "iCAF-C2" = "#8c564b",
  "iCAF-C3" = "#ff9896",
  "pericyte-C1" = "#279e68",
  "pericyte-C2" = "#d62728",
  "pericyte-C3" = "#8c6d31",
  "apCAF" = "#aa40fc",
  "VSMC" = "#B24745FF"
)

#####################################################
# B cells
#####################################################

b_level_1_color = c(
  "Naive" = "#aec7e8",
  "Mem" = "#1f77b4",
  "PC" = "#ff7f0e",
  "GCB" = "#aa40fc",
  "DN" = "#279e68",
  "CD5" = "#B24745FF"
)

b_level_2_color = c(
  "Naive-C1" = "#1f77b4",
  "Naive-C2" = "#ff7f0e",
  "Naive-C3" = "#279e68",
  "Naive-C4" = "#d62728",
  "Naive-C5" = "#aa40fc",
  "Naive-C6" = "#8c564b",
  "Naive-C7" = "#e377c2",
  "Naive-C8" = "#b5bd61",
  "Mem-C1" = "#17becf",
  "Mem-C2" = "#aec7e8",
  "Mem-C3" = "#ffbb78",
  "Mem-C4" = "#98df8a",
  "Mem-C5" = "#1f77b4",
  "Mem-C6" = "#ff7f0e",
  "Mem-C7" = "#279e68",
  "Mem-C8" = "#d62728",
  "PC-C1" = "#aa40fc",
  "PC-C2" = "#8c564b",
  "PC-C3" = "#e377c2",
  "PC-C4" = "#b5bd61",
  "GCB-C1" = "#17becf",
  "GCB-C2" = "#aec7e8",
  "GCB-C3" = "#ffbb78",
  "DN-C1" = "#98df8a",
  "DN-C2" = "#1f77b4",
  "CD5-C1" = "#ff7f0e",
  "CD5-C2" = "#279e68"
)

# un-used colors

# from scanpy
# https://github.com/theislab/scanpy/blob/master/scanpy/plotting/palettes.py
vega_10_scanpy = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc',
                   '#8c564b', '#e377c2', '#7f7f7f', '#b5bd61', '#17becf')

vega_20_scanpy = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc',
                   '#8c564b', '#e377c2', '#b5bd61', '#17becf', '#aec7e8',
                   '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94',
                   '#f7b6d2', '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31')

zeileis_28 = c("#023fa5", "#7d87b9", "#bec1d4", "#d6bcc0", "#bb7784",
               "#8e063b", "#4a6fe3", "#8595e1", "#b5bbe3", "#e6afb9",
               "#e07b91", "#d33f6a", "#11c638", "#8dd593", "#c6dec7",
               "#ead3c6", "#f0b98d", "#ef9708", "#0fcfc0", "#9cded6",
               "#d5eae7", "#f3e1eb", "#f6c4e1", "#f79cd4", '#7f7f7f',
               "#c7c7c7", "#1CE6FF", "#336600")
godsnot_102 = c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941",
                "#006FA6", "#A30059", "#FFDBE5", "#7A4900", "#0000A6",
                "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601",
                "#3B5DFF", "#4A3B53", "#FF2F80", "#61615A", "#BA0900",
                "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA",
                "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
                "#300018", "#0AA6D8", "#013349", "#00846F", "#372101",
                "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2",
                "#C2FF99", "#001E09", "#00489C", "#6F0062", "#0CBD66",
                "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459",
                "#456648", "#0086ED", "#886F4C", "#34362D", "#B4A8BD",
                "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
                "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF",
                "#3B9700", "#04F757", "#C8A1A1", "#1E6E00", "#7900D7",
                "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600",
                "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625",
                "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98",
                "#A4E804", "#324E72")

