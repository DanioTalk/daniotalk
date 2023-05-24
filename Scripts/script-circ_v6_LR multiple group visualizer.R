library(tidyverse)
library(png)
if (!require(circlize)) {
  remotes::install_github("jokergoo/circlize")
  library(circlize)
}
library(glue)


df = readxl::read_excel("Resulting-ligand-receptor-pairs.xlsx")

#########################
############# SETTINGS ##
#########################
# For cells just type "any_cell" to use all cell types
cell_lig = c("any_cell")
cell_rec = c("any_cell")

# Ranking1, Ranking2, Ranking3 
# Rank1: STRING or Human ligand-receptor score
# Rank2: Rank1 * Fold-Change/Percentile, 
# Rank3: Human ligand-receptor score * Fold-Change/Percentile)
ranking = "Ranking1"

# How many interactions to plot?
howmany = 50

# Use below lists with quotations
# Like that:
# receptors_list = c("Rec1", "Rec2") etc.
#
# Leave empty to use all receptors
receptors_list = c()

# Leave empty to use all ligands
ligands_list = c()

# Confidence level
# Highest = Score >= 900
# High Score >= 700
# Medium = Score >= 400
# Any = any score
confidence_level = "High"

# Plot settings
transparency_level = .2 # from 0.0 to 1.0
# highlight_transparency = 255 # from 0 up to 255
color_palette = "Inferno" # Run command 'hcl.pals()' in console to display all palettes
cell_group_color_palette = "Viridis" # Color palette for cell grouping adnotation
cell_group_text_color = "white" # Run command 'colors()' to see all colors

##### END OF SETTINGS

outname = glue(
  "Receptors-{cell_rec}-vs-Ligands-{cell_lig}-top-{howmany}-interactions-ranking-{ranking}-confidence-{confidence_level}.tiff"
)


final = df %>% select(
  Gene.rec,
  Gene.lig,
  `Cell type.rec`,
  `Cell type.lig`,
  Ranking1,
  Ranking2,
  Ranking3,
  ScoreForRank1,
  ScoreForRank2,
  ScoreForRank3
)  %>%
  mutate(
    V1 = "Receptor",
    V2 = "Ligand",
    Rank = case_when(
      ranking == "Ranking1" ~ Ranking1,
      ranking == "Ranking2" ~ Ranking2,
      ranking == "Ranking3" ~ Ranking3
    ),
    Score = ScoreForRank1
  ) %>%
  select(Gene.rec, Gene.lig, `Cell type.rec`, `Cell type.lig`, Score, Rank)

# Cell type filtering
if (!("any_cell" %in% cell_lig)) {
  final = final %>% filter(tolower(`Cell type.lig`) %in% tolower(cell_lig))
}
if (!("any_cell" %in% cell_rec)) {
  final = final %>% filter(tolower(`Cell type.rec`) %in% tolower(cell_rec))
}

# Ligands and receptors list filtering
if (length(receptors_list) != 0) {
  final = final %>% filter(Gene.rec %in% receptors_list)
}

if (length(ligands_list) != 0) {
  final = final %>% filter(Gene.lig %in% ligands_list)
}

# Confidence filtering
final = if (confidence_level == "Highest") {
  final %>% filter(Score >= 900)
} else if (confidence_level == "High") {
  final %>% filter(Score >= 700)
} else if (confidence_level == "Medium") {
  final %>% filter(Score >= 400)
} else {
  final
}


f = final %>% arrange(Rank) %>%
  head(howmany)
f = f %>%
  select(Gene.rec, Gene.lig, `Cell type.lig`, `Cell type.rec`, Score)

# Unique name-cell pairs
# genes = tibble(
#   cell = c(f$`Cell type.rec`, f$`Cell type.lig`),
#   gene = c(f$Gene.rec, f$Gene.lig)
# ) %>% distinct()
# 
# group = genes$cell %>% setNames(genes$gene)

f = f %>% mutate(
  Gene.rec = str_to_title(Gene.rec),
  Gene.lig = str_to_title(Gene.lig),
  `Cell type.lig` = str_to_title(`Cell type.lig`), 
  `Cell type.rec` = str_to_title(`Cell type.rec`)
)

f.mutated = f %>% mutate(
  Gene.rec = paste0(Gene.rec, "---", `Cell type.rec`),
  Gene.lig = paste0(Gene.lig, "---", `Cell type.lig`)
)

group = c(f.mutated$`Cell type.rec`, f.mutated$`Cell type.lig`) %>% 
  setNames(c(f.mutated$Gene.rec, f.mutated$Gene.lig))

group = group[unique(names(group))]


#####
tiff(outname, 6000, 4000, res = 500, compression = "lzw")

genesToAssignColors = c(f.mutated$Gene.rec, f.mutated$Gene.lig) %>% unique()
chordColors = hcl.colors(genesToAssignColors %>% length(), color_palette)

circos.clear()
circos.par(start.degree = -45)
chordDiagram(
  f.mutated %>% arrange(`Cell type.lig`, `Cell type.rec`) %>% select(Gene.rec, Gene.lig, Score),
  annotationTrack = c("grid"),
  group = group,
  transparency = transparency_level,
  grid.col = chordColors,
  preAllocateTracks = list(list(
    track.height = mm_h(6),
    track.margin = c(mm_h(18), 0)
    
  ), list(
    track.height = max(strwidth(names(group) %>% strsplit(., "---") %>% unlist() %>% .[1]))
  ))
)

cells = tibble(gene = names(group), cell = unname(group))
cells = split(cells, cells$cell)

cells_cols = hcl.colors(length(cells), cell_group_color_palette, fixup = T)

col_idx = 1
for(df in cells){
  highlight.sector(df$gene, track.index = 1, col = cells_cols[col_idx], padding = c(.001, .001, .001, .001),
                   text = df$cell %>% unique(), cex = 0.8, text.col = cell_group_text_color, niceFacing = TRUE)
  col_idx = col_idx + 1
}

lig_char = "\U25B2"
rec_char = "\U25A0"

circos.track(track.index = 2, panel.fun = function(x, y) {
  
  label = CELL_META$sector.index %>% str_split(., "---") %>% unlist() %>% .[1]
  label = paste(label, ifelse(label %in% f$Gene.rec, rec_char, lig_char))
  
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], label, 
              facing = "clockwise", cex = 1, niceFacing = TRUE, adj = c(0, 1))
  },  bg.border = NA)

dev.off()
