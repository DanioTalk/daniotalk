if(!require(tidyverse)){
  install.packages("tidyverse")
  require(tidyverse)
}
if(!require(circlize)){
  remotes::install_github("jokergoo/circlize")
  library(circlize)
}
if(!require(glue)){
  install.packages("glue")
  require(glue)
}

df = readxl::read_excel("Resulting-ligand-receptor-pairs.xlsx")

#########################
############# SETTINGS ##
#########################
# For cells just type "any_cell" to use all cell types
cell_lig = "any_cell"
cell_rec = "any_cell"

# Ranking1, Ranking2 (Rank1: STRING or IID, Rank2: Rank1 + Fold-Change data)
ranking = "Ranking1"

# How many interactions to plot?
howmany = 50

# Use below lists with quotations
# Like that:
# receptors_list = c("rec1", "rec2") etc.
#
# Leave empty to use all receptors
receptors_list = c()

# Leave empty to use all ligands
ligands_list = c()
#ligands_list = c("ndr1")

# Confidence level
# This treshold works only for Ranking1 type because
# Scores produced for Ranking2 are a in other range
# (almost always > 900 because of the FC scaling factor)
# Highest = Score >= 900
# High Score >= 700
# Medium = Score >= 400
# Any = any score
confidence_level = "Highest"

# Plot settings
transparency_level = 0.2
color_palette = "Inferno" # Run command 'hc.pals()' in console to display all palettes

# Output file name
outname = glue("Receptors-{cell_rec}-vs-Ligands-{cell_lig}-top-{howmany}-interactions-ranking-{ranking}-{confidence_level}-confidence.tiff")

##### END OF SETTINGS



final = df %>% select(Gene.rec,
                      Gene.lig,
                      `Cell type.rec`,
                      `Cell type.lig`,
                      Ranking1, 
                      Ranking2,
                      ScoreForRank1, 
                      ScoreForRank2)  %>%
  mutate(V1 = "Receptor", V2 = "Ligand", Rank = case_when(
    ranking == "Ranking1" ~ Ranking1,
    ranking == "Ranking2" ~ Ranking2
  ),
  Score = case_when(
    ranking == "Ranking1" ~ ScoreForRank1,
    ranking == "Ranking2" ~ ScoreForRank2
  )) %>%
  select(Gene.rec, Gene.lig, `Cell type.rec`, `Cell type.lig`, Score, Rank)

# Cell type filtering
if(cell_lig != "any_cell"){
  final = final %>% filter(
    `Cell type.lig` == cell_lig
  )
}
if(cell_rec != "any_cell"){
  final = final %>% filter(
    `Cell type.rec` == cell_rec
  )
}

# Ligands and receptors list filtering
if(length(receptors_list) != 0){
  final = final %>% filter(Gene.rec %in% receptors_list)
}

if(length(ligands_list) != 0){
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

if(nrow(final) == 0) {
  stop("You have an empty dataset after filtering. Check if the confidence level is not too high.")
}

f = final %>% arrange(Rank) %>% 
  head(howmany) %>% 
  select(Gene.rec, Gene.lig, `Cell type.lig`, `Cell type.rec`, Score)

f = f %>% mutate(
  Gene.rec = str_to_title(Gene.rec),
  Gene.lig = str_to_title(Gene.lig)
)

f.mutated = f %>% mutate(
  Gene.rec = paste0(Gene.rec, "---", `Cell type.rec`),
  Gene.lig = paste0(Gene.lig, "---", `Cell type.lig`)
)



tiff(outname, 6000, 4000, res = 500, compression = "lzw")

circos.clear()
chordDiagram(
  f.mutated %>% select(Gene.rec, Gene.lig, Score),
  annotationTrack = "grid",
  grid.col = hcl.colors(c(f$Gene.rec, f$Gene.lig) %>% unique() %>% length(), color_palette),
  transparency = transparency_level,
  preAllocateTracks = 1,
)
circos.track(track.index = 1, panel.fun = function(x, y) {
  
  label = CELL_META$sector.index %>% str_split(., "---") %>% unlist() %>% .[1]
  
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], label, 
              facing = "clockwise", cex = 0.8, niceFacing = TRUE, adj = c(0, 1))
},  bg.border = NA)

dev.off()

