
######TO EDIT###################
Input_filename = "Data Sheet_test.xlsx"
Expression_MIN_cutoff=0.05  #Minimum expression cutoff value. 
#########################


print("Loading packages...")
if (!require("tidyverse")) {
  install.packages("tidyverse")
  library("tidyverse")
}

if (!require("writexl")) {
  install.packages("writexl")
}

if (!require("textclean")) {
  install.packages("textclean")
  library("textclean")
}

if(!require("data.table")){
  install.packages("data.table")
  require("data.table")
}

if(!require("ontologyIndex")){
  install.packages("ontologyIndex")
}


print("Loading gene aliases...")
gene_map = read.delim("aliases.txt",
                      header = F,
                      col.names = c("Symbol", "Alias")) %>%
  tibble() %>%
  mutate(across(where(is.character), tolower))

GENE_SET = gene_map$Symbol %>% unique()
ALIAS_DICT = gene_map %>% filter(Alias != "")
ALIAS_DICT  = c(ALIAS_DICT$Symbol %>%
                  as.character()) %>%
  setNames(ALIAS_DICT$Alias %>%
             as.character())

filterSplicing = function(g) {
  if((g %in% GENE_SET) | (g %in% names(ALIAS_DICT))){
    return(g)
  }
  result = stringr::str_split(g, "\\.") %>% unlist() %>% .[[1]]
  return(result)
}

gene_map = gene_map %>% filter(!is.na(Alias))

check_in_map = function(x) {
  if (x %in% GENE_SET) {
    return(x)
  } else if (x %in% names(ALIAS_DICT)) {
    return(ALIAS_DICT[x])
  } else {
    return(x)
  }
}


print("Loading data...
       Removing (1 of many).X")

data = readxl::read_excel(Input_filename)  %>%
  mutate(across(where(is.character), tolower)) %>%
  mutate(across(
    where(is.character),
    ~ stringr::str_replace(.x, "\\(1 of many\\)", "")
  )) %>%
  mutate(Gene = purrr::map_chr(Gene, filterSplicing)) %>%
  mutate(across(where(is.character), ~ stringr::str_trim(.x) %>% tolower()))

# Add percentile
data = data %>% group_by(`Cell type`) %>% mutate(`QuasiPercentile` = rank(Expression) / length(Expression))

print("Updating gene names")
data$Gene = data$Gene %>% purrr::map_chr(check_in_map) %>% tolower()


print("Adding human orthologue gene information...")
human = read_tsv("human_orthos.txt")
colnames(human) = c("ZFIN_ID",
                    "ZFIN_SYMBOL",
                    "ZFIN_NAME",
                    "HUMAN_SYMBOL")

if(!file.exists("drug.target.interaction.tsv.gz")){
  print("Downloading drug data...")
  download.file("https://unmtid-shinyapps.net/download/DrugCentral/2021_09_01/drug.target.interaction.tsv.gz", "drug.target.interaction.tsv.gz")
}
drugs = data.table::fread("drug.target.interaction.tsv.gz") %>% 
  tibble() %>%
  mutate(GENE = strsplit(GENE, split = "|", fixed = T)) %>%
  unnest(GENE) %>%
  distinct() %>%
  select(GENE, DRUG_NAME) %>%
  mutate(GENE = tolower(GENE)) %>%
  group_by(GENE) %>% 
  summarize(DRUG_NAME = paste(DRUG_NAME, collapse = ",")) %>%
  ungroup()

data = data %>% left_join(drugs, by = c("Gene" = "GENE"))

human = human %>%
  select(ZFIN_SYMBOL, HUMAN_SYMBOL) %>%
  distinct() %>%
  mutate_all( ~ .x %>% textclean::strip(digit.remove = F, apostrophe.remove = F))
data = data %>%
  left_join(human, by = c("Gene" = "ZFIN_SYMBOL")) %>%
  plyr::rename(., c("HUMAN_SYMBOL" = "Human ortholog"))

genes = data$Gene %>% unique()

print("Fetching GO's from remote...")
valid_gos = "GO:0001664
GO:0005184
GO:0005179
GO:0042379
GO:0008009
GO:0005125
GO:0005126
GO:0031628
GO:0008083
GO:0004930
GO:0008188
GO:0042923
GO:0030594
GO:0015276
GO:0004896
GO:0001594
GO:0004985
GO:0004879" %>%
  stringr::str_split("\n") %>%
  unlist() %>%
  stringr::str_replace(" ", "")
if(!file.exists("go.obo")){
  download.file("http://purl.obolibrary.org/obo/go.obo", "go.obo")
}
bo = ontologyIndex::get_OBO("go.obo")
bo_map = bo$name %>% setNames(bo$id)


if(!file.exists("zfin.gaf.gz")){
  download.file("http://current.geneontology.org/annotations/zfin.gaf.gz", "zfin.gaf.gz")
}
gos_zfin = data.table::fread("zfin.gaf.gz", skip = "!") %>% 
  select(1:5) %>% 
  setNames(c("uniprot", "code", "Gene", "1", "GO")) %>%
  select(Gene, GO) %>%
  mutate(Gene = tolower(Gene)) %>%
  mutate(Gene = ifelse(Gene %in% names(gene_map), gene_map[Gene], Gene)) %>%
  mutate(`Type` = bo_map[GO]) %>% select(Gene, Type) %>%
  group_by(Gene) %>% summarise(Type = paste(unique(Type), collapse = ";")) %>% ungroup() # %>% filter(GO %in% valid_gos)

gos_zfin$Gene = purrr::map_chr(gos_zfin$Gene, check_in_map)
gos_zfin$Type = purrr::map_chr(gos_zfin$Type, toupper)

# gos = gos_zfin %>% group_by(Gene) %>% mutate(Type = paste0(unique(Type), collapse = " & "))
gos_map = c(gos_zfin$Type) %>% setNames(gos_zfin$Gene)

if(!file.exists("goa_human.gaf.gz")){
  download.file("http://geneontology.org/gene-associations/goa_human.gaf.gz", "goa_human.gaf.gz")
}
gos_human = data.table::fread("goa_human.gaf.gz", skip="!") %>% 
  select(1:5) %>% 
  setNames(c("uniprot", "code", "Gene", "1", "GO")) %>%
  select(Gene, GO) %>%
  mutate(Gene = tolower(Gene)) %>%
  mutate(Gene = ifelse(Gene %in% names(gene_map), gene_map[Gene], Gene)) %>%
  mutate(`Type` = bo_map[GO]) %>% select(Gene, Type) %>%
  group_by(Gene) %>% summarise(Type = paste(unique(Type), collapse = ";")) %>% ungroup() # %>% filter(GO %in% valid_gos)

gos_human$Type = purrr::map_chr(gos_human$`Type`, toupper)
# gos_human = gos_human %>% group_by(Gene) %>% mutate(Type = paste0(unique(Type), collapse = " & "))

data = data %>% 
  left_join(gos_zfin %>% select(Gene, Type), by = "Gene") %>%
  rename(`ZFIN GO Type` = `Type`) %>% 
  left_join(gos_human %>% select(Gene, `Type`), by = c("Human ortholog" = "Gene")) %>%
  rename(`Human GO Type` = `Type`) 




print("Loading database...")
database = readr::read_csv("Database.csv")  %>%
  mutate(across(where(is.character), tolower))


print("Joining data onto database. Computing pairs...")
joined = database %>%
  filter(ZFIN_SYMBOL_LIG %in% genes, ZFIN_SYMBOL_REC %in% genes) %>%
  inner_join(data,
             by = c("ZFIN_SYMBOL_LIG" = "Gene"),
             keep = T) %>%
  inner_join(
    data,
    by = c("ZFIN_SYMBOL_REC" = "Gene"),
    keep = T,
    suffix = c(".lig", ".rec")
  )

print("Loading plasma data...")
plasma_exp = readxl::read_excel("Plasma ligands_expt.xlsx", col_names = "Symbol") %>%
  mutate(across(where(is.character), tolower),
         Symbol = purrr::map_chr(Symbol, check_in_map))

plasma_pred = readxl::read_excel("Plasma ligands_predicted.xlsx", col_names = "Symbol") %>%
  mutate(across(where(is.character), tolower),
         Symbol = purrr::map_chr(Symbol, check_in_map))


print("Selecting columns, computing information about plasma ligands...")
final = joined %>%
  select(
    matches("Cell type"),
    matches("Gene\\."),
    matches("Expression\\."),
    matches("IID"),
    matches("DRUG_NAME\\."),
    matches("QuasiPercentile\\."),
    matches("Human"),
    matches("Algorithms"),
    matches("ConservationScore"),
    matches("GO Type"),
    matches("CellTalk"),
    matches("CellCell"),
    matches("Ramilowski"),
    matches("CellPhone"),
    matches("Ligand_type"),
    matches("CellChat"),
    matches("Human_Matrisome"),
    matches("ZFIN_Matrisome"),
    PHYSICAL_INTERACTION_SCORE
  )


final = final %>%
  relocate(ConservationScore_LIG, Algorithms_LIG, .after = HUMAN_SYMBOL_LIG) %>%
  relocate(ConservationScore_REC, Algorithms_REC, .after = HUMAN_SYMBOL_REC)


# IID to numeric IID
final = final %>%
  mutate(
    `Human LR pair score` = case_when(
      `CellTalkDB` == "yes" | `Ramilowski` == "yes" | `CellPhoneDB` == "yes" | `CellChat` == "yes" | !is.na(CellChat_Dimer_or_Complex) | !is.na(CellPhoneDB_Dimer_or_Complex) ~ 500,
      `IID_EVIDENCE` == "exp" | IID_EVIDENCE == "exp|pred"  | IID_EVIDENCE == "exp|ortho" | IID_EVIDENCE == "exp|ortho|pred" ~ 400,
      `CellCellInteractions` == "yes" | `IID_EVIDENCE` == "pred" | `IID_EVIDENCE` == "ortho|pred" | `IID_EVIDENCE` == "ortho" ~ 250,
      TRUE ~ 0
    )
  )

## ADDED RANKS 10-11-2022
# Ranking 1
final = final %>% mutate(
  `ScoreForRank1` = ifelse(
    is.na(PHYSICAL_INTERACTION_SCORE) |
      PHYSICAL_INTERACTION_SCORE < 200,
    `Human LR pair score`,
    PHYSICAL_INTERACTION_SCORE
  ),
  `ScoreForRank2` = QuasiPercentile.lig * QuasiPercentile.rec * ScoreForRank1,
  `ScoreForRank3` = QuasiPercentile.lig * QuasiPercentile.rec * `Human LR pair score`
) 

final = final %>% mutate(
  Ranking1 = dense_rank(desc(ScoreForRank1)),
  Ranking2 = dense_rank(desc(ScoreForRank2)),
  Ranking3 = dense_rank(desc(ScoreForRank3))
) %>% select(-c(`Human LR pair score`)) %>%
  select(
    `Cell type.lig`,
    `Cell type.rec`,
    `Gene.lig`,
    `Gene.rec`,
    `Expression.lig`,
    `Expression.rec`,
    `QuasiPercentile.lig`,
    `QuasiPercentile.rec`,
    `Ranking1`,
    `Ranking2`,
    `Ranking3`,
    `PHYSICAL_INTERACTION_SCORE`,
    `CellCellInteractions`,
    `CellChat`,
    `CellChat_Dimer_or_Complex`,
    `CellPhoneDB`,
    `CellPhoneDB_Dimer_or_Complex`,
    `CellTalkDB`,
    `IID_EVIDENCE`,
    `Ramilowski`,
    `DRUG_NAME.lig`,
    `HUMAN_SYMBOL_LIG`,
    `ConservationScore_LIG`,
    `Algorithms_LIG`,
    `HUMAN_GENE_ID_LIG`,
    `DRUG_NAME.rec`,
    `HUMAN_SYMBOL_REC`,
    `ConservationScore_REC`,
    `Algorithms_REC`,
    `HUMAN_GENE_ID_REC`,
    `ZFIN_Matrisome_annotation_LIG`,
    `ZFIN GO Type.lig`,
    `Human GO Type.lig`,
    `Human_Matrisome_annotation_LIG`,
    `Ligand_type`,
    `CellChat_annotation`,
    `ZFIN_Matrisome_annotation_REC`,
    `ZFIN GO Type.rec`,
    `Human GO Type.rec`,
    `Human_Matrisome_annotation_REC`,
    `ScoreForRank1`,
    `ScoreForRank2`,
    `ScoreForRank3`
  ) %>% mutate(
  HUMAN_SYMBOL_LIG = toupper(HUMAN_SYMBOL_LIG),
  HUMAN_SYMBOL_REC = toupper(HUMAN_SYMBOL_REC)
)


final_unfiltered = final
final_filtered = final_unfiltered %>% filter(Expression.lig > Expression_MIN_cutoff, Expression.rec > Expression_MIN_cutoff)

final_name = "Resulting-ligand-receptor-pairs.xlsx"
final_name_unfiltered = "Resulting-ligand-receptor-pairs-unfiltered.xlsx"

writexl::write_xlsx(final_filtered, final_name)
writexl::write_xlsx(final_unfiltered, final_name_unfiltered)

print(paste("Saved interactome to file:", final_name))
print(paste("Saved unfiltered interactome to file:", final_name_unfiltered))

tmp = database %>% select(matches("SYMBOL"), matches("Algorithms"), matches("Conservation"))
tmp1 = tmp %>% select(ZFIN_SYMBOL_LIG, Algorithms_LIG, ConservationScore_LIG)
tmp2 = tmp %>% select(ZFIN_SYMBOL_REC, Algorithms_REC, ConservationScore_REC)

colnames(tmp1) = c("Symbol", "Algorithms", "ConservationScore")
colnames(tmp2) = c("Symbol", "Algorithms", "ConservationScore")

genes_metadata = rbind(tmp1, tmp2)


genes_with_cell = data %>%
  mutate(
    `Gene Type` = case_when(
      Gene %in% database$ZFIN_SYMBOL_LIG ~ "Ligand",
      Gene %in% database$ZFIN_SYMBOL_REC ~ "Receptor",
      TRUE ~ "Not in db"
    ),
    `Plasma` = case_when(
      Gene %in% plasma_pred$Symbol &
        Gene %in% plasma_exp$Symbol ~ "Plasma Exp & Pred",
      Gene %in% plasma_pred$Symbol ~ "Plasma Pred",
      Gene %in% plasma_exp$Symbol ~ "Plasma Exp",
      TRUE ~ "None"
    ),
    `Plasma` = ifelse(`Gene Type` == "Ligand", `Plasma`, NA)
  ) %>% filter(`Gene Type` != "Not in db") %>%
  inner_join(genes_metadata, by = c("Gene" = "Symbol")) %>% unique()


print("Adding information about Ligand/Receptor type...")
final = genes_with_cell %>%
  mutate(
    `Ligand Type` = ifelse(Gene %in% names(gos_map) &
                             `Gene Type` == "Ligand", gos_map[Gene], NA),
    `Receptor Type` = ifelse(Gene %in% names(gos_map) &
                               `Gene Type` == "Receptor", gos_map[Gene], NA)
  ) %>% unique()

final_name_genes = "Resulting-ligand-receptor-list-plasma-and-GO.xlsx"
print(paste("Saved genes list to file:", final_name_genes))
writexl::write_xlsx(final %>%
                      mutate(`Human ortholog` = toupper(`Human ortholog`)) %>% 
                      select(-`Ligand Type`, -`Receptor Type`) %>%
                      rename(`ZFIN GO` = `ZFIN GO Type`,
                             `Human GO` = `Human GO Type`), final_name_genes)

cell_counts_unfiltered = "Resulting-cell-type-counts-unfiltered.xlsx"
cell_counts = "Resulting-cell-type-counts.xlsx"
df_unfiltered = final %>% group_by(`Cell type`) %>% summarize(
  `Ligand count` = length(unique(Gene[`Gene Type` == "Ligand"])),
  `Receptor count` = length(unique(Gene[`Gene Type` == "Receptor"]))
)
df = final %>% filter(Expression > Expression_MIN_cutoff) %>% 
  group_by(`Cell type`) %>% 
  summarize(
    `Ligand count` = length(unique(Gene[`Gene Type` == "Ligand"])),
    `Receptor count` = length(unique(Gene[`Gene Type` == "Receptor"]))
  )
writexl::write_xlsx(df, cell_counts)
writexl::write_xlsx(df_unfiltered, cell_counts_unfiltered)

ligands_and_receptors = "Resulting-ligands-receptors-each-cell.xlsx"
ligands_and_receptors_unfiltered = "Resulting-ligands-receptors-each-cell-unfiltered.xlsx"
df_list = final %>%
  filter(Expression > Expression_MIN_cutoff) %>%
  select(`Cell type`, Gene, `Gene Type`, Expression, `DRUG_NAME`,`Human ortholog`) %>%
  inner_join(genes_metadata, by = c("Gene" = "Symbol")) %>%
  unique() %>%
  mutate(`Human ortholog` = toupper(`Human ortholog`)) %>%
  split(.$`Cell type`) %>%
  purrr::map(~ .x %>% select(-`Cell type`))
df_list_unfiltered = final %>%
  select(`Cell type`, Gene, `Gene Type`, Expression, `DRUG_NAME`,`Human ortholog`) %>%
  inner_join(genes_metadata, by = c("Gene" = "Symbol")) %>%
  unique() %>%
  mutate(`Human ortholog` = toupper(`Human ortholog`)) %>%
  split(.$`Cell type`) %>%
  purrr::map(~ .x %>% select(-`Cell type`))

openxlsx::write.xlsx(df_list, ligands_and_receptors)
openxlsx::write.xlsx(df_list_unfiltered, ligands_and_receptors_unfiltered)
print(glue(
  "Saved count and list of genes to {cell_counts} and {ligands_and_receptors}"
))
