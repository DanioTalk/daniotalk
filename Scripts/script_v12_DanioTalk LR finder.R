
######TO EDIT###################
Input_filename = "Data Sheet.xlsx"
FC_MIN_cutoff=1.5
P_value_cutoff=0.05
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
  filter(`P-value` < P_value_cutoff) %>%
  mutate(across(where(is.character), tolower)) %>%
  mutate(across(
    where(is.character),
    ~ stringr::str_replace(.x, "\\(1 of many\\)", "")
  )) %>%
  mutate(Gene = purrr::map_chr(Gene, filterSplicing)) %>%
  mutate(across(where(is.character), ~ stringr::str_trim(.x) %>% stringr::str_to_title()))

print("Updating gene names")
data$Gene = data$Gene %>% purrr::map_chr(check_in_map) %>% tolower()


print("Adding human orthologue gene information...")
human = read_tsv("human_orthos.txt")
colnames(human) = c("ZFIN_ID",
                    "ZFIN_SYMBOL",
                    "ZFIN_NAME",
                    "HUMAN_SYMBOL")

print("Downloading drug data...")
drugs = data.table::fread("drug.target.interaction.tsv") %>% 
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

print("Reading GO db...")
gos = readxl::read_excel("GO ID_db.xlsx", col_names = c("Gene", "GO Term", "Type")) %>%
  mutate(across(where(is.character), tolower))
print("Updating names in GO db...")
gos$Gene = purrr::map_chr(gos$Gene, check_in_map)
gos$Type = purrr::map_chr(gos$Type, toupper)

gos = gos %>% group_by(Gene) %>% mutate(Type = paste0(unique(Type), collapse = " & "))
gos_map = c(gos$Type) %>% setNames(gos$Gene)

print("Selecting columns, computing information about plasma ligands...")
final = joined %>%
  select(
    matches("Cell"),
    matches("Gene\\."),
    matches("FC\\."),
    matches("P-value\\."),
    matches("IID"),
    matches("DRUG_NAME\\."),
    matches("Human"),
    PHYSICAL_INTERACTION_SCORE
  )

final = final %>% mutate(
  Lig.type = purrr::map_chr(Gene.lig, ~ gos_map[.x]),
  Rec.type = purrr::map_chr(Gene.rec, ~ gos_map[.x])
)


IID_mapping = c(
  "pred" = 250,
  "exp" = 500,
  "exp|pred" = 500,
  "ortho|pred" = 250,
  "exp|ortho" = 500,
  "ortho" = 250,
  "exp|ortho|pred" = 500
)

# IID to numeric IID
final = final %>% mutate(`IID_numeric` = ifelse(is.na(IID_EVIDENCE), 0, IID_mapping[IID_EVIDENCE]))

## ADDED RANKS 10-11-2022
# Ranking 1
final = final %>% mutate(
  `ScoreForRank1` = ifelse(
    is.na(PHYSICAL_INTERACTION_SCORE) |
      PHYSICAL_INTERACTION_SCORE == 0,
    IID_numeric,
    PHYSICAL_INTERACTION_SCORE
  ),
  `ScoreForRank2` = FC.lig * FC.rec * ScoreForRank1
) 

final = final %>% mutate(
  Ranking1 = dense_rank(desc(ScoreForRank1)),
  Ranking2 = dense_rank(desc(ScoreForRank2))
) %>% select(-c(IID_numeric))

final_unfiltered = final %>% 
  mutate(across(where(is.character), ~ stringr::str_to_title(.x)))
final_filtered = final_unfiltered %>% filter(FC.lig > FC_MIN_cutoff, FC.rec > FC_MIN_cutoff)

final_name = "Resulting-ligand-receptor-pairs.xlsx"
final_name_unfiltered = "Resulting-ligand-receptor-pairs-unfiltered.xlsx"
print(paste("Saved interactome to file:", final_name))
print(paste("Saved unfiltered interactome to file:", final_name_unfiltered))
writexl::write_xlsx(final_filtered, final_name)
writexl::write_xlsx(final_unfiltered, final_name_unfiltered)

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
  ) %>% filter(`Gene Type` != "Not in db")


print("Adding information about Ligand/Receptor type...")
final = genes_with_cell %>%
  mutate(
    `Ligand Type` = ifelse(Gene %in% names(gos_map) &
                             `Gene Type` == "Ligand", gos_map[Gene], NA),
    `Receptor Type` = ifelse(Gene %in% names(gos_map) &
                               `Gene Type` == "Receptor", gos_map[Gene], NA)
  )

final_name_genes = "Resulting-ligand-receptor-list-plasma-and-GO.xlsx"
print(paste("Saved genes list to file:", final_name_genes))
writexl::write_xlsx(final, final_name_genes)

cell_counts_unfiltered = "Resulting-cell-type-counts-unfiltered.xlsx"
cell_counts = "Resulting-cell-type-counts.xlsx"
df_unfiltered = final %>% group_by(`Cell type`) %>% summarize(
  `Ligand count` = sum(`Gene Type` == "Ligand"),
  `Receptor count` = sum(`Gene Type` == "Receptor")
)
df = final %>% filter(FC > FC_MIN_cutoff) %>% 
  group_by(`Cell type`) %>% 
  summarize(
  `Ligand count` = sum(`Gene Type` == "Ligand"),
  `Receptor count` = sum(`Gene Type` == "Receptor")
)
writexl::write_xlsx(df, cell_counts)
writexl::write_xlsx(df_unfiltered, cell_counts_unfiltered)

ligands_and_receptors = "Resulting-ligands-receptors-each-cell.xlsx"
ligands_and_receptors_unfiltered = "Resulting-ligands-receptors-each-cell-unfiltered.xlsx"
df_list = final %>%
  filter(FC > FC_MIN_cutoff) %>%
  select(`Cell type`, Gene, `Gene Type`, FC, `P-value`,`DRUG_NAME`) %>%
  split(.$`Cell type`) %>%
  purrr::map(~ .x %>% select(-`Cell type`))
df_list_unfiltered = final %>%
  select(`Cell type`, Gene, `Gene Type`, FC, `P-value`, `DRUG_NAME`) %>%
  split(.$`Cell type`) %>%
  purrr::map(~ .x %>% select(-`Cell type`))

openxlsx::write.xlsx(df_list, ligands_and_receptors)
openxlsx::write.xlsx(df_list_unfiltered, ligands_and_receptors_unfiltered)
print(glue(
  "Saved count and list of genes to {cell_counts} and {ligands_and_receptors}"
))
