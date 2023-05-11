#!/usr/bin/env python3
import logging

import pandas as pd
import pyreadr
import rdata
import sqlalchemy as sa
from collections import defaultdict
from functools import reduce

logging.basicConfig(encoding='utf-8', level=logging.DEBUG, format="%(asctime)s [%(levelname)s 🐠]: %(message)s", datefmt="%I:%M:%S")

print("""
______            _     _____     _ _    
|  _  \          (_)   |_   _|   | | |   
| | | |__ _ _ __  _  ___ | | __ _| | | __
| | | / _` | '_ \| |/ _ \| |/ _` | | |/ /
| |/ | (_| | | | | | (_) | | (_| | |   < 
|___/ \__,_|_| |_|_|\___/\_/\__,_|_|_|\_\
""")

def read_gene_map(file: str = "Data/aliases.txt"):
    """This function reads data from ZFIN aliases.txt file
    
    It creates the dictionary of gene aliases where:
    key - gene name (alias)
    value - newest gene name

    It also creates a set of newest gene name for checking the presence in O(1) time.

    Returns: returns set of all newest names, and dictionary described above
    """
    # Dict and set:
    # set - genes
    # dict - alias: gene
    df = (
        pd.read_csv(file, sep="\t", names=["Gene", "Alias"])
        .apply(lambda x: x.astype(str).str.lower().str.strip())
        .drop_duplicates(ignore_index=True)
    )

    s = set() 
    d_r = dict()
    for x in (tuple(arr) for arr in df.values):
        s.add(x[0])
        if not pd.isna(x[1]):
            d_r[x[1]] = x[0]
    return s, d_r


def check_in_map(x: str, gene_set: set, alias_map: dict[str, str]):
    """This function takes a gene name and returns it's newest name.
    
    Params:
        x: gene name
        gene_set: set containing all the newest names for all genes from aliases.txt computed using `read_gene_map` func
        alias_map: dictionary computed using `read_gene_map` func
    Returns:
        (str): newest gene name
    """
    
    if x in gene_set:
        return x
    elif x in alias_map:
        return alias_map[x]
    else:
        return x


def replace_names(df: pd.DataFrame, column: str, gene_set: set, alias_map: dict[str, str]):
    """This function takes dataframe with specific column and 
    replaces all the gene names in this column with their newest names"""

    series = df[column].str.lower().str.strip()
    new_series = series.apply(check_in_map, gene_set=gene_set, alias_map=alias_map)
    df[column] = new_series
    return df

# Global variables, gene set with newest names and dictionary of aliases
GENE_SET, ALIAS_MAP = read_gene_map()


def get_zfin_data(zfin_file="Data/zfin_genes.txt", ensembl_file="Data/ensembl_assoc.txt", ncbi_file="Data/ncbi_assoc.txt", uniprot_file="Data/uni_assoc.txt"):
    """This function reads all the association files downloaded from ZFIN and
    merges them into one DataFrame
    
    Params:
        Files downloaded from ZFIN
    Returns:
        Single dataframe merged from parsed and cleaned ZFIN files
    """
    zfin_df = pd.read_csv(
        zfin_file,
        sep="\t",
        names=["ZFIN_ID", "SO_ID", "ZFIN_SYMBOL", "SEQ_ACC"],
    ).drop(columns=["SO_ID", "SEQ_ACC"])

    ensembl_df = pd.read_csv(
        ensembl_file,
        sep="\t",
        names=["ZFIN_ID", "SO_ID", "ZFIN_SYMBOL", "ENSEMBL_ID"],
    ).drop(columns="SO_ID")

    ncbi_df = pd.read_csv(
        ncbi_file,
        sep="\t",
        names=["ZFIN_ID", "SO_ID", "ZFIN_SYMBOL", "NCBI_ID"],
    ).drop(columns="SO_ID")

    uniport_df = pd.read_csv(
        uniprot_file,
        sep="\t",
        names=["ZFIN_ID", "SO_ID", "ZFIN_SYMBOL", "UNIPROT_ID"],
    ).drop(columns="SO_ID")

    dfs = list(
        replace_names(df.apply(lambda x: x.astype(str).str.strip().str.lower()), "ZFIN_SYMBOL", GENE_SET, ALIAS_MAP)
        for df in [zfin_df, ensembl_df, ncbi_df, uniport_df]
    )

    merged = reduce(
        lambda left, right: pd.merge(
            left, right, on=["ZFIN_ID", "ZFIN_SYMBOL"], how="left"
        ),
        dfs,
    )
    return merged.apply(lambda x: x.astype(str).str.lower())


def get_human_orthos():
    """This file reads the human orthology file into DataFrame with specific names and cleans it"""
    human_orthos = pd.read_csv(
        "Data/human_orthos.txt",
        sep="\t",
        names=[
            "ZFIN_ID",
            "ZFIN_SYMBOL",
            "ZFIN_NAME",
            "HUMAN_SYMBOL",
            "HUMAN_NAME",
            "OMIM_ID",
            "GENE_ID",
            "HGNC_ID",
            "EVIDENCE",
            "PUB_ID",
        ],
    )
    human_orthos = human_orthos.drop(
        columns=[
            "ZFIN_NAME",
            "ZFIN_ID",
            "HUMAN_NAME",
            "OMIM_ID",
            "HGNC_ID",
            "EVIDENCE",
            "PUB_ID",
        ]
    )

    human_orthos = replace_names(human_orthos, "ZFIN_SYMBOL", GENE_SET, ALIAS_MAP)
    human_orthos = human_orthos.rename(columns={"GENE_ID": "HUMAN_GENE_ID"})

    return human_orthos.apply(lambda x: x.astype(str).str.lower())


# FROM IID DATABASE
def ligand_receptor_pairs(annotated_file):
    """This function reads ligand-receptor data from IID into DataFrame and cleans it"""
    df = pd.read_csv(annotated_file, sep="\t").apply(
        lambda x: x.astype(str).str.lower()
    )
    return df.drop_duplicates().reset_index(drop=True)


def string_ligand_receptors(string_file="Data/string-protein-links.tsv", string_info_file="Data/string-protein-info.tsv", ensembl_file="Data/ensembl_assoc.txt"):
    """This function reads ligand-receptor data from STRING
    
    It parses two STRING files and creates a dictionary in form:
    key - (prot1, prot2) | (prot2, prot1)
    value - STRING physical score
    """
    ensembl_df = (
        pd.read_csv(
            ensembl_file,
            sep="\t",
            names=["ZFIN_ID", "SO_ID", "ZFIN_SYMBOL", "ENSEMBL_ID"],
        )
        .apply(lambda x: x.astype(str).str.strip().str.lower())
        .drop(columns=["SO_ID", "ZFIN_ID"])
        .drop_duplicates()
        .reset_index(drop=True)
    )
    # Ensembl_ID: Zfin_Symbol (gene name)
    ensembl_to_name = {line[1]: line[0] for line in ensembl_df.to_records(index=False)}

    with open(string_info_file) as fh_al:
        # 7955.ensdarg...: preferred_gene_name
        aliases = dict()
        for line in fh_al:
            data = line.split("\t")
            prot = data[0].strip().lower()  # 7955.ensdarg...
            pref_name = data[1].strip().lower()  # gene_name
            aliases[prot] = (
                ensembl_to_name.get(pref_name, pref_name)
                if pref_name.startswith("ensdarg")
                else pref_name
            )

    d = dict()
    df = pd.read_csv(string_file, sep="\t").apply(lambda x: x.astype(str).str.strip().str.lower())
    for record in df.to_records(index=False):
        rec0 = record[0]
        rec1 = record[1]

        prot1 = aliases.get(rec0, rec0)  # Gets ENS_ID from aliases dict
        prot2 = aliases.get(rec1, rec1)  # Same as above

        prot1 = check_in_map(prot1, GENE_SET, ALIAS_MAP)
        prot2 = check_in_map(prot2, GENE_SET, ALIAS_MAP)

        key = (prot1, prot2)
        rev = (prot2, prot1)
        score = float(record[2])

        d[key] = score
        d[rev] = score
    return d

def agg_frame(df: pd.DataFrame, organism: str = "HUMAN"):
    """This function reduces the number of records in a DataFrame
    by aggregating all the UNIPROT, NCBI and ENSEMBL ID's into a single string (into single Excel cell)"""
    uniprot_ids = (
        df.groupby(f"{organism}_SYMBOL")
        .agg(
            {
                "UNIPROT_ID": lambda x: ";".join(set(str(y) for y in x.unique())),
                "NCBI_ID": lambda x: ";".join(set(str(y) for y in x.unique())),
                "ENSEMBL_ID": lambda x: ";".join(set(str(y) for y in x.unique())),
            }
        )
        .reset_index()
    )

    df = (
        df.drop(columns={"UNIPROT_ID", "NCBI_ID", "ENSEMBL_ID"})  # "Evidence_type"
        .drop_duplicates()
        .reset_index(drop=True)
    )
    df = reduce(
        lambda left, right: left.merge(right, on=f"{organism}_SYMBOL"),
        [df, uniprot_ids],  # evidences
    )
    return df.drop_duplicates().reset_index(drop=True)


if __name__ == "__main__":
    # Loads ZFIN data and human orthologs
    logging.info("Fetching ZF and HUMAN orthology data...")
    zfin_df = get_zfin_data()
    human_ort = get_human_orthos()
    
    # Loads zebrafish ligands and receptors from assets
    logging.info("Loading receptors and ligands...")
    zebrafish_ligands = pd.read_excel("Assets/zebrafish_ligands.xlsx").apply(
        lambda x: x.astype(str).str.lower()
    )
    zebrafish_receptors = pd.read_excel("Assets/zebrafish_receptors.xlsx").apply(
        lambda x: x.astype(str).str.lower()
    )

    # Joins receptors and ligands on ensembl, ncbi and uniprot data
    logging.info("Joining association files...")
    zebrafish_ligands = replace_names(zebrafish_ligands, "Symbol", GENE_SET, ALIAS_MAP)
    zebrafish_receptors = replace_names(
        zebrafish_receptors, "Symbol", GENE_SET, ALIAS_MAP
    )

    zebrafish_ligands = zebrafish_ligands.merge(
        zfin_df.rename(columns={"ZFIN_SYMBOL": "Symbol"}), on="Symbol", how="left"
    ).rename(columns={"Symbol": "ZFIN_SYMBOL"})
    zebrafish_ligands = zebrafish_ligands.drop_duplicates().reset_index(drop=True)

    zebrafish_receptors = zebrafish_receptors.merge(
        zfin_df.rename(columns={"ZFIN_SYMBOL": "Symbol"}), on="Symbol", how="left"
    ).rename(columns={"Symbol": "ZFIN_SYMBOL"})
    zebrafish_receptors = zebrafish_receptors.drop_duplicates().reset_index(drop=True)

    # Reduces number of records by aggregating various IDs
    logging.info("Reducing redundant records; aggregating frames...")
    zebrafish_ligands = agg_frame(zebrafish_ligands, "ZFIN")
    zebrafish_receptors = agg_frame(zebrafish_receptors, "ZFIN")

    # Joins human orthology df on zfin genes
    logging.info("Joining human orthology...")
    zfin_human_receptors = (
        zebrafish_receptors.merge(human_ort, on=["ZFIN_SYMBOL"], how="left")
        .drop_duplicates()
        .reset_index(drop=True)
    )
    zfin_human_ligands = (
        zebrafish_ligands.merge(human_ort, on=["ZFIN_SYMBOL"], how="left")
        .drop_duplicates()
        .reset_index(drop=True)
    )

    logging.info("Loading Matrisome annotation")
    parsed = rdata.parser.parse_file("Data/matrisome.list.rda")
    converted = rdata.conversion.convert(parsed)

    matrisome_zfin = converted['matrisome.list']['zebrafish'].apply(lambda x: x.astype(str).str.lower().str.strip())
    matrisome_human = converted['matrisome.list']['human'].apply(lambda x: x.astype(str).str.lower().str.strip())

    zfin_families = dict()
    human_families = dict()

    for pair in matrisome_zfin[["gene", "family"]].to_records(index=False):
        zfin_families[pair[0]] = pair[1]
    for pair in matrisome_human[["gene", "family"]].to_records(index=False):
        human_families[pair[0]] = pair[1]

    logging.info("Adding Matrisome annotation")
    zfin_human_receptors["ZFIN_Matrisome_annotation"] = (
        zfin_human_receptors["ENSEMBL_ID"].map(lambda x: zfin_families.get(x))
    )

    zfin_human_receptors["Human_Matrisome_annotation"] = (
        zfin_human_receptors["HUMAN_SYMBOL"].map(lambda x: human_families.get(x))
    )

    zfin_human_ligands["ZFIN_Matrisome_annotation"] = (
        zfin_human_ligands["ENSEMBL_ID"].map(lambda x: zfin_families.get(x))
    )

    zfin_human_ligands["Human_Matrisome_annotation"] = (
        zfin_human_ligands["HUMAN_SYMBOL"].map(lambda x: human_families.get(x))
    )

    # Loads genome alliance file, filters and adds conservation scores
    logging.info("Adding conservation scores")
    genome_alliance = pd.read_table("Data/conservation.tsv", comment="#")
    genome_alliance = genome_alliance.loc[(genome_alliance["Gene1SpeciesName"] == "Danio rerio") & (genome_alliance["Gene2SpeciesName"] == "Homo sapiens")]
    genome_alliance = genome_alliance.loc[(genome_alliance["Algorithms"].str.contains("ZFIN")) & (genome_alliance["IsBestScore"] == "Yes")]
    genome_alliance = genome_alliance[["Gene1Symbol", "Gene2Symbol", "Algorithms", "AlgorithmsMatch"]].drop_duplicates(["Gene1Symbol"]).reset_index(drop=True)
    genome_alliance["Gene2Symbol"] = genome_alliance["Gene2Symbol"].str.lower()

    zfin_human_receptors = zfin_human_receptors.merge(genome_alliance, left_on = ["ZFIN_SYMBOL", "HUMAN_SYMBOL"], right_on = ["Gene1Symbol", "Gene2Symbol"], how = "left").drop(columns=["Gene1Symbol", "Gene2Symbol"])
    zfin_human_ligands = zfin_human_ligands.merge(genome_alliance, left_on = ["ZFIN_SYMBOL", "HUMAN_SYMBOL"], right_on = ["Gene1Symbol", "Gene2Symbol"], how = "left").drop(columns=["Gene1Symbol", "Gene2Symbol"])

    zfin_human_receptors["AlgorithmsMatch"] = zfin_human_receptors["AlgorithmsMatch"].fillna(0)
    zfin_human_ligands["AlgorithmsMatch"] = zfin_human_ligands["AlgorithmsMatch"].fillna(0)

    zfin_human_receptors.rename(columns={"AlgorithmsMatch": "ConservationScore"}, inplace = True)
    zfin_human_ligands.rename(columns={"AlgorithmsMatch": "ConservationScore"}, inplace = True)


    # Reads all potential pairs for human
    logging.info("Loading data from IID...")
    df_pairs_human = ligand_receptor_pairs("Data/human_annotated_PPIs_cut.txt")

    logging.info("Creating pairs from IID data...")
    human_protein_pairs = dict()
    for record in df_pairs_human.to_records(index=False):
        key = (record[0], record[1])
        rev_key = (record[1], record[0])
        
        val = record[2]
        
        human_protein_pairs[key] = val
        human_protein_pairs[rev_key] = val

    # Creates all possible receptor-ligand pair
    logging.info("Creating every possible ligand-receptor pair...")
    human_every_pair = zfin_human_ligands.merge(
        zfin_human_receptors, how="cross", suffixes=["_LIG", "_REC"]
    )

    # Checks which of the constructed pairs are in the list of all possible interactions from IID and STRING.
    logging.info("Creating list of pairs from STRING...")
    string_zfin_pairs = string_ligand_receptors()

    logging.info("Creating list of pairs from CellTalkDB")
    celltalkdf = pyreadr.read_r("Data/LRdb.rda")["LRdb"].apply(lambda x: x.astype(str).str.strip().str.lower())
    celltalkdf = celltalkdf[celltalkdf["species"] == "human"]

    celltalk_pairs = set()
    for record in celltalkdf.to_records(index=False):
        key = (record[0], record[1])
        rev_key = (record[1], record[0])
        
        celltalk_pairs.add(key)
        celltalk_pairs.add(rev_key)

    logging.info("Creating list of pairs from CellCellInteractions")
    cellcelldf = pd.read_table("Data/cellcellinteractions.txt").apply(lambda x: x.astype(str).str.strip().str.lower())
    cellcell_pairs = set()
    for record in cellcelldf.to_records(index=False):
        key = (record[0], record[1])
        rev_key = (record[1], record[0])
        
        cellcell_pairs.add(key)
        cellcell_pairs.add(rev_key)

    logging.info("Creating list of pairs from Ramilowski dataset")
    ramilowskidf = pd.read_table("Data/ramilowski.txt").apply(lambda x: x.astype(str).str.strip().str.lower())
    ramilowski_pairs = set()
    for record in ramilowskidf.to_records(index=False):
        key = (record[0], record[1])
        rev_key = (record[1], record[0])
        
        ramilowski_pairs.add(key)
        ramilowski_pairs.add(rev_key)

    logging.info("Creating list of pairs from CellPhoneDB dataset")
    cellphonedf = pd.read_csv("Data/cellphonedb.csv")[["protein_name_a", "protein_name_b", "is_ppi"]].apply(lambda x: x.astype(str).str.strip().str.lower())
    cellphonedf = cellphonedf.loc[(cellphonedf["protein_name_a"] != "nan") & (cellphonedf["protein_name_b"] != "nan") & (cellphonedf["is_ppi"] == "true")]
    cellphone_pairs = set()

    for record in cellphonedf.to_records(index=False):
        key = (record[0], record[1])
        rev_key = (record[1], record[0])
        
        cellphone_pairs.add(key)
        cellphone_pairs.add(rev_key)

    logging.info("Creating list of dimers from CellPhoneDB")
    cellphone_partners_df = pd.read_csv("Data/cellphonedb.csv")[["protein_name_a", "protein_name_b", "partner_b", "is_ppi"]].apply(lambda x: x.astype(str).str.strip().str.lower())
    cellphone_partners_df = cellphone_partners_df.loc[(cellphone_partners_df["protein_name_a"] != "nan") & (cellphone_partners_df["protein_name_b"] == "nan") & (cellphone_partners_df["partner_b"] != "nan") & (cellphone_partners_df["is_ppi"] == "true")]
    cellphone_partners_df = cellphone_partners_df.loc[~((cellphone_partners_df["partner_b"].str.contains("enhancer")) | (cellphone_partners_df["partner_b"].str.contains("receptor")) | (cellphone_partners_df["partner_b"].str.contains("complex")))]
    cellphone_partners = set()

    l = cellphone_partners_df[["protein_name_a", "partner_b"]].to_records(index=False).tolist()
    for p in l:
        g = p[0]
        pair = tuple(p[1].split("_"))
        if len(pair) != 2:
            continue
        cellphone_partners.add((g, pair[0]))
        cellphone_partners.add((pair[0], g))
        
        cellphone_partners.add((g, pair[1]))
        cellphone_partners.add((pair[1], g))

    logging.info("Creating list of pairs from CellChat")
    parsed = rdata.parser.parse_file("Data/CellChatDB.human.rda")
    converted = rdata.conversion.convert(parsed)["CellChatDB.human"]["interaction"][["ligand", "receptor", "annotation"]].apply(lambda x: x.astype(str).str.lower().str.strip())

    cellchat = dict()
    cellchat_complex = dict()
    for pair in converted.to_records():
        idx = pair[0].lower()
        ligand, receptor = pair[1], pair[2]
        annotation = pair[3]
        genes = idx.split("_")
        try:
            if "_" in receptor:
                key1, key2 = (genes[0], genes[1]), (genes[0], genes[2])
                rev_key1, rev_key2 = (genes[1], genes[0]), (genes[2], genes[0])
                cellchat_complex[key1] = annotation
                cellchat_complex[key2] = annotation
                cellchat_complex[rev_key1] = annotation
                cellchat_complex[rev_key2] = annotation
            elif "_" in ligand:
                key1, key2 = (genes[2], genes[0]), (genes[2], genes[1])
                rev_key1, rev_key2 = (genes[0], genes[2]), (genes[1], genes[2])
                cellchat_complex[key1] = annotation
                cellchat_complex[key2] = annotation
                cellchat_complex[rev_key1] = annotation
                cellchat_complex[rev_key2] = annotation
            else:
                key = (genes[0], genes[1])
                rev_key = (genes[1], genes[0])
                cellchat[key] = annotation
                cellchat[rev_key] = annotation
        except IndexError:
            pass


    # Computing section
    logging.info("Computing pairs based on STRING data...")
    human_every_pair["exists_STRING"] = (
        human_every_pair[["ZFIN_SYMBOL_LIG", "ZFIN_SYMBOL_REC"]]
        .apply(tuple, 1)
        .isin(string_zfin_pairs)
    )

    logging.info("Computing pairs based on IID data...")
    human_every_pair["exists"] = (
        human_every_pair[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]]
        .apply(tuple, 1)
        .isin(human_protein_pairs)
    )

    logging.info("Computing pairs based on CellTalkDB data...")
    human_every_pair["exists_CellTalk"] = (
            human_every_pair[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]]
            .apply(tuple, 1)
            .isin(celltalk_pairs)
    )

    logging.info("Computing pairs based on CellCellInteractions data...")
    human_every_pair["exists_CellCellInteractions"] = (
            human_every_pair[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]]
            .apply(tuple, 1)
            .isin(cellcell_pairs)
    )

    logging.info("Computing pairs based on Ramilowski dataset...")
    human_every_pair["exists_Ramilowski"] = (
            human_every_pair[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]]
            .apply(tuple, 1)
            .isin(ramilowski_pairs)
    )

    logging.info("Computing pairs based on CellPhoneDB dataset...")
    human_every_pair["exists_CellPhone"] = (
            human_every_pair[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]]
            .apply(tuple, 1)
            .isin(cellphone_pairs)
    )

    logging.info("Computing pairs based on CellPhoneDB Dimer dataset...")
    human_every_pair["exists_Dimer"] = (
            human_every_pair[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]]
            .apply(tuple, 1)
            .isin(cellphone_partners)
    )

    logging.info("Computing pairs based on CellChat data...")
    human_every_pair["exists_CellChat"] = (
        human_every_pair[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]]
        .apply(tuple, 1)
        .isin(cellchat.keys())
    )
    logging.info("Computing pairs based on CellChat dimer data...")
    human_every_pair["exists_CellChat_dimer"] = (
        human_every_pair[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]]
        .apply(tuple, 1)
        .isin(cellchat_complex.keys())
    )


    logging.info("Filtering results")
    human_every_pair = human_every_pair.loc[
        (human_every_pair["exists"]) | \
            (human_every_pair["exists_STRING"]) | \
                (human_every_pair["exists_CellTalk"]) | \
                    (human_every_pair["exists_CellCellInteractions"]) | \
                      (human_every_pair["exists_Ramilowski"]) | \
                        (human_every_pair["exists_CellPhone"]) | \
                            (human_every_pair["exists_Dimer"]) | \
                                (human_every_pair["exists_CellChat"]) | \
                                    (human_every_pair["exists_CellChat_dimer"])
    ]

    logging.info("Adding IID evidences")
    human_every_pair["IID_EVIDENCE"] = (
        human_every_pair[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]]
        .apply(tuple, 1)
        .map(lambda x: human_protein_pairs.get(x))
    )
    logging.info("Adding STRING scores")
    human_every_pair["PHYSICAL_INTERACTION_SCORE"] = (
        human_every_pair[["ZFIN_SYMBOL_LIG", "ZFIN_SYMBOL_REC"]]
        .apply(tuple, 1)
        .map(lambda x: string_zfin_pairs.get(x))
    )

    logging.info("Adding CellTalk yes/no annotation")
    human_every_pair["CellTalkDB"] = (
        human_every_pair["exists_CellTalk"]
        .map(lambda x: "Yes" if x else "No")
    )

    logging.info("Adding CellCellInteractions yes/no annotation")
    human_every_pair["CellCellInteractions"] = (
        human_every_pair["exists_CellCellInteractions"]
        .map(lambda x: "Yes" if x else "No")
    )

    logging.info("Adding Ramilowski yes/no annotation")
    human_every_pair["Ramilowski"] = (
        human_every_pair["exists_Ramilowski"]
        .map(lambda x: "Yes" if x else "No")
    )

    logging.info("Adding CellPhoneDB yes/no annotation")
    human_every_pair["CellPhoneDB"] = (
        human_every_pair["exists_CellPhone"]
        .map(lambda x: "Yes" if x else "No")
    )

    logging.info("Adding CellPhone dimer annotation...")
    human_every_pair["CellPhoneDB_Dimer_or_Complex"] = (
        human_every_pair["exists_Dimer"]
        .map(lambda x: "LR signaling might require dimerization or complex" if x else None)
    )

    logging.info("Adding CellChat annotation")
    human_every_pair["CellChat"] = (
        human_every_pair["exists_CellChat"]
        .map(lambda x: "Yes" if x else "No")
    )

    human_every_pair["CellChat_annotation"] = (
        human_every_pair[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]]
        .apply(tuple, 1)
        .map(lambda x: cellchat.get(x, cellchat_complex.get(x)))
    )

    logging.info("Adding CellChat dimer annotation") 
    human_every_pair["CellChat_Dimer_or_Complex"] = (
        human_every_pair["exists_CellChat_dimer"]
        .map(lambda x: "LR signaling might require dimerization or complex" if x else None)
    )

    logging.info("Dropping unused columns")
    human_every_pair = (
        human_every_pair.drop(
            columns={
                "exists",
                "exists_STRING",
                "exists_CellTalk",
                "exists_CellCellInteractions",
                "exists_Ramilowski", 
                "exists_CellPhone",
                "exists_Dimer",
                "exists_CellChat",
                "exists_CellChat_dimer"
            }
        )
        .drop_duplicates()
        .reset_index(drop=True)
    )


    logging.info("Loading ligand-type annotations")
    ligand_type_zfin = pd.read_table("Data/zfin.gaf", comment="!", names=["Gene", "GO"]).apply(lambda x: x.astype(str).str.lower()).to_records(index=False).tolist()
    ligand_type_human = pd.read_table("Data/goa_human.gaf", comment="!", names=["Gene", "GO"]).apply(lambda x: x.astype(str).str.lower()).to_records(index=False).tolist()

    zfin_ligands = defaultdict(set)
    human_ligands = defaultdict(set)
    zfin_secreted = set()
    human_secreted = set()

    for ligand in ligand_type_zfin:
        zfin_ligands[ligand[0]].add(ligand[1])

    for ligand in ligand_type_human:
        human_ligands[ligand[0]].add(ligand[1])

    # Rule 1: for secreted ligands: 
    # (Our pair database Ligand gene or human ligand gene AND (((go:0005615 OR GO:0005576) NOT (go:0005886 OR go:0009986)))

    for ligand in zfin_ligands:
        values = zfin_ligands[ligand]
        condition = ("go:0005615" in values or "go:0005576" in values) and ("go:0005886" not in values and "go:0009986" not in values)
        if condition:
            zfin_secreted.add(ligand)

    for ligand in human_ligands:
        values = zfin_ligands[ligand]
        condition = ("go:0005615" in values or "go:0005576" in values) and ("go:0005886" not in values and "go:0009986" not in values)
        if condition:
            human_secreted.add(ligand)
    
    logging.info("Adding ligand-type annotations")
    human_every_pair["Ligand_type"] = (
        human_every_pair[["ZFIN_SYMBOL_LIG", "HUMAN_SYMBOL_LIG"]]
        .apply(tuple, 1)
        .map(lambda x: "Secreted" if x[0] in zfin_secreted or x[1] in human_secreted else None)
    )

    final_pairs = human_every_pair
    final_pairs[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]] = final_pairs[
        ["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]
    ].apply(lambda x: x.str.upper().str.strip())
    
    logging.info("Exporting files")

    final_column_order = [
        "ZFIN_SYMBOL_LIG",
        "ZFIN_SYMBOL_REC",
        "PHYSICAL_INTERACTION_SCORE",
        "HUMAN_SYMBOL_LIG",
        "HUMAN_SYMBOL_REC",
        "IID_EVIDENCE",
        "CellCellInteractions"
        "CellChat",
        "CellChat_Dimer_or_Complex",
        "CellPhoneDB",
        "CellPhoneDB_Dimer_or_Complex",
        "CellTalkDB",
        "Ramilowski",
        "CellChat_annotation",
        "Ligand_type",
        "ZFIN_Matrisome_annotation_LIG",
        "Human_Matrisome_annotation_LIG"
        "ZFIN_Matrisome_annotation_REC",
        "Human_Matrisome_annotation_REC",
        "ZFIN_ID_LIG",
        "UNIPROT_ID_LIG",
        "NCBI_ID_LIG",
        "ENSEMBL_ID_LIG",
        "HUMAN_GENE_ID_LIG",
        "Algorithms_LIG",
        "ConservationScore_LIG",
        "ZFIN_ID_REC",
        "UNIPROT_ID_REC",
        "NCBI_ID_REC",
        "ENSEMBL_ID_REC",
        "HUMAN_GENE_ID_REC",
        "Algorithms_REC",
        "ConservationScore_REC"
    ]

    db_name = "Database"
    conn = sa.create_engine(f"sqlite:///Database/{db_name}.sqlite")
    
    final_pairs.to_sql("PAIRS", conn, index=False, if_exists="replace")
    final_pairs.to_csv(f"Database/{db_name}.csv", index=False)
    final_pairs.to_excel(f"Database/{db_name}.xlsx", index=False)
    
    logging.info(f"Saved to {db_name}.sqlite and {db_name}.xlsx/csv")
