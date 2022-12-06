#!/usr/bin/env python3
import pandas as pd
import sqlalchemy as sa
from functools import reduce


def read_gene_map(file: str = "Data/aliases.txt"):
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
    if x in gene_set:
        return x
    elif x in alias_map:
        return alias_map[x]
    else:
        return x


def replace_names(
    df: pd.DataFrame, column: str, gene_set: set, alias_map: dict[str, str]
):
    series = df[column].str.lower().str.strip()
    new_series = series.apply(check_in_map, gene_set=gene_set, alias_map=alias_map)
    df[column] = new_series
    return df

GENE_SET, ALIAS_MAP = read_gene_map()


def get_zfin_data(
    zfin_file="Data/zfin_genes.txt",
    ensembl_file="Data/ensembl_assoc.txt",
    ncbi_file="Data/ncbi_assoc.txt",
    uniprot_file="Data/uni_assoc.txt",
):
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
    df = pd.read_csv(annotated_file, sep="\t").apply(
        lambda x: x.astype(str).str.lower()
    )
    return df.drop_duplicates().reset_index(drop=True)


def string_ligand_receptors(
    string_file="Data/string-protein-links.tsv",
    string_info_file="Data/string-protein-info.tsv",
    ensembl_file="Data/ensembl_assoc.txt",
):
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
    # ZFIN GENES ===========
    # Obtain zfin data and orthology to mouse and human
    print("Fetching ZF and HUMAN orthology data...")
    zfin_df = get_zfin_data()
    human_ort = get_human_orthos()
    # Load zebrafish ligands and receptors symbol
    print("Loading receptors and ligands...")
    zebrafish_ligands = pd.read_excel("Assets/zebrafish_ligands.xlsx").apply(
        lambda x: x.astype(str).str.lower()
    )
    zebrafish_receptors = pd.read_excel("Assets/zebrafish_receptors.xlsx").apply(
        lambda x: x.astype(str).str.lower()
    )
    # Join receptors and ligands on ensembl, ncbi and uniprot data
    print("Joining association files...")
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

    print("Reducing redundant records; aggregating frames...")
    # Reduce number of records by aggregating various IDs
    zebrafish_ligands = agg_frame(zebrafish_ligands, "ZFIN")
    zebrafish_receptors = agg_frame(zebrafish_receptors, "ZFIN")

    # ==========
    # ORTHOLOGY CREATION
    # ==========
    # Join human orthology df on zfin genes
    print("Joining human orthology...")
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

    # Pair computing
    # Read all potential pairs for human
    print("Loading data from IID...")
    df_pairs_human = ligand_receptor_pairs("Data/human_annotated_PPIs_cut.txt")

    # Symbol1, symbol2, evidence_type
    print("Creating pairs from IID data...")
    human_protein_pairs = dict()
    for record in df_pairs_human.to_records(index=False):
        key = (record[0], record[1])
        rev_key = (record[1], record[0])
        
        val = record[2]
        
        human_protein_pairs[key] = val
        human_protein_pairs[rev_key] = val

    # Create all possible receptor-ligand pair
    print("Creating every possible ligand-receptor pair...")
    human_every_pair = zfin_human_ligands.merge(
        zfin_human_receptors, how="cross", suffixes=["_LIG", "_REC"]
    )

    # Check which of the constructed pairs are in the list of all possible interactions from IID and STRING.
    print("Creating list of pairs from STRING...")
    string_zfin_pairs = string_ligand_receptors()

    print("Computing pairs based on IID data...")
    human_every_pair["exists"] = (
        human_every_pair[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]]
        .apply(tuple, 1)
        .isin(human_protein_pairs)
    )

    print("Computing pairs based on STRING data...")
    human_every_pair["exists_STRING"] = (
        human_every_pair[["ZFIN_SYMBOL_LIG", "ZFIN_SYMBOL_REC"]]
        .apply(tuple, 1)
        .isin(string_zfin_pairs)
    )

    print("Filtering results")
    human_every_pair = human_every_pair.loc[
        (human_every_pair["exists"]) | (human_every_pair["exists_STRING"])
    ]

    print("Adding IID evidences")
    human_every_pair["IID_EVIDENCE"] = (
        human_every_pair[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]]
        .apply(tuple, 1)
        .map(lambda x: human_protein_pairs.get(x))
    )
    print("Adding STRING scores")
    human_every_pair["PHYSICAL_INTERACTION_SCORE"] = (
        human_every_pair[["ZFIN_SYMBOL_LIG", "ZFIN_SYMBOL_REC"]]
        .apply(tuple, 1)
        .map(lambda x: string_zfin_pairs.get(x))
    )

    print("Dropping unused columns")
    human_every_pair = (
        human_every_pair.drop(
            columns={
                "exists",
                "exists_STRING"
            }
        )
        .drop_duplicates()
        .reset_index(drop=True)
    )

    print("Exporting files")
    final_pairs = human_every_pair
    final_pairs[["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]] = final_pairs[
        ["HUMAN_SYMBOL_LIG", "HUMAN_SYMBOL_REC"]
    ].apply(lambda x: x.str.upper().str.strip())
    
    db_name = "Database"
    conn = sa.create_engine(f"sqlite:///Database/{db_name}.sqlite")
    
    final_pairs.to_sql("PAIRS", conn, index=False, if_exists="replace")
    final_pairs.to_csv(f"Database/{db_name}.csv", index=False)
    final_pairs.to_excel(f"Database/{db_name}.xlsx", index=False)
    
    print(f"Saved to {db_name}.sqlite and {db_name}.xlsx/csv")
