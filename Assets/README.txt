ListA_Receptors: 
    First Select CellCellinteractions Receptors Type: Receptor AND ECM/RECEPTOR
    Next, to the above list, remove all CellTalkDB Human ligands
    Next, to the above list, combine with CellTalkDB Human Receptors, remove duplicates. 
    Next, to the above list, compare with Human Matrisome (without retired genes), and manually curate to remove ligands. 
    ANXA2 can act as both ligand for ROBO4 and receptor for tenascin-C

ListB_Ligands:
    CellCellinteractions Ligands Type: Ligand AND Receptor/Ligand AND ECM/Ligand AND ECM/Receptor/Ligand
    Next, to the above list, remove all ListA genes
    Next, to the above list, combine with CellTalkDB Human ligands
    Next, to the above list, combine with human matrisome without retired genes and without contaminating receptors

Zebrafish matrisome:
    Remove retired genes
    Remove all ListA genes