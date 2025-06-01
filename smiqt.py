import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy


def ligand_prep(ligmol, out):
    '''
    input:
    :ligmol: RDkit mol
    :out: output path and name
    output:
    :pdbqt: writes a pdbfile to the working directory
    '''
    mol_Hs = Chem.AddHs(ligmol, addCoords=True)
    AllChem.EmbedMolecule(mol_Hs, AllChem.ETKDG())
    AllChem.MMFFOptimizeMolecule(mol_Hs)
    preparator = MoleculePreparation(merge_these_atom_types=("H",))
    mol_setup = preparator.prepare(mol_Hs)[0]
    pdbqt_string, is_ok, error_msg = PDBQTWriterLegacy.write_string(mol_setup)
    if is_ok:
            #print(type(pdbqt_string) ,pdbqt_string, end="")
        with open(f'{out}.pdbqt', 'w') as pdbqt:
            pdbqt.writelines(pdbqt_string)
        return 1
    if not is_ok:
        return 0, str(error_msg)



print("This script converts smiles in a text or csv file to pdbqt files in folders")
text_or_csv = int(input("Are SMILES in a text file or csv? Enter 0 or 1 respectively: "))


SMILES_COL = 'SMILES'

if text_or_csv == 1:
    df = pd.read_csv(r"./ligands.csv")
    df['mol'] = df[SMILES_COL].apply(Chem.MolFromSmiles)
    mol_index = list(df.columns).index('mol')
    failed = []
    for idx,row in enumerate(df.values):
        os.mkdir(f"{idx}")
        try:
            ligand_prep(row[mol_index], out=rf"{idx}/l")
        except Exception as e:
            os.removedirs(f"{idx}")
            failed.append(f"{idx}")
elif text_or_csv == 0:
    print('text format is not availible yet')



