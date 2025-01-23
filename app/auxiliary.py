from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


def valid_smiles(smiles):
     # TODO RA: Docstring. Would be good to know what the return values mean. 
    try:
        canon = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        mol = Chem.MolFromSmiles(canon)
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        xyz = Chem.MolToXYZBlock(mol)
        if len(xyz) == 0:
            return 0
        else:
            return 1
    except:
        return 0


def valid_smiles_multiple(f):
    # TODO RA: Docstring. Would be good to know what the return values mean. 
    smiles = f.iloc[:, 1].tolist()
    smiles.append(list(f)[1])
    smile_xyz_dict = {}
    for smi in smiles:
        try:
            canon = Chem.MolToSmiles(Chem.MolFromSmiles(smi))
            mol = Chem.MolFromSmiles(canon)
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)
            xyz = Chem.MolToXYZBlock(mol)
            smile_xyz_dict[smi] = xyz

        except:
            return False, smi
    for key, value in smile_xyz_dict.items():
        if value:
            continue
        else:
            return False, key
    return True, 'okay'


def check_size(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return len(mol.GetAtoms())
    else:
        return 0


def check_size_multiple(f):
    # TODO RA: Docstring. Would be good to know what the return values mean. 

    smiles = f.iloc[:, 1].tolist()
    smiles.append(list(f)[1])
    smiles_length = []
    for smi in smiles:
        mol = Chem.MolFromSmiles(smi)
        try:
            smiles_length.append(len(mol.GetAtoms()))
        except:
            continue
    if smiles_length:
        return max(smiles_length)
    else:
        return 0


def get_image_data(smiles, sites):
    # TODO RA: Docstring. Would be good to know what the return values mean. 

    canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    mol = Chem.MolFromSmiles(canonical_smiles)
    d = rdMolDraw2D.MolDraw2DSVG(250, 250)
    for atom in sites:
        mol.GetAtomWithIdx(atom - 1).SetProp("_displayLabel", str(atom))
    d.DrawMolecule(mol)
    d.FinishDrawing()
    molecule = d.GetDrawingText()
    return molecule