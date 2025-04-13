import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


def valid_smiles(smiles: str) -> int:
    """
    Validates a given SMILES string by attempting to canonicalize it, create a molecular representation,
    embed it in 3D space, and compute its XYZ coordinates. Returns a validation status based on the
    success of these operations.

    :param smiles: A SMILES string representing a chemical structure.
    :type smiles: str
    :return: Returns 1 if the SMILES string is valid and successfully processed; otherwise, returns 0.
    :rtype: int
    """
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


def valid_smiles_multiple(f: pd.DataFrame) -> tuple[bool, str]:
    """
    Processes a DataFrame to validate and embed molecular SMILES strings, generating
    3D structures in XYZ format if valid. The function checks if each SMILES string can
    be canonicalized, embedded, and converted into an XYZ block. The SMILES and their
    corresponding XYZ blocks are stored in a dictionary, and the function indicates
    whether the process was successful for all entries in the DataFrame.

    :param f: A pandas DataFrame where SMILES strings are located in the second column.
    :return: A tuple indicating whether the validation and embedding process is
     successful for all molecules, and either a success message or the first SMILES
     string that failed validation.
    :rtype: tuple
    """
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


def check_size_multiple(f: pd.DataFrame) -> int:
    """
    Analyzes a given DataFrame to compute the maximum number of atoms in molecular
    structures represented by SMILES strings. The second column of the DataFrame
    is used to extract SMILES strings, and an additional feature name from the
    header is appended to the list of SMILES for analysis. If no valid molecules
    can be constructed from the SMILES strings, the function will return 0.

    :param f: A pandas DataFrame where the second column contains SMILES strings
        for molecular representation. The header of the second column is also
        included in the analysis.
    :type f: pandas.DataFrame
    :return: Maximum number of atoms in molecules represented by SMILES strings
        in the provided DataFrame, or 0 if no valid molecular structures are
        constructed.
    :rtype: int
    """

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


def get_image_data(smiles: str, sites: list[int]) -> str:
    """
    Generates a molecule representation in SVG format based on the provided
    SMILES string and a list of atomic sites to annotate.

    :param smiles: A SMILES (Simplified Molecular Input Line Entry System)
        string representing the molecular structure.
    :type smiles: str
    :param sites: A list of atomic sites (1-based indices) that will be
        labeled on the molecular structure.
    :type sites: list[int]
    :return: A string containing the SVG representation of the molecule
        with annotated atomic sites.
    :rtype: str
    """

    canonical_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
    mol = Chem.MolFromSmiles(canonical_smiles)
    d = rdMolDraw2D.MolDraw2DSVG(250, 250)
    for atom in sites:
        mol.GetAtomWithIdx(atom - 1).SetProp("_displayLabel", str(atom))
    d.DrawMolecule(mol)
    d.FinishDrawing()
    molecule = d.GetDrawingText()
    return molecule