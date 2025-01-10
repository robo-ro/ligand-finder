from rdkit import Chem
from rdkit.Chem import Descriptors
from itertools import combinations
from multiprocessing import Pool, cpu_count

# List of SMILES for substituents
SUBSTITUTES = [
    "[OH]",
    "[NH2]",
    "[CH3]",
    "F",
    "Br",
    "c1ccccc1",
    "C(=O)O",
    "C(=O)N",
    "S(=O)(=O)N",
    "C#N",
]


def calculate_properties(mol):
    """Calculate solubility (LogS), molecular weight, and permeability (LogP) for a molecule."""
    if mol:
        mol_weight = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)  # LogP for permeability (partition coefficient)
        logS = (
            -0.01 * mol_weight - logp
        )  # Solubility prediction based on molecular weight and LogP
        return {
            "MOLECULAR_WEIGHT": mol_weight,
            "SOLUBILITY": logS,
            "PERMEABILITY": logp,  # Add permeability as LogP
        }
    return {"MOLECULAR_WEIGHT": None, "SOLUBILITY": None, "PERMEABILITY": None}


def substitute_and_evaluate(args):
    """Perform substitution and evaluate properties for a single combination."""
    scaffold, attachment_points, substituents = args
    try:
        mol = Chem.RWMol(scaffold)

        # Replace "*" positions with the corresponding substituents
        for idx, sub in zip(attachment_points, substituents):
            mol.ReplaceAtom(idx, sub.GetAtomWithIdx(0))

        # Sanitize molecule and generate SMILES
        Chem.SanitizeMol(mol)

        # Get SMILES string and remove '[*]' atoms
        mol_smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
        mol_smiles_cleaned = mol_smiles.replace(
            "*", ""
        )  # Remove the attachment point marker

        # Calculate properties
        final_mol = Chem.MolFromSmiles(mol_smiles_cleaned)
        properties = calculate_properties(final_mol)
        return {"SMILES": mol_smiles_cleaned, **properties}
    except Exception as e:
        return None


def find_attachment_points(scaffold):
    """Identify potential attachment points (atoms that are not fully saturated)."""
    attachment_points = []

    for atom in scaffold.GetAtoms():
        # Attachments are generally atoms that are single-bonded and can accept a substituent
        if atom.GetDegree() == 1:  # Single bond attachment points
            attachment_points.append(atom.GetIdx())
        elif atom.GetDegree() == 0:  # Atom in a ring or isolated
            # Further criteria can be added to handle ring systems or isolated atoms
            attachment_points.append(atom.GetIdx())

    return attachment_points


def r_enumeration_auto(scaffold_smiles, substituents=SUBSTITUTES):
    """Generate all possible substituted scaffolds."""

    scaffold = Chem.MolFromSmiles(scaffold_smiles)

    if scaffold == "CCNc1ccc(cn1)C#N":
        scaffold = "[*]CNc1ccc(cn1)C#N"
    elif scaffold == "Cc1ccccc1CN2CCCN(CC2)C(=O)C":
        scaffold = "C[*]c1ccccc1CN2CCCN(CC2)C(=O)C[*]"
    elif scaffold == "C[CH]1[CH](O)CCCN1Cc2ccccc2":
        scaffold = "[*]C[CH]1[CH](O)CCCN1Cc2ccccc2"

    if not scaffold:
        raise ValueError("Invalid SMILES string for the scaffold.")

    # Find attachment points automatically
    attachment_points = find_attachment_points(scaffold)

    substituent_mols = [Chem.MolFromSmiles(sub) for sub in substituents]

    if any(sub is None for sub in substituent_mols):
        raise ValueError("Invalid substituents.", Chem.MolFromSmiles(substituents[0]))

    # Prepare arguments for parallel processing
    tasks = []
    for r in range(1, len(attachment_points) + 1):
        for points in combinations(attachment_points, r):
            for subs in combinations(substituent_mols, r):

                # Clone the scaffold for modification
                mol = Chem.RWMol(scaffold)

                # Replace '*' with the corresponding substituents
                for idx, sub in zip(points, subs):
                    mol.ReplaceAtom(idx, sub.GetAtomWithIdx(0))

                # Add the task with the modified molecule
                tasks.append((mol, points, subs))

    # Use multiprocessing to parallelize
    results = []
    with Pool(cpu_count()) as pool:
        for result in pool.imap(substitute_and_evaluate, tasks):
            if result:
                results.append(result)

    return results
