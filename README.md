"""
Edge Partition Based on Vertex Degrees (Chemical Graph Theory)

This script:
1. Takes a SMILES string
2. Builds molecular graph using RDKit
3. Computes vertex degrees
4. Partitions edges based on degree pairs (du, dv)

Run:
    python edge_partition.py

Change SMILES inside the script or pass via input().
"""

from rdkit import Chem
from collections import defaultdict


def edge_partition_from_smiles(smiles: str):
    """
    Compute edge partition of a molecular graph based on vertex degrees.

    Parameters:
        smiles (str): SMILES string of the molecule

    Returns:
        dict: {(du, dv): [(u, v), ...], ...}
    """
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        print("❌ Invalid SMILES string! Please check your input.")
        return {}

    # Extract edges (bonds)
    edges = []
    for bond in mol.GetBonds():
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()
        edges.append((u, v))

    # Compute degree of each vertex (atom)
    degree = {atom.GetIdx(): atom.GetDegree() for atom in mol.GetAtoms()}

    # Partition edges based on (du, dv)
    partition = defaultdict(list)
    for (u, v) in edges:
        du, dv = degree[u], degree[v]
        key = tuple(sorted((du, dv)))
        partition[key].append((u, v))

    return partition


def main():
    print("=== Edge Partition from SMILES ===\n")

    # 🔹 Default SMILES (change if needed)
    smiles_input = "CC(C)N1C(=CC=N1)C2=C(C=CC=N2)CO"

    # Uncomment below to take user input
    # smiles_input = input("Enter SMILES string: ")

    result = edge_partition_from_smiles(smiles_input)

    if not result:
        return

    print("\nEdge Partition based on degrees:\n")
    for key, edge_list in result.items():
        print(f"Degree Pair {key} : {len(edge_list)} edges -> {edge_list}")


if __name__ == "__main__":
    main()
