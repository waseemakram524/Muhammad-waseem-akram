from rdkit import Chem
from collections import defaultdict

def edge_partition_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        print("❌ Invalid SMILES string! Check your input.")
        return {}
    
    edges = []
    for bond in mol.GetBonds():
        u = bond.GetBeginAtomIdx()
        v = bond.GetEndAtomIdx()
        edges.append((u, v))
    
    degree = {}
    for atom in mol.GetAtoms():
        degree[atom.GetIdx()] = atom.GetDegree()
    
    partition = defaultdict(list)
    
    for (u, v) in edges:
        du, dv = degree[u], degree[v]
        key = tuple(sorted((du, dv)))
        partition[key].append((u, v))
    
    return partition


# 🔹 PUT SMILES 
smiles_input = "CC(C)N1C(=CC=N1)C2=C(C=CC=N2)COC3=CC=CC(=C3C=O)O"   # <-- yahan change karein

result = edge_partition_from_smiles(smiles_input)

print("\nEdge Partition based on degrees:\n")
for key, edge_list in result.items():
    print(f"Degree Pair {key} : {len(edge_list)} edges -> {edge_list}")
