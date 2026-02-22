

[Properties.xlsx](https://github.com/user-attachments/files/25464994/Properties.xlsx)
This script contains the physicochemical properties used in the associated research paper.

Edge Partition Based on Vertex Degrees (Chemical Graph Theory)

edge_partition.py this script:
1. Takes a SMILES string
2. Builds molecular graph using RDKit
3. Computes vertex degrees
4. Partitions edges based on degree pairs (du, dv)


   

    print("\nEdge Partition based on degrees:\n")
    for key, edge_list in result.items():
        print(f"Degree Pair {key} : {len(edge_list)} edges -> {edge_list}")


if __name__ == "__main__":
    main()
