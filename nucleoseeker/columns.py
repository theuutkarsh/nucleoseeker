"""
The columns that will be extracted from the PDB metadata.
We find these to be the most relevant as they capture the essence of structures in the PDB database.

Columns:
    rcsb_id: The unique identifier for the structure in the PDB database.
    
    exptl_method: The method used to determine the structure.

    release_date: The date the structure was released.

    polymer_entity_instance_count: The number of polymer entities in the structure.

    polymer_entity_count_RNA: The number of RNA polymer entities in the structure.

    resolution: The resolution of the structure.

    selected_polymer_entity_types: The types of polymer entities in the structure.
    
    pdbx_keywords: Keywords associated with the structure.

"""
COLUMNS = [
    'rcsb_id',
    'exptl_method',
    'release_date',
    'polymer_entity_instance_count',
    'polymer_entity_count_RNA',
    'resolution',
    'selected_polymer_entity_types',
    'pdbx_keywords'
]