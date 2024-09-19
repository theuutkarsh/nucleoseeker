import pandas as pd
import numpy as np
import glob
import os
import typing
import subprocess
import logging
import Bio.PDB
import Bio
import torch
import math

"""
This module contains utility functions for the nucleoseeker package.

The functions in this module are used to download PDB files, download fasta files, generate contact maps, extract sequences from fasta files, etc.
"""

def remove_backbone_contacts(contacts: np.ndarray, width: int = 0):
    if isinstance(contacts, torch.Tensor):
        contacts = contacts.numpy()
    for i in range(1, width + 1):
        np.fill_diagonal(contacts[i:], 0)
        np.fill_diagonal(contacts[:, i:], 0)
        np.fill_diagonal(contacts[:-i], 0)
        np.fill_diagonal(contacts[:, :-i], 0)

    # Make sure the diagonal is zero
    if width != 0:
        np.fill_diagonal(contacts, 0)

    return contacts

def contacts_from_pdb(
        structure: Bio.PDB.Structure.Structure,
        chain: str,
        distance_threshold: float,
        sequence_length: int,
    ):
    chunk_size = 600
        
    # assert len(structure) == 1
    # assert len(structure[0]) == 1

    # chain = structure[0].get_list()[0]
    chain = structure[0][chain]
    # print(structure.id, sequence_length)
    standard_residues = [residue for residue in chain]
    if len(standard_residues) < sequence_length:
        logging.warning(f"Sequence length {len(standard_residues)} is less than the specified sequence length {sequence_length}")
    standard_residues = standard_residues[:sequence_length]
    atom_coords = [[a.get_coord() for a in r] for r in standard_residues]
    nums_atoms = [len(atom_coords_r) for atom_coords_r in atom_coords]
    num_atoms_max = max(nums_atoms)
    num_residues = len(atom_coords)
    num_chunks = int(math.ceil(num_residues / chunk_size))

    distances = torch.zeros((num_residues, num_residues)) # ! RuntimeError: Cannot re-initialize CUDA in forked subprocess. To use CUDA with multiprocessing, you must use the 'spawn' start method 
    atom_coords_t = torch.full((num_residues, num_atoms_max, 3), torch.inf)  # [R, A, 3]
    for idx in range(num_residues):
        atom_coords_rt = torch.tensor(np.array(atom_coords[idx]))
        atom_coords_t[idx, :nums_atoms[idx], :] = atom_coords_rt

    for idx in range(num_chunks):
        residues_slice = slice(idx * chunk_size, (idx + 1) * chunk_size)

        atom_coords_t1 = atom_coords_t[residues_slice, :, :].view(-1, 1, num_atoms_max, 1, 3).expand(-1, num_residues, num_atoms_max, num_atoms_max, 3)  # [RC, R, A, A, 3]
        atom_coords_t2 = atom_coords_t.view(1, num_residues, 1, num_atoms_max, 3).expand(atom_coords_t1.shape[0], num_residues, num_atoms_max, num_atoms_max, 3)  # [RC, R, A, A, 3]

        distances_chunk = torch.linalg.vector_norm(atom_coords_t1 - atom_coords_t2, dim=-1)  # [RC, R, A, A]
        distances_chunk = torch.nan_to_num(distances_chunk, nan=torch.inf, posinf=torch.inf)
        distances_chunk = torch.amin(distances_chunk, dim=(-1, -2))  # [RC, R]
        distances[residues_slice, :] = distances_chunk
    contacts = distances < distance_threshold
    return contacts


def generate_contact_map_from_mmcif_file(
        mmcif_file: str,
        output_dir: str,
        chain: str,
        seq_len: int,
        distance_cutoff: float = 8.0,
        save: bool = True,
        width: int = 0
    ):
    """
    It generates a contact map from a mmcif file.
    It used MMCIFParser from Bio.PDB to parse the mmcif file.
    """
    from Bio.PDB import MMCIFParser

    name = os.path.basename(mmcif_file).split('.')[0]
    structure = MMCIFParser().get_structure(name, mmcif_file)

    contact_map = contacts_from_pdb(
        structure,
        chain=chain,
        distance_threshold=distance_cutoff,
        sequence_length=seq_len,
    )
    if save:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        np.save(f"{output_dir}/{name}_{chain}.npy", contact_map)
    if width != 0:
        contact_map = remove_backbone_contacts(contact_map, width)
    return contact_map

def generate_contact_map_from_pdb_file(
        pdb_file: str,
        output_dir: str,
        chain: str,
        seq_len: int,
        distance_cutoff: float = 8.0,
        save: bool = True,
        width: int = 0,
    ):
    """
    It generates a contact map from a pdb file.
    It uses the PDBParser from the Bio.PDB module to parse the pdb file.
    """
    from Bio.PDB import PDBParser

    name = os.path.basename(pdb_file).split('.')[0]
    structure = PDBParser().get_structure(name, pdb_file)
    contact_map = contacts_from_pdb(
        structure,
        chain=chain,
        distance_threshold=distance_cutoff,
        sequence_length=seq_len,
    )
    if save:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
        np.save(f"{output_dir}/{name}_{chain}.npy", contact_map)
    if width != 0:
        contact_map = remove_backbone_contacts(contact_map, width)

    return contact_map

def extract_sequence_from_combined_fasta(
        combined_fasta_file: str,
        output_dir: str,
        pdb_id: str,
    ):
    """
    It takes as input a combined fasta file and extracts the sequence of a given pdb_id.
    If the fasta header is not of the format {pdb_id}_{chain}, it will throw an error.
    Args:
        combined_fasta_file (str): The combined fasta file.
        output_dir (str): The output directory.
        pdb_id (str): The pdb_id.

    Returns:
        Saves the sequence to a file with name {output_dir}/{pdb_id}_{chain}.fa
    """
    with open(combined_fasta_file, 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            header = lines[i].strip()
            sequence = lines[i+1].strip()
            if header.startswith(f">{pdb_id}_"):
                chain = header.split('_')[1]
                with open(f"{output_dir}/{pdb_id}_{chain}.fa", 'w') as out:
                    out.write(f"{header}\n{sequence}\n")
                return
        raise ValueError(f"Could not find sequence for pdb_id {pdb_id} in {combined_fasta_file}")

def extract_sequences_for_pdb_ids(
        combined_fasta_file: str,
        output_dir: str,
        pdb_ids: list,
    ):
    """
    It takes as input a combined fasta file and extracts the sequence of a given pdb_id.
    If the fasta header is not of the format {pdb_id}_{chain}, it will throw an error.
    Args:
        combined_fasta_file (str): The combined fasta file.
        output_dir (str): The output directory.
        pdb_ids (list): The list of pdb_ids.

    Returns:
        Saves the sequence to a file with name {output_dir}/{pdb_id}_{chain}.fa
    """
    with open(combined_fasta_file, 'r') as f:
        lines = f.readlines()
        for i in range(0, len(lines), 2):
            header = lines[i].strip()
            sequence = lines[i+1].strip()
            for pdb_id in pdb_ids:
                if header.startswith(f">{pdb_id}_"):
                    chain = header.split('_')[1]
                    with open(f"{output_dir}/{pdb_id}_{chain}.fa", 'w') as out:
                        out.write(f"{header}\n{sequence}\n")


def get_final_fam_pdb_chain_csv(clean_tblout_path):
    fam_pdb_chain_csv_path = os.path.join(os.path.dirname(clean_tblout_path), 'fam_pdb_chain.csv')
    temp_file_path = os.path.join(os.path.dirname(clean_tblout_path), 'temp.csv')
    logging.debug('Getting final PDB list')
    cmd = f"awk 'NR>2{{print $2, $3, $16}}' {clean_tblout_path} > {temp_file_path}"
    subprocess.run(cmd, shell=True, check=True)
    df = pd.read_csv(temp_file_path, header=None, sep=' ')
    df.columns = ['fam_id', 'rcsb_id', 'e_value']
    min_indices = df.groupby(['rcsb_id'])['e_value'].idxmin()
    filtered_df = df.loc[min_indices]
    filtered_df['chain_id'] = filtered_df['rcsb_id'].str.split('_').str[1]
    filtered_df['rcsb_id'] = filtered_df['rcsb_id'].str.split('_').str[0]
    filtered_df = filtered_df[['fam_id', 'rcsb_id', 'chain_id']]
    filtered_df.to_csv(fam_pdb_chain_csv_path, index=False)
    os.remove(temp_file_path)
    logging.debug('Final PDB list obtained')

def download_pdb_file(pdb_id, output_dir):
    """
    This is used to download a PDB file.

    Args:
        pdb_id (str): The PDB ID.
        output_dir (str): The output directory.

    """
    os.system(f"wget -qc -nc -O {output_dir}/{pdb_id}.cif.gz 'https://files.rcsb.org/download/{pdb_id}.cif.gz'; gunzip {output_dir}/{pdb_id}.cif.gz")

def download_pdb_files(pdb_list: list, output_dir: str):
    """
    This is used to download a list of PDB files.

    Args:
        pdb_list (list): The list of PDB IDs.
        output_dir (str): The output directory.

    """
    for pdb_id in pdb_list:
        download_pdb_file(pdb_id, output_dir)

def download_fasta_file(pdb_id, output_dir):
    """
    This is used to download a fasta file.

    Args:
        pdb_id (str): The PDB ID.
        output_dir (str): The output directory.

    """
    os.system(f"wget -qc -nc -O {output_dir}{pdb_id}.fa 'https://www.rcsb.org/fasta/entry/{pdb_id}' ")

def download_fasta_files(pdb_list: list, output_dir: str):
    """
    This is used to download a list of fasta files.

    Args:
        pdb_list (list): The list of PDB IDs.
        output_dir (str): The output directory.

    """
    for pdb_id in pdb_list:
        download_fasta_file(pdb_id, output_dir)

def parse_fasta(fasta_file):
    sequences = {}
    header = ''
    with open(fasta_file, 'r') as file:
        for line in file:
            if line.startswith('>'):
                header = line.strip()[1:]
                sequences[header] = ''
            else:
                sequences[header] += line.strip()
    return sequences

def write_sequences_to_files(sequences, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    for header, sequence in sequences.items():
        with open(os.path.join(output_dir, header + '.fa'), 'w') as file:
            file.write('>' + header + '\n')
            file.write(sequence + '\n')

def save_sequences(final_fasta, output_dir):
    fasta_file = final_fasta  # Change this to your input FASTA file
    sequences = parse_fasta(fasta_file)
    write_sequences_to_files(sequences, output_dir)




    