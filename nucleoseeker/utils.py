import pandas as pd
import os
import subprocess
import logging

"""
This module contains utility functions for the nucleoseeker package.

The functions in this module are used to download PDB files, download fasta files, generate contact maps, extract sequences from fasta files, etc.
"""

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




    