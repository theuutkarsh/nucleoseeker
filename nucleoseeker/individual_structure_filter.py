from itertools import chain
import os
import pathlib
import logging
import nucleoseeker.utils as utils
from Bio.PDB.MMCIFParser import MMCIF2Dict

DATA_PATH = os.environ.get('DATA_PATH')

class IndividualStructureFilter:
    """
    This is used to apply filters on PDB files. It analyses the PDB files at the chain and polymer level.
    It checks if each chain in the PDB file is of the given polymer type and has the required sequence length.
    
    Args:
        pdb_parser (Bio.PDB.MMCIFParser): PDB parser.
        pdb_id (str): PDB ID.
        polymer_type (str): Polymer type, common possible values are 'polyribonucleotide', 'polydeoxyribonucleotide', 'polypeptide' etc.
        sequence_length (int): Sequence length.
        auto_download (bool): Auto download PDB file if not found.
        
    Attributes:
        pdb_id (str): PDB ID.
        polymer_type (str): Polymer type.
        sequence_length (int): Sequence length.
        auto_download (bool): Auto download PDB file.
        pdb_parser (Bio.PDB.MMCIFParser): PDB parser.
        pdb_file (str): Path to the PDB file.
        structure_dict (dict): Structure dictionary.
        
    Methods:
        _check_and_auto_download: Check if the PDB file exists and auto download if required.
        check_polymer_type: Check if the PDB file satisfies the given criteria.
        
    """
    def __init__(
            self,
            pdb_parser: MMCIF2Dict,
            pdb_id: str,
            polymer_type: str,
            sequence_length: int,
            auto_download: bool = False,
        ) -> None:
        self.pdb_id = pdb_id
        self.polymer_type = polymer_type
        assert sequence_length > 0, 'Sequence length should be greater than 0.'
        self.sequence_length = sequence_length
        self.auto_download = auto_download
        self.pdb_parser = pdb_parser
        self.pdb_file = self._check_and_auto_download()
        self.structure_dict = self.pdb_parser(self.pdb_file)

    def _check_and_auto_download(self):
        """
        Check if the PDB file exists and auto download if required.
        
        Returns:
            str: Path to the PDB file."""
        pdb_file = f'{self.pdb_id}.cif'
        DATA_PATH_PDB = DATA_PATH + '/pdb_files'
        path = os.path.join(DATA_PATH_PDB, pdb_file)
        if not os.path.exists(path):
            os.makedirs(DATA_PATH_PDB, exist_ok=True)
            if self.auto_download:
                logging.info(f'PDB file {pdb_file} not found. Downloading...')
                utils.download_pdb_file(self.pdb_id, DATA_PATH_PDB)
                logging.info(f'PDB file saved at {DATA_PATH_PDB}')
            else:
                raise FileNotFoundError(f'File {pdb_file} not found and auto_download is set to False')
        return path
    
    def check_polymer_type(self):
        """
        Check if the PDB file contains the given polymer type and has the required sequence length.

        Returns:
            list: List of tuples containing PDB ID, chain ID and corresponding sequence.
        """
        polymer_types = self.structure_dict.get('_entity_poly.type', [])
        chain_ids = self.structure_dict.get('_entity_poly.pdbx_strand_id', [])
        data_list = []
        if self.polymer_type is None:
            for polymer_type, polymer_type_id in zip(polymer_types, chain_ids):
                chain_id = polymer_type_id.split(',')[0]
                chain_number = chain_ids.index(polymer_type_id)
                sequence_ = self.structure_dict.get('_entity_poly.pdbx_seq_one_letter_code_can', [])[chain_number]
                if sequence_ and len(sequence_) >= self.sequence_length:
                    sequence = sequence_.replace('\n', '')
                    if self.all_characters_are_n(sequence) or self.all_characters_are_x(sequence):
                        logging.info(f"The chain {chain_id} in {self.pdb_id} has all 'N' characters. Skipping...")
                        continue
                    data_list.append((self.pdb_id, chain_id, sequence))
        else:
            for polymer_type, polymer_type_id in zip(polymer_types, chain_ids):
                if polymer_type == self.polymer_type:
                    chain_id = polymer_type_id.split(',')[0]
                    chain_number = chain_ids.index(polymer_type_id)
                    sequence_ = self.structure_dict.get('_entity_poly.pdbx_seq_one_letter_code_can', [])[chain_number]
                    if sequence_ and len(sequence_) >= self.sequence_length:
                        sequence = sequence_.replace('\n', '')
                        if self.all_characters_are_n(sequence) or self.all_characters_are_x(sequence):
                            logging.info(f"The chain {chain_id} in {self.pdb_id} has all 'N' characters. Skipping...")
                            continue
                        data_list.append((self.pdb_id, chain_id, sequence))

        if not data_list:
            logging.info(f'The {self.pdb_id} does not satisfy the given criteria.')
        return data_list
    
    def all_characters_are_n(self, s):
        """
        Check if all characters in the string are 'N'.
        This is helpful to remove chains with all 'N' characters as they can't be processed by Clustal Omega.
        
        Args:
            s (str): Input string.
            
        Returns:
                bool: True if all characters are 'N', False otherwise.
        """
        return all(char == 'N' for char in s)
    
    def all_characters_are_x(self, s):
        """
        Check if all characters in the string are 'X'.
        This is helpful to remove chains with all 'X' characters as they can't be processed by Clustal Omega.
        
        Args:
            s (str): Input string.
            
        Returns:
                bool: True if all characters are 'X', False otherwise.
        """
        return all(char == 'X' for char in s)
