import os
import pathlib
import logging
import subprocess
import re
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIF2Dict
from nucleoseeker.individual_structure_filter import IndividualStructureFilter


class StructureComparisonFilter:
    """
    This is used to apply PDBFilter for multiple PDB files.
    It analyses the PDB files at the chain and polymer level.
    It also calculates the sequence identity between the sequences of the PDB files.
    You can use two alignment tools: clustal and emboss. They should be installed on your system.
    We recommend using clustal as it is faster.

    This gives the final dataframe after applying all the filters.

    Args:
        polymer_type (str): Polymer type.
        sequence_length (int): Sequence length.
        sequence_identity (float): Sequence identity.
        auto_download (bool): Auto download PDB file.
        alignment_tool (str): Alignment tool. Default is clustal.

    Attributes:
        polymer_type (str): Polymer type.
        sequence_length (int): Sequence length.
        sequence_identity (float): Sequence identity.
        auto_download (bool): Auto download PDB file.
        alignment_tool (str): Alignment tool.
        pdb_parser (Bio.PDB.MMCIFParser): PDB parser.
        DATA_PATH (str): Data path.
        sequence_identity_mat_file (str): Sequence identity matrix file.
        sequence_identity_mat_path (str): Sequence identity matrix path.

    Methods:
        apply_filters_on_pdb_id: Apply filters on PDB ID.
        apply_filters_on_list: Apply filters on list of PDB IDs.
        create_combined_fasta_file: Create combined fasta file.
        create_sequence_identity_mat_emboss: Create sequence identity matrix using emboss.
        get_sequence_identity_df_emboss: Get sequence identity dataframe using emboss.
        create_sequence_identity_mat_clustal: Create sequence identity matrix using clustal.
        get_sequence_identity_df_clustal: Get sequence identity dataframe using clustal.
        apply_filter_on_df: Apply filter on dataframe.
        create_final_fasta_file: Create final fasta file.

    #! NOTE: When applying similarity cutoff, if the df contains resolution column with no values, it will raise an error.

    """
    def __init__(
            self,
            polymer_type: str,
            sequence_length: int,
            sequence_identity: float,
            auto_download: bool = False,
            alignment_tool: str = 'clustal',
        ) -> None:
        self.polymer_type = polymer_type
        assert sequence_length > 0, 'Sequence length should be greater than 0.'
        assert 0 <= sequence_identity <= 100, 'Sequence identity should be between 0 and 100.'
        self.sequence_length = sequence_length
        self.sequence_identity = sequence_identity
        self.auto_download = auto_download
        assert alignment_tool in ['clustal', 'emboss'], 'Alignment tool should be either clustal or emboss.'
        self.alignment_tool = alignment_tool
        
        self.pdb_parser = MMCIF2Dict
        
        self.DATA_PATH = os.environ.get('DATA_PATH')
        self.sequence_identity_mat_file = f'sequence_identity_mat_{self.alignment_tool}.csv'
        self.sequence_identity_mat_path = os.path.join(self.DATA_PATH, self.sequence_identity_mat_file)
        
        logging.basicConfig(level=logging.INFO)

    def apply_filters_on_pdb_id(self, pdb_id: str):
        """
        Using this method, you can apply filters on a single PDB ID.

        Args:
            pdb_id (str): PDB ID.

        Returns:
            list: List of tuples containing PDB ID, chain ID, and sequence.
        """
        pdb_filters = IndividualStructureFilter(
            pdb_parser=self.pdb_parser,
            pdb_id=pdb_id,
            polymer_type=self.polymer_type,
            sequence_length=self.sequence_length,
            auto_download=self.auto_download
        )
        single_pdb_data_list = pdb_filters.check_polymer_type()
        return single_pdb_data_list
    
    def apply_filters_on_list(self, pdb_list: list):
        """
        Using this method, you can apply filters on a list of PDB IDs.

        Args:
            pdb_list (list): List of PDB IDs.

        Returns:
            list: List of tuples containing PDB ID, chain ID, and sequence.
        """
        multiple_pdbs_data_list = []
        for pdb_id in pdb_list:
            data = self.apply_filters_on_pdb_id(pdb_id)
            if data:
                multiple_pdbs_data_list.append(data[0]) #append the first element of list of tuples

        if not multiple_pdbs_data_list:
            raise ValueError("No PDB files satisfy the given criteria. Please check the arguments for 'polymer_type' and 'sequence_length'.")
        return multiple_pdbs_data_list
    
    def create_combined_fasta_file(self, data_list: list):
        """
        Create a combined fasta file from the data list.
        
        Args:
            data_list (list): List of tuples containing PDB ID, chain ID, and sequence.
            
        Returns:
            None"""
        self.combined_fasta_file = os.path.join(self.DATA_PATH, f'combined.fasta')
        with open(self.combined_fasta_file, 'w') as f:
            for data in data_list:
                pdb_id, chain_id, sequence = data
                f.write(f'>{pdb_id}_{chain_id}\n{sequence}\n')

    def create_sequence_identity_mat_emboss(self, data_list: list):
        """
        Create sequence identity matrix using emboss.
        
        Args:
            data_list (list): List of tuples containing PDB ID, chain ID, and sequence.
            
        Returns:
            pd.DataFrame: Sequence identity matrix.
        """
        self.DATA_PATH_FASTA = pathlib.Path(__file__).resolve().parents[1]/ 'fasta_log'
        os.makedirs(self.DATA_PATH_FASTA, exist_ok=True)
        file1 = os.path.join(self.DATA_PATH_FASTA, 'file1.fasta')
        file2 = os.path.join(self.DATA_PATH_FASTA, 'file2.fasta')

        for i, (pdb_id, chain_id, sequence) in enumerate(data_list):
            with open(file1, 'w') as f:
                f.write(f'>{pdb_id}_{chain_id}\n{sequence}\n')

            for j, (pdb_id2, chain_id2, sequence2) in enumerate(data_list):
                with open(file2, 'w') as f:
                    f.write(f'>{pdb_id2}_{chain_id2}\n{sequence2}\n')

                command = f"needle \
                    -auto\
                    -stdout\
                    -asequence {file1}\
                    -bsequence {file2}\
                    -datafile EDNAFULL\
                    -gapopen 10.0\
                    -gapextend 0.5\
                    -endopen 10.0\
                    -endextend 0.5\
                    -aformat3 pair\
                    -snucleotide1\
                    -snucleotide2\
                    "
                try:
                    result = subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE)
                    result = result.stdout.decode('utf-8')
                    pattern = r'Identity:\s+(\d+/\d+)\s+\((\d+\.\d+)%\)'
                    match = re.search(pattern, result)
                    if match is None:
                        identity_percent = 0.0
                    else:
                        identity_percent = match.group(2)
                    sequence_identity_mat.iloc[i, j] = float(identity_percent)
                except subprocess.CalledProcessError as e:
                    logging.error(f'Error occurred: {e}')
                    raise e
        sequence_identity_mat = sequence_identity_mat.reset_index()
        sequence_identity_mat.rename(columns={'index': 'rcsb_id'}, inplace=True)
        
        sequence_identity_mat.to_csv(self.sequence_identity_mat_path, index=False)
        return sequence_identity_mat

    def get_sequence_identity_df_emboss(self, data_list: list):
        """
        Get sequence identity dataframe using emboss.
        
        Args:
            data_list (list): List of tuples containing PDB ID, chain ID, and sequence.
            
        Returns:
            pd.DataFrame: Sequence identity dataframe.
        """
        if not os.path.exists(self.sequence_identity_mat_path):
            self.create_sequence_identity_mat_emboss(data_list)
        df = pd.read_csv(self.sequence_identity_mat_path)
        return df

    def create_sequence_identity_mat_clustal(self, combined_fasta_file: str):
        """
        Create sequence identity matrix using clustal.
        
        Args:
            combined_fasta_file (str): Combined fasta file.
            
        Returns:
            None
        """
        command = f"clustalo \
            -i {combined_fasta_file} \
            --distmat-out={self.sequence_identity_mat_path}\
            --percent-id \
            --full \
            --force\
        "
        try:
            subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE)
        except subprocess.CalledProcessError as e:
            logging.error(f'Error occurred: {e}')
            raise e
        
    def get_sequence_identity_df_clustal(self, data_list: list):
        """
        Get sequence identity dataframe using clustal.
        
        Args:
            data_list (list): List of tuples containing PDB ID, chain ID, and sequence.
            
        Returns:
            pd.DataFrame: Sequence identity dataframe.
            
        """
        if not os.path.exists(self.sequence_identity_mat_path):
            self.create_combined_fasta_file(data_list)
            self.create_sequence_identity_mat_clustal(self.combined_fasta_file)
        df = pd.read_csv(self.sequence_identity_mat_path, sep=' ', header=None, skiprows=1, skipinitialspace=True)
        df.columns = ['rcsb_id'] +  df[0].to_list() 
        df = df.round(1)
        return df
    
    def apply_filter_on_df(self, df: pd.DataFrame, data_list: list):
        """
        Apply filter on dataframe.
        
        Args:
            df (pd.DataFrame): Dataframe.
            data_list (list): List of tuples containing PDB ID, chain ID, and sequence.
            
        Returns:
            pd.DataFrame: Filtered dataframe.
        """
        if 'resolution' not in df.columns:
            raise ValueError('The dataframe does not have the resolution column. Please use the StructureLevelFilter class to get the right dataframe')
        rcsb_id_with_no_resolution = df[df['resolution'].isnull()]['rcsb_id'].to_list()
        #drop data from data list based on rcsb_id which don't have resolution
        data_list = [data for data in data_list if data[0] not in rcsb_id_with_no_resolution]
        df.dropna(subset=['resolution'], inplace=True)
        df = df[df['rcsb_id'].isin([data[0] for data in data_list])].copy()
        df.reset_index(drop=True, inplace=True)
        if self.alignment_tool == 'clustal':
            similarity_df = self.get_sequence_identity_df_clustal(data_list)
        else:
            similarity_df = self.get_sequence_identity_df_emboss(data_list)
        similarity_df = similarity_df[similarity_df['rcsb_id'].str.split('_').str[0].isin(df['rcsb_id'])]
        similarity_df.dropna(inplace=True)
        similarity_df.reset_index(drop=True, inplace=True)
        df['chain_id'] = similarity_df['rcsb_id'].str.split('_').str[1]
        df['sequence'] = [data[2] for data in data_list]
        assert df['rcsb_id'].to_list() == similarity_df['rcsb_id'].str.split('_').str[0].to_list(), 'The PDB ids in the dataframes do not match'
        assert df.shape[0] == similarity_df.shape[0], f'The number of rows in the dataframes do not match. They are {df.shape[0]} and {similarity_df.shape[0]}'

        similarity_df['resolution'] = df['resolution'].astype(float)
        dissimilar_pdb_lowest_resolutions = []       
        logging.info('Applying similarity cutoff, to get the final dataframe...')
        name_list = ['_'.join([data[0], data[1]]) for data in data_list]
        for name in name_list:
            si = similarity_df[similarity_df[name] > self.sequence_identity][['rcsb_id', name, 'resolution']]
            si.reset_index(drop = True, inplace = True)
            si = si.iloc[si['resolution'].idxmin()]
            dissimilar_pdb_lowest_resolutions.append(si['rcsb_id'])
        pdb_names = [pdb_name.split('_')[0] for pdb_name in dissimilar_pdb_lowest_resolutions]
        df_final = df[df['rcsb_id'].isin(pdb_names)]
        return df_final
    
    def create_final_fasta_file(self, df_final: pd.DataFrame):
        """
        Create final fasta file.
        
        Args:
            df_final (pd.DataFrame): Final dataframe.
            
        Returns:
            None
        """
        final_fasta_file = os.path.join(self.DATA_PATH, f'final.fasta')
        with open(final_fasta_file, 'w') as f:
            for index, row in df_final.iterrows():
                f.write(f'>{row["rcsb_id"]}_{row["chain_id"]}\n{row["sequence"]}\n')
