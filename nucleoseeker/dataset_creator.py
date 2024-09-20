import os
import subprocess
import logging
import typing
import platform
import argparse
from typing import Optional, Union
import pandas as pd
from nucleoseeker.dataset_download import DatasetDownload
from nucleoseeker.metadata_filter import MetadataFilter
from nucleoseeker.structure_comparison_filter import StructureComparisonFilter
import nucleoseeker.utils as utils
from nucleoseeker.check_tool import check_tool



DATA_PATH = os.environ.get('DATA_PATH')

class DatasetCreator:
    """
    It is used to create a dataset of RNA structure from the RCSB PDB database.
    This makes use of the DatasetDownload, StructureLevelFilter and PolymerLevelFilter classes.
    It can be used to download the dataset, apply filters at the structure and polymer level,
    all the parameters can be controlled by the user.
    Make sure to download Rfam.cm file from https://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz, unzip it and use **infernal** **cmpress** to compress it.
    Desired DATA_PATH can be set as an environment variable.
    When you run this class for the first time, it will create a folder for PDB files or you can create a dir called 'pdb_files' one level below the current directory.

    Args:
        dataset_name (str): Name of the dataset.
        rfam_cm_path (str): Path to the Rfam.cm file.
        structure_determination_methodology (str): Structure determination methodology.
        rcsb_entity_polymer_type (str): Polymer type.
        dstart (int): Start index for pagination.
        dend (int): End index for pagination.
        download_all (bool): Download all entries.
        exptl_method (list): Experimental method.
        resolution (float): Resolution.
        year_range (int or list): Year range.
        polymer_entity_instance_count (int): Polymer entity instance count.
        polymer_entity_count_RNA (int): Polymer entity count RNA.
        selected_polymer_entity_types (list): Selected polymer entity types.
        pdbx_keywords (list): PDBx keywords.
        polymer_type (str): Polymer type.
        sequence_length (int): Sequence length.
        sequence_identity (float): Sequence identity.
        auto_download (bool): Auto download PDB file.
        alignment_tool (str): Alignment tool.
        e_value_cmscan (float): E-value for CMScan.
        save (bool): Save filtered data, if set to False, it will not save data at each filter level.
    
    Attributes:
        dataset_name (str): Name of the dataset.
        structure_determination_methodology (str): Structure determination methodology.
        rcsb_entity_polymer_type (str): Polymer type.
            ['Protein', 'DNA', 'RNA', 'NA-hybrid', 'Other'] a list can't be passed
        dstart (int): Start index for pagination.
        dend (int): End index for pagination.
        download_all (bool): Download all entries.
        exptl_method (list): Experimental method.
            ['X-RAY DIFFRACTION', 'SOLUTION NMR', 'ELECTRON MICROSCOPY', 'SOLID-STATE NMR'] a list can be passed
            Further details can be found at www.wwpdb.org or www.rcsb.org
        resolution (float): Resolution.
        year_range (int or list): Year range.
        polymer_entity_instance_count (int): Polymer entity instance count.
        polymer_entity_count_RNA (int): Polymer entity count RNA.
        selected_polymer_entity_types (list): Selected polymer entity types.
        pdbx_keywords (list): PDBx keywords.
            ['RNA', 'DNA/RNA', 'RIBOSOME', 'RIBOZYME'] a list can be passed
            Further details can be found at https://www.wwpdb.org or https://www.rcsb.org or in the raw data
        polymer_type (str): Polymer type.
        sequence_length (int): Sequence length.
        sequence_identity (float): Sequence identity.
        auto_download (bool): Auto download PDB file.
        alignment_tool (str): Alignment tool.
        rfam_cm_path (str): Path to the Rfam.cm file.
        e_value_cmscan (float): E-value for CMScan.
        save (bool): Save filtered data.
        dataset_files (str): Path to the dataset folder.
        final_fasta_path (str): Path to the final FASTA file.
        cmscan_path (str): Path to the CMScan output file.
        tblout (str): Path to the tblout file.
        clean_tblout_path (str): Path to the clean tblout file.
        final_fams_path (str): Path to the final families file.
        final_pdb_list_path (str): Path to the final PDB list file.
        final_chain_ids_path (str): Path to the final chain IDs file.
        cmd (str): cmscan command.
        data (DatasetDownload): DatasetDownload component.
        structure_filter (StructureLevelFilter): StructureLevelFilter component.
        polymer_level_filter (PolymerLevelFilter): PolymerLevelFilter component.
    
    Methods:
        save_filtered_data: Save PDB list and DataFrame for each level of filters.
        save_pdb_list: Save list of PDBs to a file.
        save_dataframe: Save DataFrame to a file.
        get_without_filter_df: Get dataframe without applying any filters.
        get_structure_filtered_df: Get dataframe after applying structure-level filters.
        get_polymer_filtered_df: Get dataframe after applying polymer-level filters.
        get_final_df: Get final dataframe after applying all filters.
        create_final_fasta_file: Create final FASTA file.
        run_cmscan: Run CMScan.
        get_final_families: Get final families.
        get_final_pdb_list: Get final PDB list.

    NOTE: When you use different params it is better to create a explantory dataset_name
    because internally all files will have same names and it will be hard to distinguish

    Example:
            dataset_name = 'test'
            os.environ['DATA_PATH'] = str(Path(__file__).resolve().parents[1] / f'data/{dataset_name}')
            DATA_PATH = os.environ.get('DATA_PATH')
            dc = DatasetCreator(**your_params_here**)
    """
    def __init__(
            self,
            dataset_name: str,
            rfam_cm_path: str,
            structure_determination_methodology: typing.Optional[str] = 'experimental',
            rcsb_entity_polymer_type: typing.Optional[str] = 'RNA',
            dstart: typing.Optional[int] = 0,
            dend: typing.Optional[int] = 10000,
            download_all: typing.Optional[bool] = False,
            exptl_method: typing.Optional[typing.List[str]] = ['X-RAY DIFFRACTION'],
            resolution: typing.Optional[float] = 3.6,
            year_range: typing.Optional[Union[int, typing.List[int]]] = 2024,
            polymer_entity_instance_count: typing.Optional[int] = None,
            polymer_entity_count_RNA: typing.Optional[int] = None,
            selected_polymer_entity_types: typing.Optional[typing.List[str]] = None,
            pdbx_keywords: typing.Optional[typing.List[str]] = ['RNA'],
            polymer_type: typing.Optional[str] = 'polyribonucleotide',
            sequence_length: typing.Optional[int] = 40,
            sequence_identity: typing.Optional[float] = 50.0,
            auto_download: typing.Optional[bool] = True,
            alignment_tool: typing.Optional[str] = 'clustal',
            e_value_cmscan: typing.Optional[float] = 0.0001,
            save: typing.Optional[bool] = False
        ) -> None:
        self.dataset_name = dataset_name
        self.structure_determination_methodology = structure_determination_methodology
        self.rcsb_entity_polymer_type = rcsb_entity_polymer_type
        self.dstart = dstart
        self.dend = dend
        self.download_all = download_all
        self.exptl_method = exptl_method
        self.resolution = resolution
        self.year_range = year_range
        self.polymer_entity_instance_count = polymer_entity_instance_count
        self.polymer_entity_count_RNA = polymer_entity_count_RNA
        self.selected_polymer_entity_types = selected_polymer_entity_types
        self.pdbx_keywords = pdbx_keywords
        self.polymer_type = polymer_type
        self.sequence_length = sequence_length
        self.sequence_identity = sequence_identity
        self.auto_download = auto_download
        self.alignment_tool = alignment_tool
        self.rfam_cm_path = rfam_cm_path
        assert e_value_cmscan > 0, 'E-value should be greater than 0.'
        self.e_value_cmscan = e_value_cmscan
        self.save = save

        # Ensure data directory exists
        os.makedirs(DATA_PATH, exist_ok=True)
        
        logging.debug(f'Data directory created: {DATA_PATH}')
        
        # Create a folder for this dataset
        self.dataset_files = os.path.join(DATA_PATH, f"{self.dataset_name}/files")
        os.makedirs(self.dataset_files, exist_ok=True)
        logging.debug(f'Dataset folder created: {self.dataset_files}')

        # File paths
        self.final_fasta_path = os.path.join(DATA_PATH, f"{self.dataset_name}/final.fasta")
        self.cmscan_path = os.path.join(DATA_PATH, f"{self.dataset_name}/cmscan.out")
        self.tblout = os.path.join(DATA_PATH, f"{self.dataset_name}/tblout.tblout")
        self.clean_tblout_path = os.path.join(DATA_PATH, f"{self.dataset_name}/clean_tblout.tblout")
        self.final_fams_path = os.path.join(DATA_PATH, f"{self.dataset_name}/final_families.txt")
        self.final_pdb_list_path = os.path.join(DATA_PATH, f"{self.dataset_name}/final_pdb_list.txt")
        self.final_chain_ids_path = os.path.join(DATA_PATH, f"{self.dataset_name}/final_chain_ids.txt")
        self.sequences_path = os.path.join(DATA_PATH, f"{self.dataset_name}/sequences")
        os.makedirs(self.sequences_path, exist_ok=True)

        self.use_cmscan = True
        if self.rcsb_entity_polymer_type == 'Protein':
            self.use_cmscan = False
        #cmscan commands
        standard_cmd = (
            f'cmscan --tblout {self.tblout} --notextw --cut_ga --FZ 5 --nohmmonly --cpu 4 '
            f'{self.rfam_cm_path} {self.final_fasta_path} > {self.cmscan_path}'
        )
        e_value_cmd = (
            f'cmscan --tblout {self.tblout} --notextw --FZ 5 --nohmmonly --cpu 4 '
            f'-E {self.e_value_cmscan} {self.rfam_cm_path} {self.final_fasta_path} > {self.cmscan_path}'
        )

        if self.e_value_cmscan is not None:
            self.cmd = e_value_cmd
        else:
            self.cmd = standard_cmd
        

        # Initialize components
        self.data = DatasetDownload(
            structure_determination_methodology=self.structure_determination_methodology,
            rcsb_entity_polymer_type=self.rcsb_entity_polymer_type,
            dstart=self.dstart,
            dend=self.dend,
            download_all=self.download_all
        )

        logging.debug('DatasetDownload component initialized')

        self.structure_filter = MetadataFilter(
            exptl_method=self.exptl_method,
            resolution=self.resolution,
            year_range=self.year_range,
            polymer_entity_instance_count=self.polymer_entity_instance_count,
            polymer_entity_count_RNA=self.polymer_entity_count_RNA,
            selected_polymer_entity_types=self.selected_polymer_entity_types,
            pdbx_keywords=self.pdbx_keywords
        )

        logging.debug('StructureLevelFilter component initialized')

        self.polymer_level_filter = StructureComparisonFilter(
            polymer_type=self.polymer_type,
            sequence_length=self.sequence_length,
            sequence_identity=self.sequence_identity,
            auto_download=self.auto_download,
            alignment_tool=self.alignment_tool
        )

        logging.debug('PolymerLevelFilter component initialized')

    def save_pdb_list(self, df: pd.DataFrame, filename: str) -> None:
        """Save list of PDBs to a file."""
        pdb_list = df['rcsb_id'].tolist()
        with open(os.path.join(self.dataset_files, filename), 'w') as f:
            f.write('\n'.join(pdb_list))
        logging.debug(f'PDB list saved to {filename}')

    def save_dataframe(self, df: pd.DataFrame, filename: str) -> None:
        """Save DataFrame to a file."""
        df.to_csv(os.path.join(self.dataset_files, filename), index=False)
        logging.debug(f'DataFrame saved to {filename}')
    
    def create_final_fasta_file(self, df) -> None:
        """Create final FASTA file."""
        return self.polymer_level_filter.create_final_fasta_file(df)

    def apply_structure_filters(self, df: pd.DataFrame) -> pd.DataFrame:
        """Apply structure-level filters."""
        df_structure_filtered = self.structure_filter.apply_filters(df)
        if self.save:
            self.save_pdb_list(df_structure_filtered, 'structure_filtered_pdb_list.txt')
            self.save_dataframe(df_structure_filtered, 'structure_filtered_dataframe.csv')
        return df_structure_filtered
    
    def apply_polymer_filters(self, pdb_list: list) -> typing.Tuple[pd.DataFrame, list]:
        """Apply polymer-level filters."""
        df_polymer_filtered_list = self.polymer_level_filter.apply_filters_on_list(pdb_list)
        df_polymer_filtered_df = pd.DataFrame(df_polymer_filtered_list, columns=['rcsb_id', 'chain_id', 'sequence'])
        if self.save:
            self.save_dataframe(df_polymer_filtered_df, 'polymer_filtered_dataframe.csv')
            self.save_pdb_list(df_polymer_filtered_df, 'polymer_filtered_pdb_list.txt')
        return df_polymer_filtered_df, df_polymer_filtered_list
    
    def apply_filters(self, df: pd.DataFrame, df_polymer_filtered_list: list) -> pd.DataFrame:
        """Apply all filters."""
        df_final = self.polymer_level_filter.apply_filter_on_df(df, df_polymer_filtered_list)
        if self.save:
            self.save_pdb_list(df_final, 'final_pdb_list.txt')
            self.save_dataframe(df_final, 'final_dataframe.csv')
        return df_final
        
    
    def run(self) -> None:
        """Run CMScan."""

        df = self.data.df

        df_structure_filtered = self.apply_structure_filters(df)
        pdb_id_list = df_structure_filtered['rcsb_id'].to_list()
        df_polymer_filtered_list = self.apply_polymer_filters(pdb_id_list)[1]

        df_final = self.apply_filters(df, df_polymer_filtered_list)


        self.create_final_fasta_file(df_final)

        if self.use_cmscan:
            logging.debug('Running CMScan')
            try:
                subprocess.run(self.cmd, shell=True, check=True, stdout=subprocess.PIPE)
            except subprocess.CalledProcessError as e:
                logging.error(f'Error occurred while running cmscan: {e}')
                raise e
            # clean the tblout file
            if platform.system() == "Darwin":
                os.system(
                    f"tail -r {self.tblout} | awk '!/^#/{{p=1}} p' | tail -r > {self.clean_tblout_path}"
                )
            else:
                os.system(
                    f"tac {self.tblout} | awk '!/^#/{{p=1}} p' | tac > {self.clean_tblout_path}"
                )
            logging.debug('cmscan completed successfully')
        else:
            logging.warning('Polymer type must be RNA, if you wish to generate a dataset for proteins, please hmmer instead of infernal. Look here for more details: http://eddylab.org/software/hmmer/Userguide.pdf')


def main():
    logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s - %(message)s')
    check_tool('cmscan')
    check_tool('clustalo')

    parser = argparse.ArgumentParser(description='Create a dataset of RNA structure from the RCSB PDB database.')
    parser.add_argument('--dataset_name', type=str, help='Name of the dataset', required=True)
    parser.add_argument('--rfam_cm_path', type=str, help='Path to the Rfam.cm file', required=True)
    parser.add_argument('--structure_determination_methodology', type=str, help='Structure determination methodology', default='experimental')
    parser.add_argument('--rcsb_entity_polymer_type', type=str, help='Polymer type', default='RNA')
    parser.add_argument('--dstart', type=int, help='Start index for pagination', default=0)
    parser.add_argument('--dend', type=int, help='End index for pagination', default=10000)
    parser.add_argument('--download_all', type=bool, help='Download all entries', default=False)
    parser.add_argument('--exptl_method', nargs='*', help='Experimental method', default=None) #['X-RAY DIFFRACTION', 'SOLUTION NMR', 'ELECTRON MICROSCOPY', 'SOLID-STATE NMR']
    parser.add_argument('--resolution', type=float, help='Resolution', default=None) #3.6
    parser.add_argument('--year_range', type=int, nargs='*', help='Year range', default=None) #2024 or [2024, 2025]
    parser.add_argument('--polymer_entity_instance_count', type=int, help='Polymer entity instance count', default=None)
    parser.add_argument('--polymer_entity_count_RNA', type=int, help='Polymer entity count RNA', default=None)
    parser.add_argument('--selected_polymer_entity_types', nargs='*', help='Selected polymer entity types', default=None) #[‘Nucleic acid (only)’,‘Other’, ‘Protein/NA’]
    parser.add_argument('--pdbx_keywords', nargs='*', help='PDBx keywords', default=None) #['RNA', 'DNA/RNA', 'RIBOSOME', 'RIBOZYME']
    parser.add_argument('--polymer_type', type=str, help='Polymer type', default='polyribonucleotide')
    parser.add_argument('--sequence_length', type=int, help='Sequence length', default=40)
    parser.add_argument('--sequence_identity', type=float, help='Sequence identity', default=50.0)
    parser.add_argument('--auto_download', type=bool, help='Auto download PDB file', default=True)
    parser.add_argument('--alignment_tool', type=str, help='Alignment tool', default='clustal')
    parser.add_argument('--e_value_cmscan', type=float, help='E-value for CMScan', default=10)
    parser.add_argument('--save', type=bool, help='Save filtered data', default=False)
    args = parser.parse_args()

    dataset_name = args.dataset_name
    if DATA_PATH is None:
        raise ValueError('DATA_PATH environment variable is not set.')
    if not os.path.exists(DATA_PATH):
        os.makedirs(DATA_PATH)
    os.environ['DATA_PATH'] = DATA_PATH + f'/{dataset_name}'

    dc = DatasetCreator(
        dataset_name=args.dataset_name,
        rfam_cm_path=args.rfam_cm_path,
        structure_determination_methodology=args.structure_determination_methodology,
        rcsb_entity_polymer_type=args.rcsb_entity_polymer_type,
        dstart=args.dstart,
        dend=args.dend,
        download_all=args.download_all,
        exptl_method=args.exptl_method,
        resolution=args.resolution,
        year_range=args.year_range,
        polymer_entity_instance_count=args.polymer_entity_instance_count,
        polymer_entity_count_RNA=args.polymer_entity_count_RNA,
        selected_polymer_entity_types=args.selected_polymer_entity_types,
        pdbx_keywords=args.pdbx_keywords,
        polymer_type=args.polymer_type,
        sequence_length=args.sequence_length,
        sequence_identity=args.sequence_identity,
        auto_download=args.auto_download,
        alignment_tool=args.alignment_tool,
        e_value_cmscan=args.e_value_cmscan,
        save=args.save
    )

    dc.run()
    utils.get_final_fam_pdb_chain_csv(dc.clean_tblout_path)
    utils.save_sequences(dc.final_fasta_path, dc.sequences_path)





            

if __name__ == '__main__':
    # Configure logging
    main()


