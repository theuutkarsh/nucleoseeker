import os
import typing
import pandas as pd
import requests
import logging
from nucleoseeker.columns import COLUMNS

DATA_PATH = os.environ.get('DATA_PATH')


class DatasetDownload:
    """
    Download RNA dataset from RCSB PDB database.
    It is used to download the dataset from the RCSB PDB database.
    It uses the search API to get the list of PDB IDs and then uses the GraphQL API to get data for each PDB ID.
    The data is then converted into a DataFrame with the columns specified in the COLUMNS variable.
    The dataframe can be accessed using the `df` attribute and is also saved as a CSV file.
    
    Args:
        structure_determination_methodology (str): Structure determination methodology.
        rcsb_entity_polymer_type (str): Polymer type, possible values are `"Protein", "DNA", "RNA" ,"NA-hybrid", "Other"`
        If "Protein" or any other value is used, the tool will likely work but cmscan will fail as it is designed to work with RNA.
        dstart (int): Start index for pagination.
        dend (int): End index for pagination.
        download_all (bool): Download all entries.
        
    Attributes:
        SEARCH_API_BASE_URI (str): Base URI for search API.
        DATA_API_BASE_URI_GRAPHQL (str): Base URI for data API.
        dstart (int): Start index for pagination.
        dend (int): End index for pagination.
        structure_determination_methodology (str): Structure determination methodology.
        rcsb_entity_polymer_type (str): Polymer type.
        data_path (str): Path to save the data.
        df (pd.DataFrame): Dataframe of the dataset.
        
    Methods:
        get_search_api_query: Get search API query.
        get_graphql_query: Get GraphQL query.
        get_pdb_list: Get list of PDB IDs.
        get_data_for_each_pdb: Get data for each PDB ID.
        get_data_as_df: Get data as DataFrame.
        save_data_as_csv: Save data as CSV.
    
    NOTE: If you wish to supply multiple values for structure_determination_methodology or rcsb_entity_polymer_type,
    you would have to change the search query to include multiple values.
    """
    SEARCH_API_BASE_URI = 'https://search.rcsb.org/rcsbsearch/v2/query'
    DATA_API_BASE_URI_GRAPHQL = 'https://data.rcsb.org/graphql'

    def __init__(
        self,
        structure_determination_methodology: str = 'experimental',
        rcsb_entity_polymer_type: str = 'RNA',
        dstart: int = 0,
        dend: int = 25,
        download_all: bool = False
    ):
        self.dstart = dstart
        self.dend = dend
        if download_all:
            self.dstart = 0
            self.dend = 10000
        self.structure_determination_methodology = structure_determination_methodology
        self.rcsb_entity_polymer_type = rcsb_entity_polymer_type
        self.data_path = os.environ.get('DATA_PATH')
        self.raw_data_path = os.path.join(
                self.data_path,
                f"raw_{self.structure_determination_methodology}_{self.rcsb_entity_polymer_type}_{self.dstart}_{self.dend}.csv"
            )
        

        if os.path.exists(self.raw_data_path):
            self.df = pd.read_csv(self.raw_data_path)
        else:
            if not os.path.exists(self.data_path):
                os.makedirs(self.data_path)
            self.df = self.get_data_as_df()
            self.save_data_as_csv(
                os.path.join(
                    self.data_path,
                    self.raw_data_path
                )
            )
        
    def get_search_api_query(self):
        """
        Get search API query, this is used to search for entries in the RCSB PDB database.
        It filters entries based on structure determination methodology and polymer type.
        This query is used to get the list of PDB IDs.
        
        Returns:
            dict: Search API query.
        NOTE: If you wish to supply multiple values for structure_determination_methodology or rcsb_entity_polymer_type,
        you would have to change the search query to include multiple values.
        """
        query = {
            "query": {
                "type": "group",
                "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                    "attribute": "rcsb_entry_info.structure_determination_methodology",
                    "operator": "exact_match",
                    "value": f"{self.structure_determination_methodology}"
                    }
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                    "attribute": "entity_poly.rcsb_entity_polymer_type",
                    "operator": "exact_match",
                    "negation": False,
                    "value": f"{self.rcsb_entity_polymer_type}"
                    }
                }
                ],
                "logical_operator": "and",
                "label": "text"
            },
            "return_type": "entry",
            "request_options": {
                "paginate": {
                "start": self.dstart,
                "rows": self.dend
                },
                "results_content_type": [
                "experimental"
                ],
                "sort": [
                {
                    "sort_by": "score",
                    "direction": "desc"
                }
                ],
                "scoring_strategy": "combined"
            }
        }
        return query

    def get_graphql_query(self, pdb_ids: list):
        """
        Get GraphQL query, this is used to get data for each PDB ID.
        
        Args:
            pdb_ids (list): List of PDB IDs.
            
        Returns:
            dict: GraphQL query."""
        query = """
            query GetPDBEntries($entryIDs: [String!]!) {
            entries(entry_ids: $entryIDs) 
            {
                rcsb_id
                exptl {
                method
                }
                rcsb_accession_info {
                initial_release_date
                }
                rcsb_entry_info {
                deposited_polymer_entity_instance_count
                polymer_entity_count_RNA
                resolution_combined
                selected_polymer_entity_types
                }
                struct_keywords {
                pdbx_keywords
                }
            }
        }
        """
        variables = {
            "entryIDs": pdb_ids
        }

        query_json = {
            "query": query,
            "variables": variables
        }
        return query_json
    
    def get_pdb_list(self):
        """
        This method fetches the list of PDB IDs from the RCSB PDB database using the search API.
        
        Returns:
            list: List of PDB IDs.
        """
        try:
            response = requests.post(self.SEARCH_API_BASE_URI, json=self.get_search_api_query())
            response.raise_for_status()
            logging.info('Successfully fetched PDB list')
        except requests.exceptions.RequestException as e:
            logging.error(f'Failed to fetch PDB list: {e}')
            raise
        try:
            content = response.json()
        except ValueError as e:
            raise ValueError(f"The query failed. This might be because of wrong argument for 'structure_determination_methodology' or 'rcsb_entity_polymer_type'. Please check these arguments '{self.structure_determination_methodology}' and '{self.rcsb_entity_polymer_type}'.")
        pdb_dataframe = pd.DataFrame.from_dict(content.get('result_set', []))
        pdb_list = pdb_dataframe['identifier'].tolist()
        return pdb_list
    
    def get_data_for_each_pdb(self):
        """
        This method fetches data for each PDB ID using the GraphQL API.
        
        Returns:
            json: combined data for each PDB ID.
        """
        pdb_list = self.get_pdb_list()
        graphql_query = self.get_graphql_query(pdb_list)
        try:
            response = requests.post(self.DATA_API_BASE_URI_GRAPHQL, json=graphql_query)
            response.raise_for_status()
            logging.info('Successfully fetched data for each PDB')
        except requests.exceptions.RequestException as e:
            logging.error(f'Failed to fetch data for each PDB: {e}')
            raise

        return response.json().get('data', {}).get('entries', [])
    
    def get_data_as_df(self):
        """
        This method converts the data for each PDB ID into a DataFrame.
        
        Returns:
            pd.DataFrame: Dataframe of the whole dataset obtained from the RCSB PDB database using the search and GraphQL APIs.
        """
        data = self.get_data_for_each_pdb()
        df = pd.json_normalize(data)
        df.columns = COLUMNS
        df['rcsb_id'] = df['rcsb_id'].str.lower()
        df['exptl_method'] = df['exptl_method'].apply(lambda x: x[0]['method'] if isinstance(x, list) and x else None)
        df['resolution'] = df['resolution'].apply(lambda x: float(x[0]) if isinstance(x, list) and x else None)
        df['release_date'] = df['release_date'].str.split('-').str[0]
        return df
    
    def save_data_as_csv(self, path: typing.Union[str, os.PathLike]):
        """
        Save data as CSV.
        
        Args:
            path (str, os.PathLike): Path to save the data.
        
        Returns:
            None
        """
        self.df.to_csv(path, index=False)
        logging.info(f'Data saved as CSV at {path}')


