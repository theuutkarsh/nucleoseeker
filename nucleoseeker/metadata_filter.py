import pandas as pd
import typing
from datetime import date
from nucleoseeker.columns import COLUMNS
CURRENT_YEAR = date.today().year

class MetadataFilter:
    """
    Filter the dataset at the structure level.
    Using this class we can filter the dataset based on structure level attributes.
    It works on the columns specified in the COLUMNS variable and filters the dataset based on the given criteria.
    
    Args:
        exptl_method (list): Experimental method.
        resolution (float): Resolution.
        year_range (int or list): Year range.
        polymer_entity_instance_count (int): Polymer entity instance count.
        polymer_entity_count_RNA (int): Polymer entity count RNA.
        selected_polymer_entity_types (list): Selected polymer entity types.
        pdbx_keywords (list): PDBx keywords.
        
    Attributes:
        exptl_method (list): Experimental method.
        resolution (float): Resolution.
        year_range (int or list): Year range.
        polymer_entity_instance_count (int): Polymer entity instance count.
        polymer_entity_count_RNA (int): Polymer entity count RNA.
        selected_polymer_entity_types (list): Selected polymer entity types.
        pdbx_keywords (list): PDBx keywords.
        
    Methods:
        _check_dataframe: Check if the dataframe columns match the expected columns.
        apply_exptl_method_filter: Apply filter based on experimental method.
        apply_resolution_filter: Apply filter based on resolution.
        apply_year_range_filter: Apply filter based on year range.
        apply_polymer_entity_instance_count_filter: Apply filter based on polymer entity instance count.
        apply_polymer_entity_count_RNA_filter: Apply filter based on polymer entity count RNA.
        apply_selected_polymer_entity_types_filter: Apply filter based on selected polymer entity types.
        apply_pdbx_keywords_filter: Apply filter based on PDBx keywords.
        apply_filters: Apply all the filters.
        
    """
    def __init__(
            self,
            exptl_method: typing.List[str],
            resolution: float,
            year_range: typing.List[int],
            polymer_entity_instance_count: int,
            polymer_entity_count_RNA: int,
            selected_polymer_entity_types: typing.List[str],
            pdbx_keywords: typing.List[str],

        ):
        self.exptl_method = exptl_method
        self._check_resolution_validity(resolution)
        self.resolution = resolution
        # assert len(year_range) == 2 or len(year_range) == 1, "Year range must be a list of two elements or a single element."
        if year_range is None:
            self.year_range = year_range
        elif len(year_range) == 1:
            self._check_year_validity(year_range[0])
            self.year_range = year_range
        elif isinstance(year_range, list):
            self._check_year_validity(year_range[0])
            self._check_year_validity(year_range[1])
            if year_range[0] < year_range[1]:
                self.year_range = year_range
            else:
                raise ValueError("The first element of the year_range must be smaller than the second.")

        self._check_polymer_entity_instance_count_validity(polymer_entity_instance_count)
        self.polymer_entity_instance_count = polymer_entity_instance_count
        self._check_polymer_entity_count_RNA_validity(polymer_entity_count_RNA)
        self.polymer_entity_count_RNA = polymer_entity_count_RNA
        self.selected_polymer_entity_types = selected_polymer_entity_types
        self.pdbx_keywords = pdbx_keywords

    def _check_year_validity(self, year: int):
        if year < 1971:
            raise ValueError("RCSB PDB started in 1971, so the year must be greater than that.")
        if year > CURRENT_YEAR:
            raise ValueError(f"Year must be less than or equal to the current year {CURRENT_YEAR}.")

    def _check_resolution_validity(self, resolution: float):
        if resolution is not None:
            if resolution < 0:
                raise ValueError("Resolution must be a positive number.")
            if resolution >= 50:
                raise ValueError("Resolution must be less than 50.")

    def _check_polymer_entity_instance_count_validity(self, count: int):
        if count is not None:
            if count < 0:
                raise ValueError("Polymer entity instance count must be a positive number.")

    def _check_polymer_entity_count_RNA_validity(self, count: int):
        if count is not None:
            if count < 0:
                raise ValueError("Polymer entity count RNA must be a positive number.") 

    def _check_dataframe(self, df: pd.DataFrame):
        if df.columns.tolist() != COLUMNS:
            raise ValueError(f'Columns of the dataframe do not match the expected columns: {COLUMNS}')

    def apply_exptl_method_filter(self, df: pd.DataFrame):
        """
        Apply filter based on experimental method.
        
        Args:
            df (pd.DataFrame): Input dataframe.
            
        Returns:
            pd.DataFrame: Filtered dataframe."""
        self._check_dataframe(df)
        if self.exptl_method is not None:
            if self.exptl_method not in df['exptl_method'].unique():
                raise ValueError(f"Experimental method {self.exptl_method} not found in the dataframe. Please check the spelling or use a different method.")
            else:
                df = df[df['exptl_method'].isin(self.exptl_method)]
                df.reset_index(drop=True, inplace=True)
        return df
    
    def apply_resolution_filter(self, df: pd.DataFrame):
        """
        Apply filter based on resolution.
        
        Args:
            df (pd.DataFrame): Input dataframe.
            
        Returns:
            pd.DataFrame: Filtered dataframe."""
        self._check_dataframe(df)
        if self.resolution is not None:
            if df['resolution'].isnull().all():
                raise ValueError("Resolution column is empty. This might be because the experimental method does not provide resolution information.")
            else:
                df = df[df['resolution'] <= self.resolution]
                df.reset_index(drop=True, inplace=True)
        return df
    
    def apply_year_range_filter(self, df: pd.DataFrame):
        """
        Apply filter based on year range. It can be a single year or a range of two years.
        
        Args:
            df (pd.DataFrame): Input dataframe.
            
        Returns:
            pd.DataFrame: Filtered dataframe."""
        self._check_dataframe(df)
        if self.year_range is not None:
            df = df.astype({'release_date': 'int'})
            if len(self.year_range) == 1:
                df = df[df['release_date'] <= self.year_range[0]]
            else:
                df = df[df['release_date'].isin(self.year_range)]
            df.reset_index(drop=True, inplace=True)
        return df
    
    def apply_polymer_entity_instance_count_filter(self, df: pd.DataFrame):
        """
        Apply filter based on polymer entity instance count.
        
        Args:
            df (pd.DataFrame): Input dataframe.
            
        Returns:
            pd.DataFrame: Filtered dataframe.
        """
        self._check_dataframe(df)
        if self.polymer_entity_instance_count is not None:
            df = df[df['polymer_entity_instance_count'] == self.polymer_entity_instance_count]
            df.reset_index(drop=True, inplace=True)
        return df
    
    def apply_polymer_entity_count_RNA_filter(self, df: pd.DataFrame):
        """
        Apply filter based on polymer entity count RNA.
        
        Args:
            df (pd.DataFrame): Input dataframe.
            
        Returns:
            pd.DataFrame: Filtered dataframe.
        """
        self._check_dataframe(df)
        if self.polymer_entity_count_RNA is not None:
            df = df[df['polymer_entity_count_RNA'] == self.polymer_entity_count_RNA]
            df.reset_index(drop=True, inplace=True)
        return df
    
    def apply_selected_polymer_entity_types_filter(self, df: pd.DataFrame):
        """
        Apply filter based on selected polymer entity types.
        
        Args:
            df (pd.DataFrame): Input dataframe.
            
        Returns:
            pd.DataFrame: Filtered dataframe.
        """
        self._check_dataframe(df)
        if self.selected_polymer_entity_types is not None:
            if self.selected_polymer_entity_types not in df['selected_polymer_entity_types'].unique():
                raise ValueError(f"Selected polymer entity types {self.selected_polymer_entity_types} not found in the dataframe. Please check the spelling or use a different polymer type.")
            else:
                df = df[df['selected_polymer_entity_types'].isin(self.selected_polymer_entity_types)]
                df.reset_index(drop=True, inplace=True)
        return df
    
    def apply_pdbx_keywords_filter(self, df: pd.DataFrame):
        """
        Apply filter based on PDBx keywords.
        
        Args:
            df (pd.DataFrame): Input dataframe.
            
        Returns:
            pd.DataFrame: Filtered dataframe.
        """
        self._check_dataframe(df)
        if self.pdbx_keywords is not None:
            if self.pdbx_keywords not in df['pdbx_keywords'].unique():
                raise ValueError(f"PDBx keywords {self.pdbx_keywords} not found in the dataframe. Please check the spelling or use a different keyword.")
            else:
                df = df[df['pdbx_keywords'].isin(self.pdbx_keywords)]
                df.reset_index(drop=True, inplace=True)
        return df
    
    def apply_filters(self, df: pd.DataFrame):
        """
        Apply all the filters.
        
        Args:
            df (pd.DataFrame): Input dataframe.
            
        Returns:
            pd.DataFrame: Filtered dataframe.
        """
        df = self.apply_exptl_method_filter(df)
        df = self.apply_pdbx_keywords_filter(df)
        df = self.apply_resolution_filter(df)
        df = self.apply_year_range_filter(df)
        df = self.apply_polymer_entity_instance_count_filter(df)
        df = self.apply_polymer_entity_count_RNA_filter(df)
        df = self.apply_selected_polymer_entity_types_filter(df)
        return df
