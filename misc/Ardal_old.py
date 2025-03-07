import pandas as pd
import dask.dataframe as dd
import csv
import numpy as np
import re
from sys import stdout, stderr
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict, namedtuple

from ArdalParser import ArdalParser

class Ardal(object):
    """ Class for handling Ardal-neighbourhood objects. """


    def __init__( self, data_source : str, ref : bool=False, file_format : str=None ):
        """ Ardal constructor
        """
        super(Ardal, self).__init__()

        self.matrix = None
        self.ref = ref

        ## if the csv argument is a dataframe the just assign it as the matrix
        if isinstance(data_source, pd.DataFrame):
            self.matrix = data_source

        ## otherwise, assume it is the path to a csv, parquet, or json/npy pair.
        else:
            # if db.endswith(".csv"):
            #     self.matrix = self._parseMatrix(db, ref)
            # elif db.endswith("parquet"):
            #     self.matrix = self._parseMatrix(db, ref, file_format="parquet")

            parser = ArdalParser(data_source, file_format)

            if parser.matrix is not None:
                self.matrix = parser.matrix
            else:
                ## raise an error if parsing fails to prevent unexpected behaviour down the line
                raise ValueError(f"Failed to parse data from: {data_source}") 


    def _parseMatrix(self, matrix_path: str, ref: bool = False, file_format="csv", chunk_size=10000) -> pd.DataFrame:
        """Parse a sparse presence absence matrix in CSV or Parquet format."""

        if file_format not in ("csv", "parquet"):
            raise ValueError("Invalid file format. Must be 'csv' or 'parquet'.")

        if file_format == "csv":
           wg_matrix = self._csvToDf(matrix_path)

        elif file_format == "parquet":
            wg_matrix = pd.read_parquet(matrix_path)

        if ref:
            if "REF" not in wg_matrix.index:  # Check before assigning
                wg_matrix.loc["REF"] = 0
            else:
                print("REF already in index, skipping assignment.")

        if wg_matrix.index.has_duplicates:
            raise KeyError(f"Multiple instances of GUIDs found in {matrix_path}.")

        ## Convert to uint8 after concat and duplicate check.  Only do this for numeric types.
        numeric_cols = wg_matrix.select_dtypes(include='number').columns
        wg_matrix[numeric_cols] = wg_matrix[numeric_cols].astype(np.uint8)

        return wg_matrix
        
    
    def _csvToDf(self, filepath: str, delimiter: str = ',', quotechar: str = '"', header: int = 0, index_col: int = 0) -> pd.DataFrame:
        """Reads a CSV file and returns a Pandas DataFrame.

        INPUT:
            filepath <str>: The path to the CSV file.
            delimiter <str, optional>: The delimiter used in the CSV. Defaults to ','.
            quotechar <str, optional>: The quote character used in the CSV. Defaults to '"'.
            header <int or list, optional>: Row number(s) to use as the column names, 
                                        or a list of column names, or None for no header. Defaults to 0 (first row).
            index_col <int or str or None, optional>: Column number or name to use as
                        the row index. If None, a default index is used. Defaults to 0.

        OUTPUT:
            df <pandas.DataFrame>: The parsed DataFrame, or None if an error occurs.
        """

        try:
            with open(filepath, 'r', newline='') as csvfile:
                reader = csv.reader(csvfile, delimiter=delimiter, quotechar=quotechar)

                ## Handle header
                if isinstance(header, int):
                    column_names = next(reader) if header is not None else None  # Read the specified header row
                    data_start = header + 1
                    if header > 0:
                        ## insert dummy headers 
                        column_names = ["Column " + str(i+1) for i in range(len(column_names))]

                elif isinstance(header, list):
                    ## pass column names by header argument
                    column_names = header
                    data_start = 0
                    
                elif header is None:  
                    ## No header is set in the header argument
                    column_names = None
                    data_start = 0 

                else:
                    raise TypeError("Invalid header argument. Must be int, list, or None.")

                data = []
                row_count = 0
                for row in reader:
                    if row_count >= data_start or header is None or isinstance(header, list):  # Start appending after header row
                        data.append(row)

                    row_count += 1

            df = pd.DataFrame(data, columns=column_names)

            if index_col is not None and len(df.columns) > index_col:
                df = df.set_index(df.columns[index_col])
            
            return df

        except Exception as e:  # Catch any potential errors during file reading or processing
            print(f"Error reading CSV: {e}")
            return None
    

    def subsetbyGUID( self, guid_list : list ):
        """ Take a list of GUIDs and subset the allele matrix to include only these GUIDs, allowing for standard operations.
        Returns an Ardal object with the subsetted matrix.
        """
        ## Filter rows based on the provided list of indexes
        subset_matrix = self.matrix[self.matrix.index.isin(guid_list)]
    
        ## Remove alleles which do not exist within the subset
        subset_matrix = subset_matrix.loc[:, (subset_matrix != 0).any(axis=0)]

        ## return a new Ardal object
        return Ardal(subset_matrix)

    
    def neighbourhood( self, guid : str, n : int ) -> set:
        """ get the SNP neighbourhood of a GUID
        """
        neighbour = namedtuple( "neighbour", ["guid", "dist"] )

        if guid not in self.matrix.index:
            return f"Variant ID '{guid}' not found in the matrix."

        query_alleles = self.matrix.loc[guid]

        close_variants = set()
        # for variant, alleles in self.matrix.iterrows():
        #     dist = Ardal._absD(query_alleles, alleles)

        for variant in self.matrix.index:
            intersected_snps = self._getAlleleSignature(self.matrix.T[guid].T) ^ self._getAlleleSignature(self.matrix.T[variant].T)
            dist = len(intersected_snps)
            
            if dist <= n:
                close_variants.add(neighbour(variant, dist))

        return set(close_variants)
    

    def toDict( self ) -> dict:
        """ Return a dictionary containing present allele IDs mapped to their guid.
        """
        allele_mapping = defaultdict(list)

        for variant, alleles in self.matrix.iterrows():
            allele_signature = Ardal._getAlleleSignature(alleles)
            allele_mapping[variant] = allele_signature

        return allele_mapping
    

    def unique( self, guids : list ) -> set:
        """ Take a set of guids and return alleles unique to this subset.
        """
        intersection, symmetric = self._getIntersectionAndSymmetric(guids)

        return intersection.difference(symmetric)
    

    def common( self, guids : list ) -> set:
        """ Take a set of guids and return alleles common to this subset.
        """
        intersection, symmetric = self._getIntersectionAndSymmetric(guids)

        return intersection
    

    def allele( self, alleles : list ) -> set:
        """ Take a set of alleles and return all GUIDs containing all of those alleles.
        """

        ## check guids are present within the matrix and construct a list of present guids to proceed with
        present_alleles = self._checkAlleles(alleles)
        allele_count = len(present_alleles)

        ## initialise a dictionary for storing present alleles
        allele_dict = { guid : [] for guid in self.matrix.index.values }

        for allele_id in present_alleles:
            mask = self.matrix[allele_id].astype(bool)
            for guid in self.matrix.loc[mask].index.values:
                allele_dict[guid].append(allele_id)

        return set([guid for guid, alleles in allele_dict.items() if len(alleles) == allele_count])
    

    def old_pairwise(self, guids: list = [], metric: str = "euclidean", chunk_size : int=100 ) -> pd.DataFrame:
        """ Takes a set of GUIDs and calculates the pairwise distance between them, returning a distance matrix.
        Pairwise distance can be calculated using Jaccard, Euclidean, or absolute distance functions.
        If an empty list if provided (as by default) then the pairwise distance of all samples within the matrix will be calculated.
        """

        ## check the specified distance function is valid
        accepted_dist_functions = ["absolute", "euclidean", "jaccard"]
        if metric not in accepted_dist_functions:
            raise ValueError(f"{metric} not an accepted distance function. Accepted distance functions: {accepted_dist_functions}")

        present_guids = self._checkGUIDs(guids)

        if len(guids) == 0:
            subset_matrix = self.matrix
            present_guids = self.matrix.index.values
        elif len(present_guids) == 0:
            raise KeyError(f"No input GUIDs found within the matrix.")
        elif len(present_guids) == 1:
            raise ValueError(f"Cannot construct a distance matrix using fewer than 2 valid GUIDs.")
        else:
            subset_matrix = self.matrix.loc[guids]

        ## initialize an empty distance matrix
        dist_matrix = np.zeros((len(present_guids), len(present_guids)))

        ## Calculate the number of chunks
        n_chunks = len(subset_matrix) // chunk_size + 1

        ## Calculate pairwise distances using vectorized operations in chunks
        for chunk_idx in range(n_chunks):

            ## track progress
            # stdout.write(f"Processing chunk {chunk_idx}/{n_chunks}\r")
            # stdout.flush()
            
            start_idx = chunk_idx * chunk_size
            end_idx = min((chunk_idx + 1) * chunk_size, len(subset_matrix))
            chunk_data = subset_matrix.iloc[start_idx:end_idx]

            if metric == "absolute":
                dist_chunk = np.sum(np.abs(chunk_data.values[:, None] - subset_matrix.values), axis=2)
            elif metric == "euclidean":
                dist_chunk = np.sqrt(((chunk_data.values[:, None, :] - subset_matrix.values) ** 2).sum(axis=2))
            elif metric == "jaccard":
                ## TODO: implement your Jaccard distance function using vectorized operations
                
                for guid_i, a_i in chunk_data.iterrows():
                    for guid_j, a_j in chunk_data.iterrows():
                        print(Ardal._absD(a_i, a_j))
                        
                # intersection = chunk_data.values[:, None] & chunk_data.values[:, None]  # Count the number of common 1s
                # union = chunk_data.values[:, None] | chunk_data.values[:, None]        # Count the number of total 1s
                # jaccard_similarity = intersection / union
                # dist_chunk = 1 - jaccard_similarity

            dist_matrix[start_idx:end_idx] = dist_chunk

        return pd.DataFrame(dist_matrix, columns=present_guids, index=present_guids)


    def pairwise(self, guids: list = [], metric: str = "euclidean", chunk_size : int=100 ) -> pd.DataFrame:
        """ Takes a set of GUIDs and calculates the pairwise distance between them, returning a distance matrix.
        Pairwise distance can be calculated using Jaccard, Euclidean, or absolute distance functions.
        If an empty list if provided (as by default) then the pairwise distance of all samples within the matrix will be calculated.
        """

        ## check the specified distance function is valid
        accepted_dist_functions = ["absolute", "euclidean", "jaccard", "hamming"]
        if metric not in accepted_dist_functions:
            raise ValueError(f"{metric} not an accepted distance function. Accepted distance functions: {accepted_dist_functions}")

        present_guids = self._checkGUIDs(guids)

        if len(guids) == 0:
            subset_matrix = self.matrix
            present_guids = self.matrix.index.values
        elif len(present_guids) == 0:
            raise KeyError(f"No input GUIDs found within the matrix.")
        elif len(present_guids) == 1:
            raise ValueError(f"Cannot construct a distance matrix using fewer than 2 valid GUIDs.")
        else:
            subset_matrix = self.matrix.loc[guids]

        ## calculate distance between two samples using the given distance function
        if metric == "absolute":
            dist_matrix = Ardal._absD(self.matrix)
        elif metric == "euclidean":
            dist_matrix = Ardal._eD(self.matrix)
        elif metric == "jaccard":
            dist_matrix = Ardal._jD(self.matrix)
        elif metric == "hamming":
            dist_matrix = Ardal._hD(self.matrix)

        return pd.DataFrame(dist_matrix, columns=present_guids, index=present_guids)
    

    def removeGUID( self, guid : str ):
        """ Removes a sample from the matrix, defined by a unique GUID
        """
        None


    def matchAlleleNames(self, expression: str) -> list:
        """ Return all allele names that match the given expression with wildcards.
        """
        if not isinstance(expression, str):
            raise ValueError("Expression must be a string.")

        pattern = re.compile(expression.replace('*', '.*'))
        return set([allele for allele in self.matrix.columns if pattern.match(allele)])
    
    
    def _getIntersectionAndSymmetric( self, guids : list ) -> list:
        """ Take a set of guids and return both the allelic intersection and the symmetric set of the union.
        """

        present_guids = self._checkGUIDs(guids)

        if len(present_guids) == 0:
            present_guids = list(self.matrix.index)

        ## get allele signature for the first guid to use as a reference
        intersection = set(Ardal._getAlleleSignature(self.matrix.loc[present_guids[0]]))

        ## list to store (A U B) / A
        symmetric = []

        for variant, alleles in self.matrix.iterrows():
            if variant in present_guids:
                intersection.intersection_update(Ardal._getAlleleSignature(alleles))
            else:
                symmetric.append(Ardal._getAlleleSignature(alleles))
        
        symmetric = set([x for xs in symmetric for x in xs])
        
        return intersection, symmetric
    

    def _checkGUIDs( self, guids : list ) -> list:
        """ Check guids are present within the matrix and construct a list of present guids to proceed with
        """

        present_guids = []
        for id in guids:
            if self._checkIndex(id):
                present_guids.append(id)
            else:
                print(f"{id} not present in allele matrix.")

        return present_guids
    

    def _checkAlleles( self, alleles : list ) -> list:
        """ Check alleles are present within the matrix and construct a list of present alleles to proceed with
        """
        present_alleles = []
        for allele_id in alleles:
            if self._checkCol(allele_id):
                present_alleles.append(allele_id)
            else:
                print(f"{allele_id} not present in allele matrix.")

        return present_alleles

    
    @staticmethod
    def _absD( df : pd.DataFrame ) -> np.array:
        """ return the squareform absolute SNP distance matrix for an input table
        """
        return squareform(pdist(df, metric='cityblock'))

    @staticmethod
    def _hD( df : pd.DataFrame ) -> np.array:
        """ return the squareform Hamming distance matrix for an input table
        """
        return squareform(pdist(df, metric='hamming'))


    @staticmethod
    def _eD( df : pd.DataFrame ) -> np.array:
        """ return the squareform Euclidean distance matrix for an input table
        """
        return squareform(pdist(df, metric='euclidean'))


    @staticmethod
    def _jD( df : pd.DataFrame ) -> np.array:
        """ return the squareform Jaccard distance matrix for an input table
        """
        return squareform(pdist(df, metric='jaccard'))
    

    def _jaccard_distance(X, Y):
        """
        Calculate the Jaccard distance between two sets of binary vectors.

        Parameters:
        - X: 2D numpy array, shape (n_samples, n_features)
            Binary feature matrix for the first set of samples.
        - Y: 2D numpy array, shape (n_samples, n_features)
            Binary feature matrix for the second set of samples.

        Returns:
        - distances: 1D numpy array, shape (n_samples,)
            Jaccard distances between corresponding samples in X and Y.
        """
        intersection = np.sum(X & Y, axis=1)  # Count the number of common 1s
        union = np.sum(X | Y, axis=1)         # Count the number of total 1s
        jaccard_similarity = intersection / union
        distances = 1 - jaccard_similarity    # Convert similarity to distance
        return distances


    @staticmethod
    def _getAlleleSignature( df : pd.DataFrame ) -> set:
        """ Use boolean indexing to filter alleles by presence and return the allele IDs
        """
        return set(df[df.eq(1)].index)
    

    def _checkIndex( self, id : str ) -> bool:
        """ take an index ID and check it exists within the matrix.
        """
        return id in list(self.matrix.index)
    

    def _checkCol( self, id : str ) -> bool:
        """ take a column ID and check it exists within the matrix.
        """
        return id in list(self.matrix.columns.values)