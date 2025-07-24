import sys
import pandas as pd
import numpy as np
import scipy as sp
import json


class ArdalParser:
    """ Parses Ardal data from CSV, Parquet, or NPY/JSON file pairs.
    """

    def __init__( self,
                 input_data_structure,
                 file_format=None,
                 prefix=None ):
        """ ArdalParser constructor
        """
        super(ArdalParser, self).__init__()

        self.input_data = input_data_structure
        self.file_format = file_format.lower() if file_format else None ## added handling of lowercase and none.
        self.matrix = None
        self.headers = None

        self._parse()  ## Parse the data upon object creation


    def _parse( self ) -> int:
        """ Parses the data based on the specified file format.
        """
        
        ## try to infer file format from extensions
        if self.file_format == None:

            if isinstance(self.input_data, list):
                if len(self.input_data) == 2:

                    if isinstance(self.input_data[0], str) and isinstance(self.input_data[1], str):
                        ## sort inputs by extension
                        self.input_data = sorted(self.input_data, key=lambda x: x.split(".")[-1])
                        formats = [i.split(".")[-1] for i in self.input_data]
                        if formats == ["json", "npy"] or formats == ["json", "npz"]:
                            self.file_format = formats[1]
                        else:
                            raise ValueError("List must contain json and npy/npz file paths.")
                    
                    ## handle data passed as variables
                    elif isinstance(self.input_data[0], np.ndarray) and isinstance(self.input_data[1], dict):
                        self.matrix = np.ascontiguousarray(self.input_data[0])
                        self.headers = self.input_data[1]
                        return 0
                    elif isinstance(self.input_data[0], dict) and isinstance(self.input_data[1], np.ndarray):
                        self.headers = self.input_data[0]
                        self.matrix = np.ascontiguousarray(self.input_data[1])
                        return 0
                    
                    else:
                        raise ValueError("List must contain json and npy/npz file paths.")

                else:
                    raise ValueError("List must contain json and npy/npz file paths.")
            
            else:
                self.file_format = self.input_data.split(".")[-1]

        ## check given format
        if self.file_format not in ["csv", "parquet", "json", "npy", "npz"]:
            raise ValueError("Provided format must be csv, parque, json, npy or npz.")

        if self.file_format == "csv":
            self._csvLoader()

        elif self.file_format == "parquet":
            self._parquetLoader()

        ## handle NPY/JSON pair
        elif self.file_format == "npy":
            self._npyLoader()
        
        elif self.file_format == "npz":
            self._npzLoader()

        return 0


    def _csvLoader(self, delimiter=","):
        """ Parse CSV files.
        """
        print(f"Loading '{self.input_data}' as a csv.")
        pd_df = pd.read_csv(self.input_data, index_col=0, delimiter=delimiter)
        self.matrix = pd_df.values
        self.headers = {
            "guids" : list(pd_df.index),
            "alleles" : list(pd_df.columns)
        }

    
    def _parquetLoader(self):
        """ Parse parquet files.
        """
        print(f"Loading '{self.input_data}' as a parquet.")
        pq_df = pd.read_parquet(self.input_data, engine="fastparquet")
        self.matrix = pq_df.values
        self.headers = {
            "guids" : list(pq_df.index),
            "alleles" : list(pq_df.columns)
        }


    def _npyLoader(self):
        """ Parse npy JSON pairs.
        """
        print(f"Loading '{self.input_data}' as a npy/JSON pair.")
        headers_json, matrix_npy = self.input_data
    
        try:
            self.matrix = np.ascontiguousarray(np.load(matrix_npy))

        except FileNotFoundError:  ## handle missing matrix file
            print(f"Error: Matrix file '{matrix_npy}' not found.")
            sys.exit(101)
            
        except ValueError:  ## handle value errors
            print(f"Error: Invalid data in matrix file '{matrix_npy}'.")
            sys.exit(101)
        
        except Exception as e:  ## catch generic exceptions
            print(e)
            sys.exit(101)
        
        try:
            with open(headers_json, "r") as f:
                self.headers = json.load(f)
        
        except FileNotFoundError:  ## handle missing headers file
            print(f"Error: Header file '{headers_json}' not found.")
            sys.exit(101)

        except TypeError:  ## handle type errors
            print(f"Error: Invalid data type in headers file '{headers_json}'.")
            sys.exit(101)

        except json.JSONDecodeError:  ## handle JSON errors
            print(f"Error: Invalid JSON in headers file '{headers_json}'.")
            sys.exit(101)
        
        except Exception as e:  ## catch generic exceptions
            print(e)
            

        if len(self.matrix) != len(self.headers["guids"]) or len(self.matrix[0]) != len(self.headers["alleles"]):
            raise ValueError(f"Dimension mismatch between matrix array {self.matrix.shape} and headers (rows (guids): {len(self.headers['guids'])}, cols (alleles): {len(self.headers['alleles'])}).")
        

    def _npzLoader(self):
        """ Parse npz JSON pairs.
        """
        print(f"Loading '{self.input_data}' as a npz/JSON pair.")
        headers_json, matrix_npz = self.input_data

        try:
            sparse_matrix = np.ascontiguousarray(sp.sparse.load_npz(matrix_npz).toarray())

            ## check the matrix is C contiguous in memory
            ## this should not happen here since it is explicity loaded as contiguous above
            ## however the criticality of this necessitates a doubly defensive approach
            if sparse_matrix.flags['C_CONTIGUOUS']:
                self.matrix = sparse_matrix
            else:
                self.matrix = np.ascontiguousarray(sparse_matrix)

        except FileNotFoundError:  ## handle missing matrix file
            print(f"Error: Matrix file '{matrix_npz}' not found.")
            sys.exit(102)
        except ValueError as e:  ## handle value errors
            print(f"Error loading matrix from '{matrix_npz}': {e}")
            sys.exit(103)
        except Exception as e:  ## catch generic exceptions
            print(f"An unexpected error occurred while loading '{matrix_npz}': {e}")
            sys.exit(104)

        try:
            with open(headers_json, "r") as f:
                self.headers = json.load(f)

        except FileNotFoundError:  ## handle missing headers file
            print(f"Error: Header file '{headers_json}' not found.")
            sys.exit(101)
        except TypeError:  ## handle type errors
            print(f"Error: Invalid data type in headers file '{headers_json}'.")
            sys.exit(101)
        except json.JSONDecodeError:  ## handle JSON errors
            print(f"Error: Invalid JSON in headers file '{headers_json}'.")
            sys.exit(101)
        except Exception as e:  ## catch generic exceptions
            print(f"An unexpected error occurred while loading '{headers_json}': {e}")
            sys.exit(101)

        if len(self.matrix) != len(self.headers["guids"]) or len(self.matrix[0]) != len(self.headers["alleles"]):
            raise ValueError(f"Dimension mismatch between matrix array {self.matrix.shape} and headers (rows (guids): {len(self.headers['guids'])}, cols (alleles): {len(self.headers['alleles'])}).")

