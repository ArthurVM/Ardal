import sys
import pandas as pd
import numpy as np
import json


class ArdalParser:
    """ Parses Ardal data from CSV, Parquet, or NPY/JSON file pairs.
    """

    def __init__(self, input_file_structure, file_format=None, prefix=None):
        """ ArdalParser constructor
        """
        super(ArdalParser, self).__init__()

        self.input_file = input_file_structure
        self.file_format = file_format.lower() if file_format else None ## added handling of lowercase and none.
        self.matrix = None
        self.headers = None

        self._parse()  ## Parse the data upon object creation


    def _parse(self) -> int:
        """ Parses the data based on the specified file format.
        """
        
        ## try to infer file format from extensions
        if self.file_format == None:

            if isinstance(self.input_file, list):
                if len(self.input_file) == 2:

                    if isinstance(self.input_file[0], str) and isinstance(self.input_file[1], str):
                        formats = [i.split(".")[-1] for i in self.input_file]
                        if all(x in ["json", "npy"] for x in formats):
                            self.file_format = formats[0]
                    
                    ## handle data passed as variables
                    elif isinstance(self.input_file[0], np.ndarray) and isinstance(self.input_file[1], dict):
                        self.matrix = np.ascontiguousarray(self.input_file[0])
                        self.headers = self.input_file[1]
                        return 0
                    elif isinstance(self.input_file[0], dict) and isinstance(self.input_file[1], np.ndarray):
                        self.headers = self.input_file[0]
                        self.matrix = np.ascontiguousarray(self.input_file[1])
                        return 0
                    
                    else:
                        raise ValueError("List must contain json and npy file paths.")

                else:
                    raise ValueError("List must contain json and npy file paths.")
            
            else:
                self.file_format = self.input_file.split(".")[-1]

        ## check given format
        if self.file_format not in ["csv", "parquet", "json", "npy"]:
            raise ValueError("Provided format must be csv, parque, json or npy")

        if self.file_format == "csv":
            self._csvLoader()

        elif self.file_format == "parquet":
            self._parquetLoader()

        ## handle NPY/JSON pair
        elif self.file_format == "npy" or self.file_format == "json":
            self._npyLoader()

        return 0


    def _csvLoader(self, delimiter=","):
        """ Parse CSV files.
        """
        print(f"Loading '{self.input_file}' as a csv.")
        pd_df = pd.read_csv(self.input_file, index_col=0, delimiter=delimiter)
        self.matrix = pd_df.values
        self.headers = {
            "guids" : list(pd_df.index),
            "alleles" : list(pd_df.columns)
        }

    
    def _parquetLoader(self):
        """ Parse parquet files.
        """
        print(f"Loading '{self.input_file}' as a parquet.")
        pq_df = pd.read_parquet(self.input_file, engine="fastparquet")
        self.matrix = pq_df.values
        self.headers = {
            "guids" : list(pq_df.index),
            "alleles" : list(pq_df.columns)
        }


    def _npyLoader(self):
        """ Parse npy JSON pairs.
        """
        print(f"Loading '{self.input_file}' as a npy/JSON pair. GOT HERE")
        if self.file_format == "npy":
                matrix_npy, headers_json = self.input_file
        elif self.file_format == "json":
            headers_json, matrix_npy = self.input_file
    
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