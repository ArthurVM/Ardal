import os
import sys
import argparse
import pysam
import yaml
import json
import numpy as np
import pandas as pd
from collections import defaultdict, namedtuple


class VCF(object):
    """ Class for handling VCF objects. """


    def __init__( self, vcf_path : str ):
        """ VCF constructor
        """
        super(VCF, self).__init__()

        self.__vcf = pysam.VariantFile(vcf_path)
        self.__samples = list(self.__vcf.header.samples)
        self.__records = self.__vcf.fetch()


    def __parse( self, qual : int=30, pos_collapse : bool=False, het_collapse : bool=True ):
        """
        """
        allele_dict = defaultdict(dict)

        for record in self.__vcf.fetch():
            for sample_id in self.__samples:
                try:
                    call = record.samples[sample_id]
                except KeyError:
                    continue

                if call['GT'] == (0, 0):
                    allele_dict[f"{record.chrom}.{record.pos}.{record.ref}.{record.alt}"][sample_id] = 0 

                ## Handles multi allelic sites
                for allele_index in call['GT']:
                    if allele_index != 0 and allele_index != None and len(record.alleles[allele_index]) == 1 and record.qual >= qual:
                        allele_id = f"{record.chrom}.{record.pos}.{record.ref}.{record.alleles[allele_index]}"

                        if allele_id not in allele_dict:
                            for s in self.__samples:
                                allele_dict[allele_id][s] = 0 

                        allele_dict[allele_id][sample_id] = 1