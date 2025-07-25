""" utilities.py
This module provides utility functions for working with GUIDs and alleles in the Ardal framework.
"""


def checkGUIDs( guids : list,
                headers : dict,
                filter : bool = False ) -> list:
    """ Check guids are present within the matrix and construct a list of present guids to proceed with
    """
    absent_guids = []
    present_guids = []
    for id in guids:
        if id not in headers["guids"]:
            absent_guids.append(id)
        else:
            present_guids.append(id)

    if filter:
        return present_guids

    if len(absent_guids) > 0:
        raise ValueError(f"guids {absent_guids} not found in allele matrix.")


def checkAlleles( alleles : list,
                  headers : dict,
                  filter : bool = False ) -> list:
    """ Check alleles are present within the matrix and construct a list of present alleles to proceed with
    """
    absent_alleles = []
    present_alleles = []
    for allele_id in alleles:
        if allele_id not in headers["alleles"]:
            absent_alleles.append(allele_id)
        else:
            present_alleles.append(allele_id)

    if filter:
        return present_alleles
    
    if len(absent_alleles) > 0:
        raise ValueError(f"alleles {absent_alleles} not found in allele matrix.")


def encodeGuid( guid : str,
                headers : dict ) -> int:
    """ Encode a GUID to its corresponding index in the headers dictionary.
    """
    return headers["guids"].index(guid)


def decodeGuid( row_coord : int,
                headers : dict ) -> str:
    """ Decode a row coordinate to its corresponding GUID in the headers dictionary.
    """
    return headers["guids"][row_coord]


def encodeAllele( allele : str,
                 headers : dict ) -> int:
    """ Encode an allele to its corresponding index in the headers dictionary.
    """
    return headers["alleles"].index(allele)


def decodeAllele( col_coord : int,
                 headers : dict ) -> str:
    """ Decode a column coordinate to its corresponding allele in the headers dictionary.
    """
    return headers["alleles"][col_coord]