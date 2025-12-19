import numpy as np
import pandas as pd
import json
import time
import random
import ardal


def generateMatrix(sample_n=1000, allele_n=50000):

    ## 1. Create the Matrix (random 0s and 1s)
    npy_matrix = np.ascontiguousarray(np.random.randint(0, 2, size=(sample_n, allele_n), dtype="uint8"))  # Adjust as needed for your desired distribution
    
    ## 2. GUIDs (Row Labels) - Simple numerical IDs
    guids = [f"sample_{i}" for i in np.arange(sample_n)]

    ## 3. Alleles (Column Labels) - Simple numerical IDs
    alleles = [f"allele_{i}" for i in np.arange(allele_n)]

    headers_json = {"guids" : guids, "alleles" : alleles}

    return npy_matrix, headers_json


def randcoords(npy_matrix, ncoords=10):
    # npy_matrix = np.load("./data/test_csv_matrix.npy")
    # test_matrix = np.array([[0, 1, 1, 0],[1, 0, 0, 1], [0, 1, 0, 1], [1, 0, 0, 1]])

    # coords = np.array([[0, 0], [0, 1], [0, 2], [0, 3], [2, 1]], dtype=int)

    n_rows, n_cols = len(npy_matrix)-1, len(npy_matrix[0])-1
    coords = []

    # print(npy_matrix, n_rows, n_cols)


    for _ in range(ncoords):
            col = random.randint(0, n_cols - 1)
            row = random.randint(0, n_rows - 1)
            coords.append([row, col])

    return np.array(coords, dtype=int)

read_s = time.time()
# npy_matrix = np.ascontiguousarray(np.load("./data/BG_pan_matrix.npy"))
npy_matrix, json_headers = generateMatrix(sample_n=100, allele_n=2000)
read_e = time.time()

# coords = np.array([[0, 1], [4, 1], [5, 1], [9, 1], [0, 6], [4, 6], [5, 6], [9, 6], [0, 11], [4, 11], [5, 11], [9, 11]], dtype=int)
coords = randcoords(npy_matrix, 1000000)

print(f"Matrix size: {len(npy_matrix[0])}x{len(npy_matrix)}")
print(f"Queries: {len(coords)}")
print(f"Read time: {read_e - read_s}")
print()

# print("NUMPY:")
npy_s = time.time()
npy_out = npy_matrix[coords[:, 0], coords[:, 1]]
# print(npy_out)
npy_e = time.time()
npy_t = npy_e - npy_s

# print("C++:")
c_s = time.time()
result_bool = ardal.access(coords)
c_e = time.time()
c_t = c_e - c_s
# print("Result shape:", result_bool.shape)
# print("Result data:\n", result_bool)

print(f"npy: {npy_t}")
print(f"c++: {c_t}")
print(f"delta: {c_t - npy_t}")

print()

# print(result_bool)
# print(npy_out)
if (result_bool == npy_out).all():
    print("checks PASSED")
else:
    print("checks FAILED")