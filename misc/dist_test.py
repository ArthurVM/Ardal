import ardal
import time
import numpy as np
import json
import pandas as pd
from scipy.spatial.distance import pdist, squareform

def npy_hamming(npy_matrix, f):
    # return np.array(squareform([int(i/f) for i in pdist(npy_matrix, metric='hamming')]))
    m=len(npy_matrix[0])
    n=len(npy_matrix)

    dist_np = np.empty([n,n])
    for i, row_i in enumerate(npy_matrix):
        for j, row_j in enumerate(npy_matrix):
            dist_np[i,j] = int(m-np.sum(row_i == row_j))
    return dist_np


def dneigh(npy_matrix, q, epsilon):
    q_row = npy_matrix[q]

    return np.array([i for i, d in enumerate(q_row) if i != q and d <= epsilon])


npy_matrix = np.ascontiguousarray(np.load("./data/BG_pan_matrix.npy"))
with open("./data/BG_pan_headers.json", 'r') as fin:
    json_headers = json.load(fin)
# npy_matrix = np.ascontiguousarray([[0, 1, 1, 0],[1, 0, 0, 1], [0, 1, 0, 1], [1, 0, 0, 1]])
# npy_matrix, json_headers = generateMatrix(sample_n=1000, allele_n=500)
# npy_matrix = np.array([[0, 1, 1, 0],[1, 0, 0, 1], [0, 1, 0, 1], [1, 0, 0, 1]])
# npy_matrix = np.array(\
#     [[1, 1, 1, 1, 1], 
#      [0, 1, 1, 1, 1], 
#      [0, 0, 1, 1, 1], 
#      [0, 0, 0, 1, 1], 
#      [0, 0, 0, 0, 1], 
#      [0, 0, 0, 0, 0], 
#      [0, 1, 0, 1, 0], 
#      [1, 0, 1, 0, 1]])

_n = len(npy_matrix)

# print(npy_matrix)
print(_n*(_n-1)/2)

## NUMPY TEST
# npy_s = time.time()
# npy_d_out = npy_hamming(npy_matrix, 1/len(npy_matrix[0]))
# npy_e = time.time()
# print(f"NUMPY\nt={npy_e-npy_s}\n")
# print(npy_d_out)

## ARDAL TEST
ardmat = ardal.AlleleMatrix(npy_matrix)
ard_s = time.time()
ard_d_out = np.array(squareform(ardmat.hamming()))
ard_e = time.time()
print(f"ARDAL\nt={ard_e-ard_s}\n")
# print(ard_d_out)

# if (ard_d_out == npy_d_out).all():
#     print("checks PASSED")
# else:
#     print("checks FAILED")
#     print(ard_d_out == npy_d_out)


## compare a manual extraction from the hamming dist matrix to the neighbourhood method

## batch comparison
# truth_box = []
# for epsilon in range(0, 60000, 10000):
#     print(f"epsilon: {epsilon}")
#     for q in range(len(ard_d_out)):
#         loc_n = dneigh(ard_d_out, q, epsilon)
#         ard_n = ardmat.neighbourhood(q, epsilon)
#         if len(loc_n) == 0 and len(ard_n) == 0:
#             truth_box.append(True)
#         elif len(loc_n) != len(ard_n):
#             truth_box.append(False)
#         else:
#             truth_box.append(np.array_equal(loc_n, ard_n))

# print(all(i for i in truth_box))


## SAVE CSVs
# ard_d_df = pd.DataFrame(ard_d_out, columns=json_headers["guids"], index=json_headers["guids"])
# npy_d_df = pd.DataFrame(npy_d_out, headers=json_headers["guids"], index=json_headers["guids"])
# ard_d_df.to_csv("BG_ardalC_hamming.csv")