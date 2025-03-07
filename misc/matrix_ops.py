import numpy as np
import random
import ctypes

# Load the C++ shared library
lib = ctypes.CDLL('./matrix_ops.so')

# Define the argument and return types for the C++ function
lib.access_matrix.argtypes = [
    ctypes.POINTER(ctypes.c_double),  # Matrix data pointer
    ctypes.c_size_t,  # Number of rows
    ctypes.c_size_t,  # Number of cols
    ctypes.POINTER(ctypes.c_size_t),  # Coordinate vector (flattened)
    ctypes.c_size_t,  # Number of coordinates
    ctypes.POINTER(ctypes.c_double),  # Output vector
]
lib.access_matrix.restype = None  # No return value

def access_matrix(matrix, coordinates):
    rows, cols = matrix.shape
    num_coords = len(coordinates)

    print(matrix)

    # Flatten the coordinates for C++ (as [x1, y1, x2, y2, ...])
    flat_coords = np.array(coordinates, dtype=np.uint64).flatten()

    # Allocate output array
    output = np.zeros(len(coordinates), dtype=np.float64)

    # Call the C++ function
    lib.access_matrix(
        matrix.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
        rows,
        cols,
        flat_coords.ctypes.data_as(ctypes.POINTER(ctypes.c_size_t)),
        num_coords,
        output.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
    )

    return output


def generateMatrix(sample_n=1000, allele_n=50000):

    ## 1. Create the Matrix (random 0s and 1s)
    npy_matrix = np.ascontiguousarray(np.random.randint(0, 2, size=(sample_n, allele_n)))  # Adjust as needed for your desired distribution

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

matrix, headers = generateMatrix(sample_n=10000, allele_n=50000)
coords = randcoords(matrix, 1000)

print(f"Matrix size: {len(matrix[0])}x{len(matrix)}")
print(f"Queries: {len(coords)}")
print()

result = access_matrix(matrix, coords)
# print(result)