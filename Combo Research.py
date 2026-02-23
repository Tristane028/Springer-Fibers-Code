from itertools import permutations
import numpy as np
import sympy as sp

#The following line is for running the script on my computer specifically, omit this line when running this
#/usr/local/bin/python3.9 "/Users/christan065/Springer Fibers/Combo Research.py"

lambda_func = lambda d: list(range(1, int(d) + 1))

def create_matrix(num_str):
    """
    Creates a matrix where each digit represents the number of columns in each row.
    """
    matrix = [lambda_func(d) for d in reversed(num_str)]
    return matrix


# Ask user for input
num_str = input("Enter the number string: ")
result = create_matrix(num_str)


# Calculate the sum of digits
n = sum(int(digit) for digit in num_str)
def nilpotent_jordan_from_digits(num_str):
    """
    Create a nilpotent Jordan matrix where each digit of num_str
    corresponds to a Jordan block of that size.
   
    Example:
        num_str = 221  -> blocks of sizes [2, 2, 1]
    """
    digits = [int(d) for d in str(num_str)]
    n = sum(digits)
   
    J = np.zeros((n, n), dtype=int)
   
    idx = 0
    for k in digits:
        for i in range(k - 1):
            J[idx + i, idx + i + 1] = 1
        idx += k
   
    return J
J = nilpotent_jordan_from_digits(num_str)
print(J)
# Create a flattened matrix and fill by column
numbers = list(range(1, n + 1))
filled_matrix = [row[:] for row in result]
col_idx = 0
num_idx = 0


for col in range(max(len(row) for row in filled_matrix)):
    for row in range(len(filled_matrix)):
        if col < len(filled_matrix[row]) and num_idx < n:
            filled_matrix[row][col] = numbers[num_idx]
            num_idx += 1



def is_valid_matrix(matrix):
    for col in range(len(matrix[0])):
        column_values = [matrix[row][col] for row in range(len(matrix)) if col < len(matrix[row])]
        if sorted(column_values) != column_values or len(set(column_values)) != len(column_values):
            return False
    return True


def is_column_increasing(matrix):
    max_cols = max(len(row) for row in matrix)
    for col in range(max_cols):
        for row in range(1, len(matrix)):
            if col < len(matrix[row]) and col < len(matrix[row - 1]):
                if matrix[row][col] <= matrix[row - 1][col]:
                    return False
    return True


def generate_matrices(filled_matrix, n):
    unique_numbers = list(range(1, n + 1))
    for perm in permutations(unique_numbers):
        new_matrix = []
        start_index = 0
        for row in filled_matrix:
            new_row = list(perm[start_index:start_index + len(row)])
            new_matrix.append(new_row)
            start_index += len(row)
       

generate_matrices(filled_matrix, n)
def count_inversions(matrix):
    """
    Count inversions in the matrix:
    1. Same row: b left of a, b > a
    2. Different rows: c below d, c in column left of d, c < d
    """
    inversions = 0
   
    # Flatten matrix with position info: (value, row, col)
    positions = []
    for row_idx, row in enumerate(matrix):
        for col_idx, val in enumerate(row):
            positions.append((val, row_idx, col_idx))
   
    # Count inversions
    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            val_i, row_i, col_i = positions[i]
            val_j, row_j, col_j = positions[j]
           
            # Same row: i is to left of j and val_i > val_j
            if row_i == row_j and col_i < col_j and val_i > val_j:
                inversions += 1
           
            # Different rows: j below i, j in column left of i, val_j < val_i
            if row_j > row_i and col_j < col_i and val_j < val_i:
                inversions += 1
   
    return inversions


def generate_matrices(filled_matrix, n):
    unique_numbers = list(range(1, n + 1))
    for perm in permutations(unique_numbers):
        new_matrix = []
        start_index = 0
        for row in filled_matrix:
            new_row = list(perm[start_index:start_index + len(row)])
            new_matrix.append(new_row)
            start_index += len(row)
        if is_valid_matrix(new_matrix) and is_column_increasing(new_matrix):
            print("\nValid matrix:")
            for row in new_matrix:
                print(row)
            inv_count = count_inversions(new_matrix)
            print(f"Inversions: {inv_count}")


generate_matrices(filled_matrix, n)
# Collect matrices by inversion count
max_inversions = 0
matrices_by_inversions = {}


unique_numbers = list(range(1, n + 1))
for perm in permutations(unique_numbers):
    new_matrix = []
    start_index = 0
    for row in filled_matrix:
        new_row = list(perm[start_index:start_index + len(row)])
        new_matrix.append(new_row)
        start_index += len(row)
    if is_valid_matrix(new_matrix) and is_column_increasing(new_matrix):
        inv_count = count_inversions(new_matrix)
        max_inversions = max(max_inversions, inv_count)
        if inv_count not in matrices_by_inversions:
            matrices_by_inversions[inv_count] = []
        matrices_by_inversions[inv_count].append(new_matrix)


# Print results
for i in range(max_inversions + 1):
    print(f"\na{i} (matrices with {i} inversions):")
    if i in matrices_by_inversions:
        print("\n")
        for matrix in matrices_by_inversions[i]:
            for row in matrix:
                print(row)
            print("\n")


temporary = []


for i in range(max_inversions + 1):
    if i in matrices_by_inversions:
        for matrix in matrices_by_inversions[i]:
            matrix_string = ""
            for row in matrix:
                for value in row:
                    matrix_string += str(value)
            temporary.append(matrix_string)


schubert_cells = []


for string in temporary:
    n = len(string)
    matrix = [[0] * n for _ in range(n)]
   
    for col, char in enumerate(string):
        row = int(char) - 1
        matrix[row][col] = 1
   
    schubert_cells.append(matrix)

# After all schubert_cells are processed and symbolic 'a's added

print("Schubert Cells and Chi * C products:")

for C_index, matrix in enumerate(schubert_cells):
    print(f"\nMatrix {C_index + 1}:")
    for row in matrix:
        print(row)
    
    violations = []

    # Find positions where free symbolic values are needed
    for row_idx in range(len(matrix)):
        for col_idx in range(len(matrix[row_idx])):
            if matrix[row_idx][col_idx] == 0:
                has_left = any(matrix[row_idx][c] == 1 for c in range(col_idx))
                has_above = any(matrix[r][col_idx] == 1 for r in range(row_idx))
                if not has_left and not has_above:
                    violations.append((row_idx, col_idx))

    # Fill in symbolic 'a' values
    symbol_dict = {}
    for idx, (row_idx, col_idx) in enumerate(violations):
        symbol_name = f"a{idx + 1}"
        matrix[row_idx][col_idx] = symbol_name
        symbol_dict[symbol_name] = sp.symbols(symbol_name)

    if violations:
        print(f"\nMatrix {C_index + 1} violations at positions: {violations}")
        print("Modified matrix:")
        for row in matrix:
            print(row)
    else:
        print(f"\nMatrix {C_index + 1}: No violations")
        print("Modified matrix (no changes):")
        for row in matrix:
            print(row)


