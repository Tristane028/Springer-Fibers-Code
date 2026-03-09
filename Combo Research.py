from itertools import permutations
import numpy as np
import sympy as sp # For symbolic computations (used in Mathematica)

#The following line is for running the script on my computer specifically, omit this line when running this
#/usr/local/bin/python3.9 "/Users/christan065/Springer Fibers/Combo Research.py"

# Young Tableau: reversed digits, each row has d columns.
def young_tableaux_from_digits(num_str):
    return [list(range(1, int(d)+1)) for d in reversed(num_str)]

# Nilpotent Jordan matrix with blocks given by digits.
def jordan_from_digits(num_str):
    digits = [int(d) for d in num_str]
    n = sum(digits)
    J = np.zeros((n, n), dtype=int)
    idx = 0
    for k in digits:
        for i in range(k-1):
            J[idx+i, idx+i+1] = 1
        idx += k
    return sp.Matrix(J)

# Check if a matrix is a valid Young tableau (columns strictly increase and no duplicates in columns)
def is_valid_matrix(matrix):
    for col in range(len(matrix[0])):
        col_vals = [matrix[r][col] for r in range(len(matrix)) if col < len(matrix[r])]
        if sorted(col_vals) != col_vals or len(col_vals) != len(set(col_vals)):
            return False
    return True

# Verify columns are strictly increasing (for Schubert cells, we want columns to be strictly increasing)
def is_column_increasing(matrix):
    max_cols = max(len(r) for r in matrix)
    for c in range(max_cols):
        for r in range(1, len(matrix)):
            if c < len(matrix[r]) and c < len(matrix[r-1]):
                if matrix[r][c] <= matrix[r-1][c]:
                    return False
    return True

# Counts number of inversions (if b is below and to the left of a, and b<a, that's an inversion)
def count_inversions(matrix):
    pos = [(val, r, c) for r,row in enumerate(matrix) for c,val in enumerate(row)]
    inv = 0
    for i in range(len(pos)):
        for j in range(i+1, len(pos)):
            vi,ri,ci = pos[i]
            vj,rj,cj = pos[j]
            if ri==rj and ci<cj and vi>vj: inv+=1
            if rj>ri and cj<ci and vj<vi: inv+=1
    return inv

# Generate all valid fillings of the Young tableau shape.
def generate_valid_fillings(tableaux):
    n = sum(len(r) for r in tableaux)
    nums = list(range(1,n+1))
    valid = []
    for perm in permutations(nums):
        m = []
        idx = 0
        for row in tableaux:
            m.append(list(perm[idx:idx+len(row)]))
            idx += len(row)
        if is_valid_matrix(m) and is_column_increasing(m):
            valid.append(m)
    return valid

# Convert valid fillings to Schubert cell matrices (1 for positions of numbers, 0 elsewhere)
def fillings_to_schubert_cells(valid):
    cells = []
    for mat in valid:
        flat = [v for row in mat for v in row]
        n = len(flat)
        S = [[0]*n for _ in range(n)]
        for col,val in enumerate(flat):
            S[val-1][col] = 1
        cells.append(S)
    return cells

# Insert free a-variables not hit by the "death rays of 1"
def fill_symbolic_entries(cell):
    a_symbols = {}
    idx = 1
    M = [row[:] for row in cell]
    for r in range(len(M)):
        for c in range(len(M)):
            if M[r][c] == 0:
                has_left = any(M[r][cc]==1 for cc in range(c))
                has_above = any(M[rr][c]==1 for rr in range(r))
                if not has_left and not has_above:
                    name = f"a{idx}"
                    M[r][c] = name
                    a_symbols[name] = sp.symbols(name)
                    idx += 1
    return M, a_symbols

# Span Checking function for each column of XC against C (tracking new relations among a-variables).
def springer_span_checks(C_sym, XC_sym):
    a_vars = sorted(
        {s for s in C_sym.free_symbols if s.name.startswith("a")},
        key=lambda x: int(x.name[1:])
    )
    relations = {}
    results = []

    for i in range(C_sym.cols):
        C_sub = C_sym.subs(relations)
        XC_sub = XC_sym.subs(relations)
        C_cols = C_sub[:, :i+1]
        target = XC_sub[:, i]

        alphas = sp.symbols(f"alpha0:{i+1}")
        unknowns = list(alphas) + a_vars
        eqs = list(C_cols*sp.Matrix(alphas) - target)
        sol = sp.solve(eqs, unknowns, dict=True)

        if not sol:
            results.append(("NO", {}))
            continue

        s = sol[0]
        new_rel = {}
        for v in a_vars:
            if v in s:
                expr = s[v]
                if any(a in expr.free_symbols for a in alphas):
                    new_rel[v] = "free"
                else:
                    new_rel[v] = expr

        for k,v in new_rel.items():
            if v!="free":
                relations[k]=v

        results.append(("YES", new_rel))

    constrained = set(relations.keys())
    free_vars = [v for v in a_vars if v not in constrained]
    return results, free_vars, relations

num_str = input("Enter number string: ")
tableaux = young_tableaux_from_digits(num_str)
J = jordan_from_digits(num_str)

valid = generate_valid_fillings(tableaux)
cells = fillings_to_schubert_cells(valid)

print("\n=== Schubert Cells and Springer Span Checks ===")

for idx, cell in enumerate(cells, start=1):
    print(f"\n--- Cell {idx} ---")
    for row in cell: print(row)

    # Insert a-variables
    M, a_dict = fill_symbolic_entries(cell)

    print("\nWith symbolic entries:")
    for row in M: print(row)

    # Implement symbolic matrix multiplication to get Chi-C matrix
    C_sym = sp.Matrix([
        [a_dict.get(val,val) for val in row]
        for row in M
    ])
    XC_sym = J * C_sym

    # Zero out rows where J is zero
    for r in range(J.rows):
        if all(J[r,c]==0 for c in range(J.cols)):
            for c in range(XC_sym.cols):
                XC_sym[r,c] = 0

    print("\nChi * C:")
    sp.pprint(XC_sym)

    # Checks Span (whether columns of Chi-C can be written as a linear combination of all columns leading up to column position in Chi-C)
    checks, free_vars, relations = springer_span_checks(C_sym, XC_sym)

    print("\nSpringer Span Checks:")
    print("{:<8} {:<8} {:<30}".format("Check","Y/N","Relations"))
    print("-"*40)
    for i,(yn,new_rel) in enumerate(checks, start=1):
        if new_rel:
            rel = ", ".join(
                f"{k}=free" if v=="free" else f"{k}={v}"
                for k,v in new_rel.items()
            )
        else:
            rel = "None"
        print("{:<8} {:<8} {:<30}".format(i, yn, rel))

    print("\nFinal Free Vars:", ", ".join(str(v) for v in free_vars))
    print("Final Relations:", ", ".join(f"{k}={v}" for k,v in relations.items()))
