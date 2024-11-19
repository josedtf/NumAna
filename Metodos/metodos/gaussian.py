import numpy as np
import pandas as pd

def eliminacion_gaussiana_simple(A, b):
    try:
        n = len(A)
        # Forward elimination
        for k in range(n - 1):
            for i in range(k + 1, n):
                if A[k][k] == 0:
                    raise ValueError("Zero pivot encountered in Gaussian Elimination.")
                factor = A[i][k] / A[k][k]
                for j in range(k, n):
                    A[i][j] -= factor * A[k][j]
                b[i] -= factor * b[k]

        # Back substitution
        x = [0] * n
        for i in range(n - 1, -1, -1):
            x[i] = (b[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))) / A[i][i]

        return {
            'solucion': True,
            'triangular': pd.DataFrame(A).to_html(),
            'solucion_tabla': pd.DataFrame(x, columns=["Solution"]).to_html(),
            'mensaje': 'Gaussian elimination completed successfully.'
        }
    except Exception as e:
        return {
            'solucion': False,
            'mensaje': f'Error during Gaussian Elimination: {str(e)}'
        }

def eliminacion_gaussiana_pivoteo_parcial(A, b):
    try:
        n = len(A)
        for k in range(n - 1):
            # Partial pivoting
            max_row = max(range(k, n), key=lambda i: abs(A[i][k]))
            if A[max_row][k] == 0:
                raise ValueError("Matrix is singular.")
            if max_row != k:
                A[[k, max_row]] = A[[max_row, k]]
                b[k], b[max_row] = b[max_row], b[k]

            # Forward elimination
            for i in range(k + 1, n):
                factor = A[i][k] / A[k][k]
                for j in range(k, n):
                    A[i][j] -= factor * A[k][j]
                b[i] -= factor * b[k]

        # Back substitution
        x = [0] * n
        for i in range(n - 1, -1, -1):
            x[i] = (b[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))) / A[i][i]

        return {
            'solucion': True,
            'triangular': pd.DataFrame(A).to_html(),
            'solucion_tabla': pd.DataFrame(x, columns=["Solution"]).to_html(),
            'mensaje': 'Gaussian elimination with partial pivoting completed successfully.'
        }
    except Exception as e:
        return {
            'solucion': False,
            'mensaje': f'Error during Gaussian Elimination with Partial Pivoting: {str(e)}'
        }

def eliminacion_gaussiana_pivoteo_total(A, b):
    try:
        n = len(A)
        indices = list(range(n))
        for k in range(n - 1):
            # Total pivoting
            max_index = max(((i, j) for i in range(k, n) for j in range(k, n)), key=lambda x: abs(A[x[0]][x[1]]))
            if A[max_index[0]][max_index[1]] == 0:
                raise ValueError("Matrix is singular.")
            i_max, j_max = max_index
            if i_max != k:
                A[[k, i_max]] = A[[i_max, k]]
                b[k], b[i_max] = b[i_max], b[k]
            if j_max != k:
                A[:, [k, j_max]] = A[:, [j_max, k]]
                indices[k], indices[j_max] = indices[j_max], indices[k]

            # Forward elimination
            for i in range(k + 1, n):
                factor = A[i][k] / A[k][k]
                for j in range(k, n):
                    A[i][j] -= factor * A[k][j]
                b[i] -= factor * b[k]

        # Back substitution
        x = [0] * n
        for i in range(n - 1, -1, -1):
            x[i] = (b[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))) / A[i][i]

        # Rearrange x to match original variable order
        x_final = [0] * n
        for i, index in enumerate(indices):
            x_final[index] = x[i]

        return {
            'solucion': True,
            'triangular': pd.DataFrame(A).to_html(),
            'solucion_tabla': pd.DataFrame(x_final, columns=["Solution"]).to_html(),
            'mensaje': 'Gaussian elimination with total pivoting completed successfully.'
        }
    except Exception as e:
        return {
            'solucion': False,
            'mensaje': f'Error during Gaussian Elimination with Total Pivoting: {str(e)}'
        }
