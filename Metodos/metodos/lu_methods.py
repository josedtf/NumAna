import numpy as np
import pandas as pd

def metodo_lu_simple(A, b):
    try:
        # LU Decomposition
        n = len(A)
        L = np.zeros((n, n))
        U = np.zeros((n, n))

        for i in range(n):
            L[i][i] = 1
            for j in range(i, n):
                U[i][j] = A[i][j] - sum(L[i][k] * U[k][j] for k in range(i))
            for j in range(i + 1, n):
                L[j][i] = (A[j][i] - sum(L[j][k] * U[k][i] for k in range(i))) / U[i][i]

        # Solve Ly = b
        y = np.zeros(n)
        for i in range(n):
            y[i] = b[i] - sum(L[i][k] * y[k] for k in range(i))

        # Solve Ux = y
        x = np.zeros(n)
        for i in range(n - 1, -1, -1):
            x[i] = (y[i] - sum(U[i][k] * x[k] for k in range(i + 1, n))) / U[i][i]

        return {
            'solucion': True,
            'L': pd.DataFrame(L).to_html(),
            'U': pd.DataFrame(U).to_html(),
            'solucion_tabla': pd.DataFrame(x, columns=["Solution"]).to_html(),
            'P': None,
            'mensaje': 'LU decomposition was successful.'
        }
    except Exception as e:
        return {
            'solucion': False,
            'mensaje': f'Error during LU Simple calculation: {str(e)}'
        }

def metodo_lu_pivoting(A, b):
    try:
        # LU Decomposition with Pivoting
        n = len(A)
        P = np.eye(n)
        L = np.zeros((n, n))
        U = np.copy(A)

        for i in range(n):
            # Pivoting
            max_row = np.argmax(abs(U[i:, i])) + i
            if i != max_row:
                U[[i, max_row]] = U[[max_row, i]]
                P[[i, max_row]] = P[[max_row, i]]
                if i > 0:
                    L[[i, max_row], :i] = L[[max_row, i], :i]

            L[i][i] = 1
            for j in range(i + 1, n):
                L[j][i] = U[j][i] / U[i][i]
                U[j, i:] -= L[j][i] * U[i, i:]

        # Solve Ly = Pb
        Pb = P @ b
        y = np.zeros(n)
        for i in range(n):
            y[i] = Pb[i] - sum(L[i][k] * y[k] for k in range(i))

        # Solve Ux = y
        x = np.zeros(n)
        for i in range(n - 1, -1, -1):
            x[i] = (y[i] - sum(U[i][k] * x[k] for k in range(i + 1, n))) / U[i][i]

        return {
            'solucion': True,
            'L': pd.DataFrame(L).to_html(),
            'U': pd.DataFrame(U).to_html(),
            'P': pd.DataFrame(P).to_html(),
            'solucion_tabla': pd.DataFrame(x, columns=["Solution"]).to_html(),
            'mensaje': 'LU decomposition with pivoting was successful.'
        }
    except Exception as e:
        return {
            'solucion': False,
            'mensaje': f'Error during LU Pivoting calculation: {str(e)}'
        }
