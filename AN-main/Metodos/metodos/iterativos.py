import numpy as np
import pandas as pd

def metodo_jacobi(a,x0,b,tol,niter):
    respuesta=MatJacobiSeid(x0,a,b,tol,niter,0)
    if respuesta['solucion']:
        respuesta['df']=respuesta['tabla'].to_html()
        respuesta['nombre_metodo']='Jacobi'
        respuesta['T']=respuesta['T'].to_html()
        respuesta['C']=respuesta['C'].to_html()
    return respuesta


def metodo_gauss_seidel(a,x0,b,tol,niter):
    respuesta=MatJacobiSeid(x0,a,b,tol,niter,1)
    if respuesta['solucion']:
        respuesta['df']=respuesta['tabla'].to_html()
        respuesta['nombre_metodo']='Gauss-Seidel'
        respuesta['T']=respuesta['T'].to_html()
        respuesta['C']=respuesta['C'].to_html()
    return respuesta



def MatJacobiSeid(x0, A, b, Tol, niter, met):
    c = 0
    error = Tol + 1
    D = np.diag(np.diag(A))
    L = -np.tril(A, -1)
    U = -np.triu(A, +1)

    tabla = { 'Error': [], 'Vector xi': []}
    respuesta={'solucion':False}
    #verificamos que D tenga inversa
    if np.linalg.det(D) == 0:
        mensaje="The matrix is ​​not invertible"
        respuesta['mensaje']= mensaje
        return respuesta
    tabla['Error'].append(0)
    tabla['Vector xi'].append(x0)
    while error > Tol and c < niter:
        if met == 0:
            T = np.linalg.inv(D) @ (L + U)
            C = np.linalg.inv(D) @ b
            x1 = T @ x0 + C
        elif met == 1:
            T = np.linalg.inv(D - L) @ U
            C = np.linalg.inv(D - L) @ b
            x1 = T @ x0 + C
        E = np.linalg.norm(x1 - x0, np.inf)
        tabla['Error'].append(E)
        tabla['Vector xi'].append(x1)
        error = E
        x0 = x1
        c += 1
    respuesta['T']=pd.DataFrame(T)
    respuesta['C']=pd.DataFrame(C)
    respuesta['radio_esp']=radio_espectral(T)
    respuesta['tabla'] = pd.DataFrame(tabla)
    if error < Tol:
        s = x0
        n = c
        mensaje='An approximation of the system solution that meets tolerance was achieved ='+str(Tol)
        respuesta['solucion']=True
        respuesta['mensaje']=mensaje
    else:
        s = x0
        n = c
        mensaje='Failed in '+str(niter)+' iterations'
        respuesta['mensaje']=mensaje
    return respuesta

def radio_espectral(matriz):
    eigenvalores = np.linalg.eigvals(matriz)
    radio_espectral = np.abs(max(eigenvalores, key=abs))
    return radio_espectral

def metodo_sor(A,x0,b,Tol,niter,w):
    c = 0
    error = Tol + 1
    D = np.diag(np.diag(A))
    L = -np.tril(A, -1)
    U = -np.triu(A, +1)
    tabla = { 'Error': [], 'Vector xi': []}
    respuesta={'solucion':False}
    if np.linalg.det(D) == 0:
        mensaje="The matrix is ​​not invertible"
        respuesta['mensaje']= mensaje
        return respuesta
    tabla['Error'].append(0)
    tabla['Vector xi'].append(x0)
    while error > Tol and c < niter:
        T = np.linalg.inv(D - w * L) @ ((1 - w) * D + w * U)
        C = w * np.linalg.inv(D - w * L) @ b
        x1 = T @ x0 + C
        E = np.linalg.norm(x1 - x0, np.inf)
        tabla['Error'].append(E)
        tabla['Vector xi'].append(x1)
        error = E
        x0 = x1
        c += 1
    print("salimos del while ")
    respuesta['T']=pd.DataFrame(T)
    respuesta['C']=pd.DataFrame(C)
    respuesta['radio_esp']=radio_espectral(T)
    respuesta['tabla'] = pd.DataFrame(tabla)
    if error < Tol:
        s = x0
        n = c
        mensaje='An approximation of the system solution that meets tolerance was achieved ='+str(Tol)
        respuesta['solucion']=True
        respuesta['mensaje']=mensaje
    else:
        s = x0
        n = c
        mensaje='Failed in '+str(niter)+' iterations'
        respuesta['mensaje']=mensaje
    print("let's try to see if solution")
    
    if respuesta['solucion']:
        respuesta['df']=respuesta['tabla'].to_html()
        respuesta['nombre_metodo']='SOR'
        respuesta['T']=respuesta['T'].to_html()
        respuesta['C']=respuesta['C'].to_html()
    print("if ready")
    return respuesta


def eliminacion_gaussiana_simple(A, b):
    """
    Realiza la eliminación gaussiana simple para resolver el sistema de ecuaciones Ax = b.
    """
    n = len(A)
    A = A.astype(float)
    b = b.astype(float)
    respuesta = {'solucion': False}
    for i in range(n - 1):
        if A[i, i] == 0:
            mensaje = "No es posible continuar, elemento diagonal nulo."
            respuesta['mensaje'] = mensaje
            return respuesta
        for j in range(i + 1, n):
            m = A[j, i] / A[i, i]
            A[j, i:] = A[j, i:] - m * A[i, i:]
            b[j] = b[j] - m * b[i]
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i]
    respuesta['solucion'] = True
    respuesta['mensaje'] = "El sistema fue resuelto exitosamente con eliminación gaussiana simple."
    respuesta['x'] = x
    respuesta['A_final'] = pd.DataFrame(A)
    respuesta['b_final'] = pd.DataFrame(b, columns=['b'])
    return respuesta


def eliminacion_gaussiana_pivoteo_parcial(A, b):
    """
    Realiza la eliminación gaussiana con pivoteo parcial para resolver el sistema de ecuaciones Ax = b.
    """
    n = len(A)
    A = A.astype(float)
    b = b.astype(float)
    respuesta = {'solucion': False}
    for i in range(n - 1):
        # Pivoteo parcial
        max_index = np.argmax(np.abs(A[i:, i])) + i
        if A[max_index, i] == 0:
            mensaje = "No es posible continuar, matriz singular detectada."
            respuesta['mensaje'] = mensaje
            return respuesta
        if max_index != i:
            # Intercambiar filas
            A[[i, max_index]] = A[[max_index, i]]
            b[[i, max_index]] = b[[max_index, i]]
        for j in range(i + 1, n):
            m = A[j, i] / A[i, i]
            A[j, i:] = A[j, i:] - m * A[i, i:]
            b[j] = b[j] - m * b[i]
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i]
    respuesta['solucion'] = True
    respuesta['mensaje'] = "El sistema fue resuelto exitosamente con pivoteo parcial."
    respuesta['x'] = x
    respuesta['A_final'] = pd.DataFrame(A)
    respuesta['b_final'] = pd.DataFrame(b, columns=['b'])
    return respuesta


def eliminacion_gaussiana_pivoteo_total(A, b):
    """
    Realiza la eliminación gaussiana con pivoteo total para resolver el sistema de ecuaciones Ax = b.
    """
    n = len(A)
    A = A.astype(float)
    b = b.astype(float)
    indices = np.arange(n)
    respuesta = {'solucion': False}
    for i in range(n - 1):
        # Pivoteo total
        max_index = np.unravel_index(np.argmax(np.abs(A[i:, i:])), A[i:, i:].shape)
        max_index = (max_index[0] + i, max_index[1] + i)
        if A[max_index] == 0:
            mensaje = "No es posible continuar, matriz singular detectada."
            respuesta['mensaje'] = mensaje
            return respuesta
        if max_index[0] != i:
            # Intercambiar filas
            A[[i, max_index[0]]] = A[[max_index[0], i]]
            b[[i, max_index[0]]] = b[[max_index[0], i]]
        if max_index[1] != i:
            # Intercambiar columnas
            A[:, [i, max_index[1]]] = A[:, [max_index[1], i]]
            indices[[i, max_index[1]]] = indices[[max_index[1], i]]
        for j in range(i + 1, n):
            m = A[j, i] / A[i, i]
            A[j, i:] = A[j, i:] - m * A[i, i:]
            b[j] = b[j] - m * b[i]
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i]
    # Reorganizar las soluciones según los cambios de columnas
    x_final = np.zeros_like(x)
    x_final[indices] = x
    respuesta['solucion'] = True
    respuesta['mensaje'] = "El sistema fue resuelto exitosamente con pivoteo total."
    respuesta['x'] = x_final
    respuesta['A_final'] = pd.DataFrame(A)
    respuesta['b_final'] = pd.DataFrame(b, columns=['b'])
    return respuesta
