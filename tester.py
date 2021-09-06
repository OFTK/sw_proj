import sys
import numpy as np
from scipy.sparse.csgraph import laplacian
import math

global EPSILON
EPSILON = 10^(-6)


def read_inputs(fp):
    datapoints = []

    f = open(fp, 'r')

    # Reading the first input separately in order to set some needed parameters such as dimension
    first_line = f.readline().split(',')
    dim = len(first_line)

    new_point = []
    for val in first_line:
        try:
            new_point.append(float(val))
        except(ValueError):
            print(f"first datapoint's value wasn't able to be parsed as an int, exiting...")
            return ERROR

    datapoints.append(new_point)

    for line in f.readlines():
        splitted_line = line.split(',')

        if len(splitted_line) != dim:
            print(f"All datapoints should be of dimension {dim}, received a bad input of dim {len(splitted_line)}, exiting...")
            return ERROR
        
        new_point = []
        for val in splitted_line:
            try:
                new_point.append(float(val))
            except(ValueError):
                print(f"some datapoints' value wasn't able to be parsed as an int, exiting...")
                return ERROR

        datapoints.append(new_point)

    return datapoints


def create_wam(dps):
	n = len(dps)
	datapoints = np.array(dps)

	W = np.full((n, n), 0, dtype=np.float64)
	for i in range(n):
		W[i] = np.exp(-(np.linalg.norm(datapoints - datapoints[i], axis=1)) / 2)
	np.fill_diagonal(W, val=0)
	return W

def create_ddg(wam):
	n = len(wam)
	ddg = np.eye(n)
	for i in range(n):
		ddg[i][i] = (np.sum(wam[i])) ** (-.5)
	return ddg

def create_laplacian(wam):
	return laplacian(wam, True)

def print_matrix(matrix):
    for dp in matrix:
        for j in range(len(dp) - 1):
            print('{:.3f},'.format(dp[j]), end='')
        print('{:.3f}'.format(dp[-1]))

def build_P_matrix(N, i, j, c, s):
    result = np.eye(N)
    result[i][i] = c
    result[j][j] = c
    result[i][j] = s
    result[j][i] = -s
    return result


def calc_off(matrix):
    result = 0
    for i in range(len(matrix)):
        for j in range(i+1, len(matrix)):
            result += matrix[i][j]**2 + matrix[j][i]**2
    return result


def find_off_diag_max(matrix):
    matrix_max = -1
    max_i = -1
    max_j = -1
    for i in range(len(matrix)):
        for j in range(i+1, len(matrix)):
            if abs(matrix[i][j]) > matrix_max:
                matrix_max = abs(matrix[i][j])
                max_i = i
                max_j = j
    return max_i, max_j


def calc_norm(x1, x2):
    return np.sqrt(np.sum((x1 - x2)**2, axis=2))


def check_equality(matrix1, matrix2):
    if len(matrix1) != len(matrix2):
        return False
    for i in range(len(matrix1)):
        if len(matrix1[i]) != len(matrix2[i]):
            return False
        for j in range(len(matrix1[0])):
            if abs(matrix1[i][j] - matrix2[i][j]) > EPSILON:
                return False
    return True


def build_jacobi_matrix(lmatrix):
    jacobi_values = lmatrix
    n_size = len(jacobi_values)
    jacobi_vectors = np.eye(n_size)
    cur_f_norm = calc_off(jacobi_values)
    for iter in range(100):
        i, j = find_off_diag_max(jacobi_values)
        psi = (jacobi_values[j][j] - jacobi_values[i][i]) / (2 * jacobi_values[i][j])
        t = (-1 if psi < 0 else 1) / (abs(psi) + math.sqrt(psi ** 2 + 1))
        c = 1 / math.sqrt(1 + t ** 2)
        s = t * c
        P_matrix = build_P_matrix(n_size,i,j,c,s)
        jacobi_values = P_matrix.transpose() @ jacobi_values @ P_matrix
        jacobi_vectors = jacobi_vectors @ P_matrix
        prev_f_norm = cur_f_norm
        cur_f_norm = calc_off(jacobi_values)
        if abs(cur_f_norm - prev_f_norm) < EPSILON:
            break
    return np.ndarray.tolist(jacobi_values), np.ndarray.tolist(jacobi_vectors)




if __name__ == '__main__':
	datapoints = read_inputs("input.csv") 
	# TODO: generate inputs
	print_matrix(datapoints)
	print("")

	wam = create_wam(datapoints)
	print_matrix(wam)
	print("")

	ddg = create_ddg(wam)
	print_matrix(ddg)
	print("")

	lnorm = create_laplacian(wam)
	print_matrix(lnorm)
	print("")

	eigenvalues, eigenvectors = build_jacobi_matrix(lnorm)
	print_matrix(eigenvalues)
	print("")
	print_matrix(eigenvectors)
	print("")
