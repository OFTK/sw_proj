#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

const static int MAX_ITER = 300;

/**
 * Parses a python list into a C type matrix.
 */
void py_arr_to_mtx(PyObject* py_arr, 
                  F_TYPE** arr_ptr, F_TYPE* arr_mem,
                  int i_max, int j_max)
{
    int i = 0;
    int j = 0;

    PyObject* double_obj;

    for (i = 0; i < i_max; i++)
    {
        arr_ptr[i] = arr_mem + (i * j_max);

        for (j = 0; j < j_max; j++)
        {
            double_obj = PyList_GetItem(py_arr, i*j_max+j);
            arr_mem[(i*j_max) + j] = PyFloat_AsDouble(double_obj);
        }
    }
}

/**
 * Executes the kmeans algorithm requested in the final exercise
 *  (Python section).
 * 
 * Should receive from python:
 * - dim (int), The dimension number of each point in the dataframe.
 * - k (int), The number of centroids received.
 * - centroids_list_obj (arr of type double), The calculated initial centroids.
 * - points_list_obj (arr of type double), The received dataframe (points).
 * 
 * Returns nothing, prints the centroids found.
 */
static PyObject* fit(PyObject* self, PyObject* args)
{
    int number_of_datapoints = 0;

    int dim = 0, k = 0;
    int i = 0, j = 0;

    int* output_cluster_assign = NULL;

    PyObject* centroids_list_obj = NULL;
    PyObject* points_list_obj = NULL;
    PyObject* float_obj = NULL;

    F_TYPE** centroids_arr_ptr = NULL;
    F_TYPE* centroids_arr_mem = NULL;
    F_TYPE** datapoints_arr_ptr = NULL;
    F_TYPE* datapoints_arr_mem = NULL;

    enum status status;

    if (!PyArg_ParseTuple(args, "iiOO", &dim, &k, &centroids_list_obj, 
                &points_list_obj))
    {
        PRINT_ERROR();
        return Py_BuildValue("");
    }

    number_of_datapoints = PyList_Size(points_list_obj) / dim;

    // Allocating according to the received length
    centroids_arr_ptr = calloc(sizeof(F_TYPE*), k);
    centroids_arr_mem = calloc(sizeof(F_TYPE), k * dim);
    datapoints_arr_ptr = calloc(sizeof(F_TYPE*), number_of_datapoints);
    datapoints_arr_mem = calloc(sizeof(F_TYPE),number_of_datapoints * dim);
    output_cluster_assign = calloc(sizeof(int), number_of_datapoints);
    if (NULL == datapoints_arr_mem || NULL == datapoints_arr_ptr || 
        NULL == centroids_arr_mem || NULL == centroids_arr_ptr ||
        NULL == output_cluster_assign) {
        free(datapoints_arr_mem); free(datapoints_arr_ptr); 
        free(centroids_arr_mem); free(centroids_arr_ptr);
        free(output_cluster_assign);
		PRINT_ERROR();
        return Py_BuildValue("");
	}

    /* Parsing initialized centroids... */
    py_arr_to_mtx(centroids_list_obj, centroids_arr_ptr, centroids_arr_mem,
                k, dim);

    // Parsing datapoints...
    py_arr_to_mtx(points_list_obj, datapoints_arr_ptr, datapoints_arr_mem, 
                number_of_datapoints, dim);

    /* THE KMEANS PROCEDURE CALL */
	status = kmeans(
		datapoints_arr_ptr, number_of_datapoints, dim, 
		centroids_arr_ptr, output_cluster_assign,
		k, MAX_ITER);

	if (Success == status) 
		print_matrix(centroids_arr_ptr, k, dim);
	else
		PRINT_ERROR();

    free(centroids_arr_ptr);
    free(centroids_arr_mem);
    free(datapoints_arr_ptr);
    free(datapoints_arr_mem);
    free(output_cluster_assign);
    return Py_BuildValue("");
}

/**
 * Executes an spkmeans algorithm goal according to a received 'goal' 
 * parameter.
 * 
 * Should receive from python:
 * - dim (int), The dimension number of each point in the dataframe.
 * - k (int), The number of centroids received.
 * - goal (int), The requested goal to execute [spk, wam, ddg, lnorm, jacobi].
 * - points_list_obj (arr of type double), The received dataframe (points).
 * 
 * Returns according to the received goal.
 *  If goal equals spk, returns the calculated T matrix as
 *  described in the specification document. 
 *  Otherwise, returns nothing and prints according to the goal's
 *  request.
 */
static PyObject* perform_subtask(PyObject* self, PyObject* args)
{
    int number_of_datapoints = 0;

    int dim = 0, k = 0, goal = 0;
    int i = 0, j = 0;

    PyObject* points_list_obj = NULL;
    PyObject* float_obj = NULL;
    PyObject* num = NULL;
    PyObject* return_value = PyList_New(0);

    F_TYPE** datapoints_arr_ptr = NULL;
    F_TYPE* datapoints_arr_mem = NULL;
    
    F_TYPE** o_mtx = NULL;
    F_TYPE* o_mtx_mem = NULL;

    F_TYPE* eigenvalues = NULL;

    enum status status;

    if (!PyArg_ParseTuple(args, "iiiO", &dim, &k, &goal, &points_list_obj))
        return return_value;

    number_of_datapoints = PyList_Size(points_list_obj) / dim;

    /* Allocating according to the received length */
    datapoints_arr_ptr = calloc(number_of_datapoints, sizeof(F_TYPE*));
    datapoints_arr_mem = calloc(number_of_datapoints * dim, sizeof(F_TYPE));
    if ((datapoints_arr_ptr == NULL) || (datapoints_arr_mem == NULL))
    {
        PRINT_ERROR();
        free(datapoints_arr_ptr); free(datapoints_arr_mem);
        return return_value;
    }

    /* Parsing datapoints... */
    py_arr_to_mtx(points_list_obj, datapoints_arr_ptr, datapoints_arr_mem, 
                number_of_datapoints, dim);

    /* Allocate memory */
	o_mtx = calloc(number_of_datapoints, sizeof(F_TYPE*));
	o_mtx_mem = calloc(number_of_datapoints*number_of_datapoints, 
            sizeof(F_TYPE));
    eigenvalues = calloc(number_of_datapoints, sizeof(F_TYPE));
	if ((NULL == o_mtx) || (NULL == o_mtx_mem) || (NULL == eigenvalues)) {
		PRINT_ERROR();
		goto finish;
	}

	for (i = 0; i < number_of_datapoints; ++i)
		o_mtx[i] = o_mtx_mem + (i*number_of_datapoints);

    /* Reseting errno */
	errno = 0;

	/* Performing subtask */
    status = spkmeans_preperations(
    datapoints_arr_ptr, number_of_datapoints, dim, &k, goal,
    o_mtx_mem, eigenvalues);

	if (Error == status) {
		PRINT_ERROR();
        goto finish;
	}

	/* Print output matrix (and if needed - eigenvalues) */
	if (SPK != goal) {
		if (JACOBI == goal) 
            print_matrix(&eigenvalues, 1, number_of_datapoints);
		
        print_matrix(o_mtx, number_of_datapoints, number_of_datapoints);
	}
    else /* Returning the received value to the python module */
    {
        Py_DECREF(return_value); /* Because we've already assigned a 
                                    py object... */
        return_value = PyList_New(k * number_of_datapoints);

    /* If NULL */
    if (!return_value)
    {
        Py_DECREF(return_value);
        PRINT_ERROR();
        goto finish;
    }

    for (i = 0; i < k * number_of_datapoints; i++)
    {
        num = PyFloat_FromDouble(o_mtx_mem[i]);

        if (!num) {
            Py_DECREF(return_value);
            goto finish;
        }

        PyList_SET_ITEM(return_value, i, num);
        }
    }

finish:

    free(datapoints_arr_ptr);
    free(datapoints_arr_mem);

    free(o_mtx); 
    free(o_mtx_mem); 
    free(eigenvalues);

    return return_value;
}

static PyMethodDef methods[] = {
        {"perform_subtask", perform_subtask, METH_VARARGS,
        "Executes a sub-spkmeans function based on a received 'goal' 
        parameter"},
        {"fit", fit, METH_VARARGS,
        "Executes kmeans algorithm"},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
        PyModuleDef_HEAD_INIT,
        "myspkmeanssp",
        "Performs a chosen function from [spk, wam, ddg, lnorm, jacobi]",
        -1,
        methods
};

PyMODINIT_FUNC PyInit_myspkmeanssp(void)
{
    return PyModule_Create(&module);
}