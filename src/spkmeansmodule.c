#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

const static int MAX_ITER = 300;

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

#ifdef DEBUG
    printf("Got to here! In fit\n");
#endif

    if (!PyArg_ParseTuple(args, "iiOO", &dim, &k, &centroids_list_obj, &points_list_obj))
    {
        PRINT_ERROR();
        return Py_BuildValue("");
    }

#ifdef DEBUG
    printf("Passed parseargs!\n");

    printf("d:%d\n", dim);
    printf("k:%d\n", k);
#endif

    number_of_datapoints = PyList_Size(points_list_obj) / dim;

#ifdef DEBUG
    printf("n:%d\n", number_of_datapoints);
#endif

    // Allocating according to the received length
    centroids_arr_ptr = calloc(sizeof(F_TYPE*), k);
    centroids_arr_mem = calloc(sizeof(F_TYPE), k * dim);
    datapoints_arr_ptr = calloc(sizeof(F_TYPE*), number_of_datapoints);
    datapoints_arr_mem = calloc(sizeof(F_TYPE),number_of_datapoints * dim);
    output_cluster_assign = calloc(sizeof(int), dim);
    if (NULL == datapoints_arr_mem || NULL == datapoints_arr_ptr || 
        NULL == centroids_arr_mem || NULL == centroids_arr_ptr ||
        NULL == output_cluster_assign) {
        free(datapoints_arr_mem); free(datapoints_arr_ptr); 
        free(centroids_arr_mem); free(centroids_arr_ptr);
        free(output_cluster_assign);
		PRINT_ERROR();
        Py_BuildValue("");
	}

#ifdef DEBUG
    printf("Start parsing centroids!\n");
#endif

    // Parsing initialized centroids...
    for (i = 0; i < k; i++)
    {
        centroids_arr_ptr[i] = centroids_arr_mem + (i * dim);

        for (j = 0; j < dim; j++)
        {
            float_obj = PyList_GetItem(centroids_list_obj, i*dim+j);
            centroids_arr_mem[(i*dim) + j] = PyFloat_AsDouble(float_obj);
#ifdef DEBUG
            printf("centroid: %f\n", centroids_arr_mem[(i*dim) + j]);
#endif
        }
    }

#ifdef DEBUG
    printf("Start parsing dps!\n");
#endif

    // Parsing datapoints...
    for (i = 0; i <number_of_datapoints; i++)
    {
        datapoints_arr_ptr[i] = datapoints_arr_mem + (i * dim);

        for (j = 0; j < dim; j++)
        {
            float_obj = PyList_GetItem(points_list_obj, i*dim+j);
            datapoints_arr_mem[(i*dim) + j] = PyFloat_AsDouble(float_obj);
#ifdef DEBUG
            printf("dp: %f\n", datapoints_arr_mem[(i*dim) + j]);
#endif
        }
    }

    /* THE KMEANS PROCEDURE CALL */
	status = kmeans(
		datapoints_arr_ptr, dim, k, 
		centroids_arr_ptr, output_cluster_assign,
		k, MAX_ITER);

	if (Success == status) 
		print_matrix(centroids_arr_ptr, k, dim);
	else
		PRINT_ERROR();

#ifdef DEBUG
    printf("Finished! :D\n");
#endif

    free(centroids_arr_ptr);
    free(centroids_arr_mem);
    free(datapoints_arr_ptr);
    free(datapoints_arr_mem);
    free(output_cluster_assign);
    Py_BuildValue("");
}

static PyObject* perform_subtask(PyObject* self, PyObject* args)
{
    int number_of_datapoints = 0;

    int dim = 0, k = 0, goal = 0;
    int i = 0, j = 0;

    PyObject* points_list_obj = NULL;
    PyObject* float_obj = NULL;
    PyObject* return_value = PyList_New(0);
    PyObject* num = NULL;

    F_TYPE** datapoints_arr_ptr = NULL;
    F_TYPE* datapoints_arr_mem = NULL;
    
    F_TYPE** o_mtx = NULL;
    F_TYPE* o_mtx_mem = NULL;

    F_TYPE* eigenvalues = NULL;

    enum status status;

#ifdef DEBUG
    printf("Got to here! in subtask\n");
#endif

    if (!PyArg_ParseTuple(args, "iiiO", &dim, &k, &goal, &points_list_obj))
        return return_value;

#ifdef DEBUG
    printf("Passed parseargs!\n");

    printf("dim:%d\n", dim);
    printf("k:%d\n", k);
    printf("goal:%d\n", goal);
#endif

    number_of_datapoints = PyList_Size(points_list_obj) / dim;

#ifdef DEBUG
    printf("n:%d\n", number_of_datapoints);
#endif

    // Allocating according to the received length
    datapoints_arr_ptr = calloc(number_of_datapoints, sizeof(F_TYPE*));
    datapoints_arr_mem = calloc(number_of_datapoints * dim, sizeof(F_TYPE));
    if ((datapoints_arr_ptr == NULL) || (datapoints_arr_mem == NULL))
    {
        PRINT_ERROR();
        free(datapoints_arr_ptr); free(datapoints_arr_mem);
        return return_value;
    }

#ifdef DEBUG
    printf("Start parsing dps!\n");
#endif

    // Parsing datapoints...
    for (i = 0; i <number_of_datapoints; i++)
    {
        datapoints_arr_ptr[i] = datapoints_arr_mem + (i * dim);

        for (j = 0; j < dim; j++)
        {
            float_obj = PyList_GetItem(points_list_obj, i*dim+j);
            datapoints_arr_mem[(i*dim) + j] = PyFloat_AsDouble(float_obj);
        }
    }

    #ifdef DEBUG
            printf("Parsed dps:\n");
            print_matrix(datapoints_arr_ptr, number_of_datapoints, dim);
    #endif

    /* Allocate memory */
	o_mtx = calloc(number_of_datapoints, sizeof(F_TYPE*));
	o_mtx_mem = calloc(number_of_datapoints*number_of_datapoints, sizeof(F_TYPE));
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

    #ifdef DEBUG
        printf("Finished prep!\n");
    #endif

	/* Print output matrix (and if needed - eigenvalues) */
	if (SPK != goal) {
		if (JACOBI == goal) print_matrix(&eigenvalues, 1, number_of_datapoints);
		print_matrix(o_mtx, number_of_datapoints, number_of_datapoints);
	}
    else /* Returning the received value to the python module */
    {   
        Py_DECREF(return_value); // Because we've already assigned a py object...
        return_value = PyList_New(k * number_of_datapoints);

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
        "Executes a sub-spkmeans function based on a received 'goal' parameter"},
        {"fit", fit, METH_VARARGS,
        "Executes the spkmeans algorithm given"},
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