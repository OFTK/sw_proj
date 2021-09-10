#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

const static int MAX_ITER = 300;

static PyObject* spkmeans_fit(PyObject* self, PyObject* args)
{
    int number_of_datapoints;

    int dim, k;

    PyObject* centroids_list_obj = NULL;
    PyObject* points_list_obj = NULL;
    PyObject* float_obj = NULL;

    F_TYPE** centroids_arr_ptr = NULL;
    F_TYPE* centroids_arr_mem = NULL;
    F_TYPE** datapoints_arr_ptr = NULL;
    F_TYPE* datapoints_arr_mem = NULL;

#ifdef DEBUG
    printf("Got to here!\n");
#endif

    if (!PyArg_ParseTuple(args, "iiOO", &dim, &k, &centroids_list_obj, &points_list_obj))
        return NULL;

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
    assert(centroids_arr_ptr != NULL);
    assert(centroids_arr_mem != NULL);
    assert(datapoints_arr_ptr != NULL);
    assert(datapoints_arr_mem != NULL);

#ifdef DEBUG
    printf("Start parsing centroids!\n");
#endif

    // Parsing initialized centroids...
    for (int i = 0; i <K; i++)
    {
        centroids_arr_ptr[i] = centroids_arr_mem + (i * dim);

        for (int j = 0; j < dim; j++)
        {
            float_obj = PyList_GetItem(centroids_list_obj, i*dim+j);
            centroids_arr_mem[(i*dim) + j] = PyFloat_AsDouble(float_obj);
#ifdef DEBUG
            printf("centroid: %Lf\n", centroids_arr_mem[(i*dim) + j]);
#endif
        }
    }

#ifdef DEBUG
    printf("Start parsing dps!\n");
#endif

    // Parsing datapoints...
    for (int i = 0; i <number_of_datapoints; i++)
    {
        datapoints_arr_ptr[i] = datapoints_arr_mem + (i * dim);

        for (int j = 0; j < dim; j++)
        {
            float_obj = PyList_GetItem(points_list_obj, i*dim+j);
            datapoints_arr_mem[(i*dim) + j] = PyFloat_AsDouble(float_obj);
#ifdef DEBUG
            printf("dp: %Lf\n", datapoints_arr_mem[(i*dim) + j]);
#endif
        }
    }

    // TODO: Perform steps 1-7 ...

    if (kmeans(dim, k, MAX_ITER, number_of_datapoints, datapoints_arr_ptr, centroids_arr_ptr) != 0)
        return NULL;

    PyObject* centroids = PyList_New(k * dim);

    if (!centroids)
        return NULL;

    for (int i = 0; i < K * dim; i++) {
        PyObject *num = PyFloat_FromDouble(centroids_arr_mem[i]);

        if (!num) {
            Py_DECREF(centroids);
            return NULL;
        }

        PyList_SET_ITEM(centroids, i, num);
    }

#ifdef DEBUG
    printf("Finished! :D\n");
#endif

    free(centroids_arr_ptr);
    free(centroids_arr_mem);
    free(datapoints_arr_ptr);
    free(datapoints_arr_mem);

    return centroids;
}

static PyObject* perform_subtask(PyObject* self, PyObject* args)
{
    int number_of_datapoints;

    int dim, k, goal;

    PyObject* points_list_obj = NULL;
    PyObject* float_obj = NULL;

    F_TYPE** datapoints_arr_ptr = NULL;
    F_TYPE* datapoints_arr_mem = NULL;

#ifdef DEBUG
    printf("Got to here!\n");
#endif

    if (!PyArg_ParseTuple(args, "iiiO", &dim, &k, &goal &points_list_obj))
        return NULL;

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
    datapoints_arr_ptr = calloc(sizeof(F_TYPE*), number_of_datapoints);
    datapoints_arr_mem = calloc(sizeof(F_TYPE),number_of_datapoints * dim);
    assert(datapoints_arr_ptr != NULL);
    assert(datapoints_arr_mem != NULL);

#ifdef DEBUG
    printf("Start parsing dps!\n");
#endif

    // Parsing datapoints...
    for (int i = 0; i <number_of_datapoints; i++)
    {
        datapoints_arr_ptr[i] = datapoints_arr_mem + (i * dim);

        for (int j = 0; j < dim; j++)
        {
            float_obj = PyList_GetItem(points_list_obj, i*dim+j);
            datapoints_arr_mem[(i*dim) + j] = PyFloat_AsDouble(float_obj);
#ifdef DEBUG
            printf("dp: %Lf\n", datapoints_arr_mem[(i*dim) + j]);
#endif
        }
    }

    PyObject* return_value; // TODO: If the values getting trashed after leaving the cases' scope, allocate according to goal above

    // Performing the requested task
    switch (goal)
    {
        case(2):
        {
            break;
        }
        case(3):
        {
            break;
        }
        case(4):
        {
            break;
        }
        case(5):
        {
            break;
        }
        default:
        {

            break;
        }
    }

    free(datapoints_arr_ptr);
    free(datapoints_arr_mem);

    return return_value;
}

static PyMethodDef methods[] = {
        {"perform_subtask", perform_subtask, METH_VARARGS,
        "Executes a sub-spkmeans function based on a received 'goal' parameter"},
        {"spkmeans_fit", spkmeans_fit, METH_VARARGS,
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