#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#define TRUE (1)
#define FALSE (0)
#define F_TYPE long double
#define F_TYPE_FORMAT_SPEC ("%Lf")
#define F_TYPE_OUTPUT_FORMAT_SPEC ("%.4Lf")
#define F_TYPE_MAX (LDBL_MAX)

int kmeans(int dim, int k, int max_iter, int datapoint_num,
                F_TYPE** datapoints_arr_ptr, F_TYPE** centroids_arr_ptr);

// #define DEBUG

static PyObject* fit(PyObject* self, PyObject* args)
{
    int number_of_datapoints;

    int dim;
    int K;
    int max_iter;

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

    if (!PyArg_ParseTuple(args, "iiiOO", &dim, &K, &max_iter, &centroids_list_obj, &points_list_obj))
        return NULL;

#ifdef DEBUG
    printf("Passed parseargs!\n");

    printf("d:%d\n", dim);
    printf("k:%d\n", K);
    printf("mi:%d\n", max_iter);
#endif

    number_of_datapoints = PyList_Size(points_list_obj) / dim;

#ifdef DEBUG
    printf("n:%d\n", number_of_datapoints);
#endif

    // Allocating according to the received length
    centroids_arr_ptr = calloc(sizeof(F_TYPE*), K);
    centroids_arr_mem = calloc(sizeof(F_TYPE), K * dim);
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


    if (kmeans(dim, K, max_iter, number_of_datapoints, datapoints_arr_ptr, centroids_arr_ptr) != 0)
        return NULL;

    PyObject* centroids = PyList_New(K * dim);

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

static PyMethodDef methods[] = {
        {"fit", fit, METH_VARARGS, "Executes the kmeans algo with the given parameters"},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef module = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        "Performs kmeans algo",
        -1,
        methods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void)
{
    return PyModule_Create(&module);
}




/* calcs the distance between to vectors of dimension dim */
/* the vectors are represented as dim long arrays of F_TYPE */ 
F_TYPE calc_vec_dist(F_TYPE* a, F_TYPE* b, int dim)
{
    F_TYPE ret = 0;
    int i = 0;
    for (i = 0; i < dim; ++i)
    {
        ret += (a[i] - b[i]) * (a[i] - b[i]);
    }
    return ret;
}

/* vector sum: dst = a+b */
void vec_sum(F_TYPE* dst, F_TYPE* a, F_TYPE* b, int dim)
{
    int i = 0;
    for (i = 0; i < dim; ++i)
        dst[i] = a[i] + b[i];
}

/* divide a vector by scalar: dst = a/s */
void vec_div_by_scalar(F_TYPE* dst, F_TYPE* a, F_TYPE s, int dim)
{
    int i = 0;
    for (i = 0; i < dim; ++i)
        dst[i] = a[i] / s;
}

/*  Returns TRUE if a's values are equal to b's values
    of course a and be must be of the same size, which is dim*vec_num . */
int cmp_matrices(F_TYPE** a, F_TYPE** b, int dim, int vec_num)
{
    int i = 0;
    int j = 0;
    for (i = 0; i < vec_num; ++i)
    {
        for (j = 0; j < dim; ++j)
        {
            if (a[i][j] != b[i][j]) return FALSE;
        }
    }
    return TRUE;
}


int assign_to_cluster(F_TYPE* dp, F_TYPE** centroids,
                      int dim, int k)
{
    F_TYPE min_dist = F_TYPE_MAX;
    int min_dist_centrd_idx = 0;
    F_TYPE curr_dist;

    int i;

    for (i = 0; i < k; ++i)
    {
        curr_dist = calc_vec_dist(dp, centroids[i], dim);
        if (curr_dist < min_dist) {
            min_dist = curr_dist;
            min_dist_centrd_idx = i;
        }
    }
    
    return min_dist_centrd_idx;
}

int kmeans(int dim, int k, int max_iter, int datapoint_num,
                F_TYPE** datapoints_arr_ptr, F_TYPE** centroids_arr_ptr)
{
    

    int status = 0;

    int iter_num = 0;

    /* variables to hold the last iteration centroids */ 
    F_TYPE* last_iter_centrds_mem;
    F_TYPE** last_iter_centrds;

    /* temp variables */
    int i = 0;
    int curr_assigned_clstr;


    F_TYPE* centrds_sum_mem = NULL;
    F_TYPE** centrds_sum = NULL;
    int* centrds_ref_cnt = NULL;


    /*===================*/
    /* memory allocation */
    /*===================*/
    last_iter_centrds_mem = calloc(sizeof(F_TYPE), k*dim);
    assert(last_iter_centrds_mem != NULL);
    last_iter_centrds = calloc( sizeof(F_TYPE*) , k);
    assert(last_iter_centrds != NULL);
    for (i = 0; i < k; ++i)
    {
        last_iter_centrds[i] = last_iter_centrds_mem + (i*dim);
    }


    centrds_sum_mem = calloc(sizeof(F_TYPE), k*dim);
    assert(centrds_sum_mem != NULL);
    centrds_sum = calloc(sizeof(F_TYPE*), k);
    assert(centrds_sum != NULL);
    for (i = 0; i < k; ++i)
    {
        centrds_sum[i] = centrds_sum_mem + (i*dim);
    }
    centrds_ref_cnt = calloc(sizeof(int), k);
    assert(centrds_ref_cnt != NULL);


    // NO NEED TO INTIALIZE CENTROIDS - THIS IS DONE BY KMEANS++


    /*======================*/
    /* algorithm iterations */
    /*======================*/
    while (TRUE) {


        /* check max_iter break condition */
        /*--------------------------------*/
        /* the other condition to finish the run appears later */
        if ((max_iter > 0) && (iter_num >= max_iter)) break;


        /* initialize centroids ref count and sum */
        /*----------------------------------------*/        
        memset(centrds_ref_cnt, 0, sizeof(int)*k);
        for (i = 0; i < k; ++i)
        {
            memset(centrds_sum[i], 0, sizeof(F_TYPE)*dim);
        }       


        /* assign each datapoint to the closest centroid */
        /*-----------------------------------------------*/     
        for (i = 0; i < datapoint_num; ++i)
        {
            curr_assigned_clstr = assign_to_cluster(datapoints_arr_ptr[i], centroids_arr_ptr, dim, k);
            // output_cluster_assign[i] = curr_assigned_clstr;

            centrds_ref_cnt[curr_assigned_clstr]++;

            vec_sum(centrds_sum[curr_assigned_clstr], 
                centrds_sum[curr_assigned_clstr], datapoints_arr_ptr[i], dim);
        }


        /* update centroids */
        /*------------------*/
        for (i = 0; i < k; ++i)
        {
            vec_div_by_scalar(centroids_arr_ptr[i], centrds_sum[i], centrds_ref_cnt[i], dim);
        }


        /* check unchanging centroids end condition */
        /*------------------------------------------*/
        if ((iter_num != 0) && /* to avoid that the default values are received 
                                  in the first iteration. This extra condition 
                                  can cause 1 redundant iteration at worst... */
                (cmp_matrices(centroids_arr_ptr, last_iter_centrds, dim, k))) {
            break;
        }


        /* update last iteration centroids */
        /*---------------------------------*/
        for (i = 0; i < k; ++i)
        {
            memcpy(last_iter_centrds[i], centroids_arr_ptr[i], sizeof(F_TYPE)*dim);
        }
        iter_num++;
    }


    free(last_iter_centrds);
    free(last_iter_centrds_mem);
    free(centrds_sum_mem);
    free(centrds_sum);
    free(centrds_ref_cnt);

    return status;

}