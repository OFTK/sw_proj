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



/*=======*/
/* debug */
/*=======*/
#define DEBUG

#ifdef DEBUG
#	define DEBUG_PRINT(x) printf(x) 
#else
#	define DEBUG_PRINT(x) do {} while (0)
#endif


/*======*/
/* misc */
/*======*/
/* TODO: for an error in the input file, should i print error or invalid_input? */
/* TODO: replace all asserts for error prints */
#define ERROR_PRINT() printf("An Error Has Occured")
#define INVALID_INPUT_PRINT() printf("Invalid Input!")

/* goal enum:
	0	: spk
	1	: wam
	2	: ddg
	3	: lnorm
	4 	: jacobi
	(-1): BAD ARGUMENT */
int goal_enum(const char* goal_str) {
	if 		(0 == strcmp("spk", goal_str)) return 0;
	else if (0 == strcmp("wam", goal_str)) return 1;
	else if (0 == strcmp("ddg", goal_str)) return 2;
	else if (0 == strcmp("lnorm", goal_str)) return 3;
	else if (0 == strcmp("jacobi", goal_str)) return 4;
	return -1;
}



/*============*/
/* math utils */
/*============*/
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

/* 	Returns TRUE if a's values are equal to b's values
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


int kmeans(int dim,
				F_TYPE** input_dps, int dp_num,
				F_TYPE** output_centrds, int* output_cluster_assign, int k,
				int max_iter)
{
	int status = 0;

	int iter_num = 0;

	/* variables to hold the last iteration centroids */ 
	F_TYPE* last_iter_centrds_mem;
	F_TYPE** last_iter_centrds;

	/* temp variables */
	int i = 0;
	/* int j = 0; */
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



	/*==========================*/
	/* centroids initialization */
	/*==========================*/
	
	for (i = 0; i < k; ++i)
	{
		/* copy the i'th vector to the i'th centroid */
		memcpy(output_centrds[i], input_dps[i], sizeof(F_TYPE)*dim);
	}



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
		for (i = 0; i < dp_num; ++i)
		{
			curr_assigned_clstr = assign_to_cluster(input_dps[i], output_centrds, dim, k);
			output_cluster_assign[i] = curr_assigned_clstr;

			centrds_ref_cnt[curr_assigned_clstr]++;

			vec_sum(centrds_sum[curr_assigned_clstr], 
				centrds_sum[curr_assigned_clstr], input_dps[i], dim);
		}


		/* update centroids */
		/*------------------*/
		for (i = 0; i < k; ++i)
		{
			vec_div_by_scalar(output_centrds[i], centrds_sum[i], centrds_ref_cnt[i], dim);
		}


		/* check unchanging centroids end condition */
		/*------------------------------------------*/
		if ((iter_num != 0) && /* to avoid that the default values are received 
								  in the first iteration. This extra condition 
								  can cause 1 redundant iteration at worst... */
				(cmp_matrices(output_centrds, last_iter_centrds, dim, k))) {
			break;
		}


		/* update last iteration centroids */
		/*---------------------------------*/
		for (i = 0; i < k; ++i)
		{
			memcpy(last_iter_centrds[i], output_centrds[i], sizeof(F_TYPE)*dim);
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


void print_centroids(F_TYPE** centroids, int dim, int k)
{
	int i = 0;
	int j = 0;

	for (i = 0; i < k; ++i)
	{
		for (j = 0; j < dim; ++j)
		{
			printf(F_TYPE_OUTPUT_FORMAT_SPEC, centroids[i][j]);
			if (j != dim-1) printf(",");
		}
		printf("\n");
	}
}

/* 	write F_TYPE to x*.
	return 0 when scanning F_TYPE then ","
	return 1 when scanning F_TYPE then "\n"
	return 2 when scanning eof
	return -1 otherwise */
int scan_next_val(F_TYPE* x, FILE* finput) {

	int value_exists = fscanf(finput, F_TYPE_FORMAT_SPEC, x);
	char c = 0;
	int next_char = fscanf(finput, "%c", &c);

	if ((1 == value_exists) && (1 == next_char) && (44 == c)) /* num with a comma */
		return 0;

	if ((1 == value_exists) && (1 == next_char) && (10 == c)) /* num with LF */
		return 1;

	if (EOF == value_exists) /* EOF */
		return 2;

	return -1;
}

int main(int argc, char const *argv[])
{

	/* return value of main is this status: 
	-1	means error
	 0 	means ok
	*/
	int status = 0;


	/* input arguments */
	int k = 0;
	int max_iter = -1;  /* -1 symbolized no upper bound for max_iter. */
 	int goal = -1;
 	FILE* finput;


	/* for allocation of datapoints */
	F_TYPE** input_dps = NULL;
	F_TYPE* input_dp_mem = NULL;
	int dim = 0;
	int dp_num = 0;


	/* memory for returned centroids */
	F_TYPE** output_centroids = NULL;
	F_TYPE* output_centroids_mem = NULL;

	/* memory for output cluster assignments */
	int* output_cluster_assign = NULL;

	/* temporary variables */
	int i = 0;
	F_TYPE scanned_num = 0;
	int curr_dim = 0;
	F_TYPE* curr_vector = NULL;
	int scan_status = 0;

	/*===================*/
	/* parsing arguments */
	/*===================*/

	if (3 != argc) {
		DEBUG_PRINT("main: bad num of args, expecting 3\n");
		INVALID_INPUT_PRINT();
		status = -1;
		goto finalize;
	}

	k = atoi(argv[1]);
	goal = goal_enum(argv[2]);
	if (-1 == goal) {
		DEBUG_PRINT("bad goal value: expected one of {spk, wam, ddg, lnorm, jacobi}");
		INVALID_INPUT_PRINT();
		goto finalize;
	}

	/* TODO: should error here be invalid input or error? */
	finput = fopen(argv[3], "r");
	if (NULL == finput) {
		DEBUG_PRINT("fopen returned with error");
		INVALID_INPUT_PRINT();
		goto finalize;
	}

	/*=================*/
	/* scanning inputs */
	/*=================*/

	/* scan stdin and build the datapoints in a matrix */
	while (TRUE) {
		scan_status = scan_next_val(&scanned_num, finput);

		if (scan_status < 0) {
			DEBUG_PRINT("bad input file format\n");
			INVALID_INPUT_PRINT();
			status = -1;
			goto close;
		}

		/* when done reading the file */
		if (scan_status == 2) {
			if (dim == 0) {
				DEBUG_PRINT("not able to read any num from file.\n");
				INVALID_INPUT_PRINT();
				status = -1;
				goto close;

			}
			if (curr_dim != 0) {
				#ifdef DEBUG
				printf("file terminated in the middle of a vector (dp num %d, in idx %d)\n",
					dp_num, curr_dim);
				#endif
				INVALID_INPUT_PRINT();
				status = -1;
				goto close;
			}

			/* if scanned succesfully - allocate pointers array */
			input_dps = realloc(input_dps, (sizeof(F_TYPE*)*dp_num));
			assert(input_dps != NULL);
			for (i = 0; i < dp_num; ++i)
			{
				input_dps[i] = input_dp_mem+(i*dim);
			}
			break;
		}


		/* when allocating the first vector: reallocate entire vector */
		if (dp_num == 0) {
			dim++;
			curr_vector = realloc(curr_vector, (sizeof(F_TYPE)*dim));
			if NULL
			assert(curr_vector != NULL);
		}

		/* add the scanned F_TYPE to the current vector */
		curr_dim++;
		if (curr_dim > dim) {
			#ifdef DEBUG
			printf("bad input vector length: first vector is of length %d while %d'th vector is of length %d\n", 
				dim, (dp_num+1), curr_dim);
			#endif
			INVALID_INPUT_PRINT();
			status = -1;
			goto close;
		}
		curr_vector[curr_dim-1] = scanned_num;


		/* 	when scanned full vector length: reallocate the
			matrix and memcpy the vector to the matrix */
		if (scan_status == 1) {
			if (curr_dim != dim) {
				printf("bad input vector length: first vector is of length %d while %d'th vector is of length %d\n", 
					dim, (dp_num+1), curr_dim);
				status = -1;
				goto close;
			}

			dp_num++;

			input_dp_mem = realloc(input_dp_mem, (sizeof(F_TYPE)*(dp_num*dim)));
			assert(input_dp_mem != NULL);
			memcpy(input_dp_mem+((dp_num-1)*dim), curr_vector, sizeof(F_TYPE)*dim);

			curr_dim = 0;
		}
	}




	/*=====================*/
	/* algorithm procedure */
	/*=====================*/
	/* allocate memory for output centroids */
	output_centroids_mem = calloc(sizeof(F_TYPE), k*dim);
	assert(output_centroids_mem != NULL);
	output_centroids = calloc(sizeof(F_TYPE*), k);
	assert(output_centroids != NULL);
	for (i = 0; i < k; ++i)
	{
		output_centroids[i] = output_centroids_mem + (i*dim);
	}

	/* allocate memory for datapoint assignment to cluster */
	output_cluster_assign = calloc(sizeof(int), dp_num);
	assert(output_cluster_assign != NULL);

	
	/* THE KMEANS PROCEDURE CALL */
	status = kmeans(dim, input_dps,
		dp_num, output_centroids, output_cluster_assign,
		k, max_iter);


	/* print or return the output centroids */
	print_centroids(output_centroids, dim, k);


	/* freeing centroids mem */
	free(output_centroids_mem);
	free(output_centroids);
	free(output_cluster_assign);

	close:
		fclose(finput);

	finalize:	
		free(curr_vector);
		free(input_dp_mem);
		free(input_dps);

	return status;
}