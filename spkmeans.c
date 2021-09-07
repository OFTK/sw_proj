#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <float.h>

#define TRUE (1)
#define FALSE (0)
#define F_TYPE double
#define F_TYPE_FORMAT_SPEC ("%lf")
#define F_TYPE_OUTPUT_FORMAT_SPEC ("%.4f")
#define F_TYPE_MAX (DBL_MAX)
#define F_TYPE_ABS(x) (fabs(x))
#define JACOBI_CONVERGENCE_SIGMA (0.001)

typedef struct f_and_idx {
	F_TYPE f;
	int idx;
} f_and_idx;

enum status {
	Error,
	Success,
	Finish
};


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
#define PRINT_ERROR() (printf("An Error Has Occured"))
#define PRINT_INVALID_INPUT() (printf("Invalid Input!"))

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

/* inverse sqare-root */
/* we use a define to save the redundant function-call run-time */
#define INV_SQRT(x) (1 / (sqrt(x)))

/* calcs the distance between to vectors of dimension dim */
/* the vectors are represented as dim long arrays of F_TYPE */ 
F_TYPE calc_l2_norm(F_TYPE* a, F_TYPE* b, int dim)
{
	F_TYPE ret = 0;
	int i = 0;
	for (i = 0; i < dim; ++i)
		ret += (a[i] - b[i]) * (a[i] - b[i]);

	return sqrt(ret);
}

/* vector sum: dst = a+b */
void vec_sum(F_TYPE* dst, F_TYPE* a, F_TYPE* b, int dim)
{
	int i = 0;
	for (i = 0; i < dim; ++i)
		dst[i] = a[i] + b[i];
}

/* sum of elements in a vector */
F_TYPE element_sum(F_TYPE* vec, int dim)
{
	int i = 0;
	F_TYPE sum = 0;
	for (i = 0; i < dim; ++i)
		sum = sum + vec[i];

	return sum;
}

/* divide a vector by scalar: dst = a/s */
void vec_div_by_scalar(F_TYPE* dst, F_TYPE* a, F_TYPE s, int dim)
{
	int i = 0;
	for (i = 0; i < dim; ++i)
		dst[i] = a[i] / s;
}

/* 	Returns TRUE if a's values are equal to b's values
	of course a and b must be of the same size, which is dim*vec_num. */
int cmp_matrices(F_TYPE** a, F_TYPE** b, int dim, int vec_num)
{
	int i = 0;
	int j = 0;
	for (i = 0; i < vec_num; ++i)
	{
		for (j = 0; j < dim; ++j)
			if (a[i][j] != b[i][j]) return FALSE;
	}
	return TRUE;
}


/* 	Transoposes a matrix. T( mtx(nXm) ) => o_t_mtx(mXn) */
void transpose_matrix(F_TYPE** mtx, F_TYPE**o_t_mtx, int n, int m)
{
	int i = 0;
	int j = 0;
	for (i = 0; i < n; ++i) {
		for (j = 0; j < m; ++j)
			o_t_mtx[j][i] = mtx[i][j];
	}
}


/*	Calculates the multiplication of square metrices
	a and b. Then insert it's values to mul. */
/* 	TODO: if matrices are diagonal or almost diagonal 
	you can calculate in linear time */
void mul_square_matrices(
	F_TYPE** a, F_TYPE** b, 
	int n, F_TYPE** mul)
{	
	int i = 0, j = 0, k = 0;
	
	for(i=0; i < n; i++) {    
		for(j=0; j < n; j++) {

			mul[i][j]=0;
			for(; k < n; k++)    
				mul[i][j] += a[i][k]*b[k][j];    

		}   
	}
}



/*===========*/
/* algorithm */
/*===========*/

enum status calc_weighted_adjacency_matrix(
	int dim,
	F_TYPE** input_dps, int dp_num,
	F_TYPE** wam)
{
	enum status s = Success;
	int i, j;

	/* calc weight */
	/*-------------*/
	for (i = 0; i < dp_num; ++i) {
		for (j = 0; j < i; ++j) {
			wam[i][j] = exp((-0.5) * calc_l2_norm(input_dps[i], input_dps[j], dim));
			if (0 != errno) s = Error;
		}
	}

	/* make matrix symmetric */
	/*-----------------------*/
	for (i = dp_num-1; i >= 0; --i) {
		for (j = dp_num-1; j>=i; --j) {
			wam[i][j] = wam[j][i];
		}
	}

	return s;
}

enum status calc_diagonal_degree_matrix(
	F_TYPE** wam, int dp_num,
	F_TYPE** ddg)
{
	int i;
	for (i = 0; i < dp_num; ++i)
		ddg[i][i] = element_sum(wam[i], dp_num);

	return Success;
}


enum status calc_normalized_graph_laplacian(
	F_TYPE** wam, F_TYPE** ddg, int dp_num,
	F_TYPE** lnorm)
{
	int i, j;
	enum status status = Success;
	F_TYPE eye = 0;

	/* calc inverse-sqare-root of ddg */
	/*--------------------------------*/

	/* memory allocation */
	F_TYPE** ddg_invsqrt = NULL;
	F_TYPE* ddg_invsqrt_mem = NULL;
	ddg_invsqrt_mem = calloc(sizeof(F_TYPE), dp_num*dp_num);
	assert(ddg_invsqrt_mem != NULL);
	ddg_invsqrt = calloc(sizeof(F_TYPE*), dp_num);
	assert(ddg_invsqrt != NULL);
	for (i = 0; i < dp_num; ++i)
		ddg_invsqrt[i] = ddg_invsqrt_mem + (i*dp_num);

	/* inverse square root of matrix main diagonal */
	for (i = 0; i < dp_num; ++i)
		ddg_invsqrt[i][i] = INV_SQRT(ddg[i][i]);

	if (0 != errno) status = Error;


	/* calc lnorm */
	/*------------*/

	for (i = 0; i < dp_num; ++i) {
		for (j = 0; j < i; ++j) {

			/* calc eye matrix */
			if (i == j) eye = 1; else eye = 0;
			/* calc lnorm */
			lnorm[i][j] = eye - (ddg_invsqrt[i][i]*wam[i][j]*ddg_invsqrt[j][j]);
		}
	}

	free(ddg_invsqrt);
	free(ddg_invsqrt_mem);

	return status;
}



int assign_to_cluster(
	F_TYPE* dp, F_TYPE** centroids,
	int dim, int k)
{
	F_TYPE min_dist = F_TYPE_MAX;
	int min_dist_centrd_idx = 0;
	F_TYPE curr_dist;

	int i;

	for (i = 0; i < k; ++i)
	{
		curr_dist = pow(calc_l2_norm(dp, centroids[i], dim), 2);
		if (curr_dist < min_dist) {
			min_dist = curr_dist;
			min_dist_centrd_idx = i;
		}
	}
	
	return min_dist_centrd_idx;
}



/**
 * @details Calculates the (Frobenius Norm(mtx))^2 - sum((diagonal values)^2), 
 * 			which is called 'off(mtx)^2' of mtx in the specification document.
 */
int calc_off(F_TYPE** mtx, size_t mtx_columns_num)
{
	int off = 0;
	size_t i = 1, j = 0;

	for (;i < mtx_columns_num; i++)
	{
		for (;j < i; j++)
			off += 2 * pow(mtx[i][j], 2);
	}

	return off;
}

/**
 * @details Calculates the next A' and P, given the current A.
 */ 
enum status calc_jacobi_iteration(
	F_TYPE** io_mtx_A, int dp_num,
	F_TYPE** o_mtx_P)
{	
	/* TODO: Why were these declared size_t? */
	/*size_t i = 0, j = 0, max = 0, k = 0, l = 0;*/
	int k = 0;
	int l = 0;
	int i = 0;
	int j = 0;
	F_TYPE max = 0;

	F_TYPE theta = 0, t = 0, c = 0, s = 0, temp = 0;

	/* TODO : Find the largest off-diagonal, absolute value A_ij in io_mtx_A */
	for (k = 1; k < dp_num; k++)
	{
		for (; l < k; l++)
		{
			if (F_TYPE_ABS(io_mtx_A[k][l]) > max)
			{
				i = k;
				j = l;
				max = F_TYPE_ABS(io_mtx_A[k][l]);
			}
		}
	}
	
	if (max != 0)
	{
		/* Calculate Theta, t, c */
		theta = (io_mtx_A[j][j] - io_mtx_A[i][i]) / (2*io_mtx_A[i][j]);

		if (theta > 0)
			t = 1 / theta + sqrt(pow(theta, 2) + 1);

		else
			t = -1 / -theta + sqrt(pow(theta, 2) + 1);

		c = 1 / sqrt(pow(t, 2) + 1);
		s = t*c;

		/* Calculate the matrix rotation matrix mtx_P */
		memset(o_mtx_P, 0, dp_num*dp_num*sizeof(F_TYPE));

		for (k = 0; k < dp_num; k++)
			o_mtx_P[k][k] = 1;

		o_mtx_P[i][i] = c;
		o_mtx_P[j][j] = c;
		o_mtx_P[i][j] = s;
		o_mtx_P[j][i] = -s;

		/* Performing step as described in sub-paragraph 6 */
		for (k = 0; k < dp_num; k++)
		{
			if (k != i && k != j)
			{
				io_mtx_A[k][i] = c*io_mtx_A[k][i] - s*io_mtx_A[k][j];
				io_mtx_A[i][k] = io_mtx_A[k][i];

				io_mtx_A[k][j] = c*io_mtx_A[k][j] + s*io_mtx_A[k][i];
				io_mtx_A[j][k] = io_mtx_A[k][j];
			}
		}

		temp = io_mtx_A[i][i] - io_mtx_A[j][j];

		io_mtx_A[i][i] = pow(c,2)*io_mtx_A[i][i] + pow(s,2)*io_mtx_A[j][j]
						- 2*s*c*io_mtx_A[i][j];

		io_mtx_A[j][j] = pow(s,2)*io_mtx_A[i][i] + pow(c,2)*io_mtx_A[j][j]
						+ 2*s*c*io_mtx_A[i][j];

		io_mtx_A[i][j] = (pow(c,2) - pow(s,2))*io_mtx_A[i][j] + s*c*temp;

		return Success;
	}

	return Finish;
}

/**
 * @details This method received a symmetric matrix, A. It returns a status value,
 * 			 an eignvalues matrix pointed by mtx_A (changes the received matrix) and
 * 			 an eignvectors matrix pointed by mtx_V.
 * @note
 * 		- This method perform in-place work on mtx_A. If you wish to save the values
 * 		  of that matrix, you should save it before using this function.
 * 		- This method does not check the validity of mtx_A, this method should be
 * 		  used cautiously only on symmetric matrices.
 */ 
enum status find_eigenvalues_jacobi(
	F_TYPE** io_mtx_A, int dp_num,
	F_TYPE** o_mtx_V)
{	
	/* Declare locals */
	F_TYPE** mtx_P = NULL, **mtx_V_temp = NULL;
	enum status status = Success;
	int prev_off = 0, curr_off = 0;

	/* Allocate locals */
	mtx_P = (F_TYPE**)calloc(dp_num*dp_num, sizeof(F_TYPE));
	assert(mtx_P != NULL);
	mtx_V_temp = (F_TYPE**)calloc(dp_num*dp_num, sizeof(F_TYPE));
	assert(mtx_V_temp != NULL);

	/* Doing the first loop's step in order to save the first P in o_mtx_V */
	prev_off = calc_off(io_mtx_A, dp_num);

	/* Calculating the first A' and P matrices */
	status = calc_jacobi_iteration(io_mtx_A, dp_num, o_mtx_V);

	curr_off = calc_off(io_mtx_A, dp_num);

	/* Calc the diagonal A' matrix */
	while (JACOBI_CONVERGENCE_SIGMA > curr_off - prev_off && status != Finish)
	{
		/* Calculating the next A and P matrices */
		status = calc_jacobi_iteration(io_mtx_A, dp_num, mtx_P);
		
		/* Calculating the next V matrix */
		memcpy(mtx_V_temp, o_mtx_V, dp_num*dp_num*sizeof(F_TYPE));
		mul_square_matrices(mtx_V_temp, mtx_P, dp_num, o_mtx_V);

		/* Calculating the function 'off^2' of the new matrix */
		prev_off = curr_off;
		curr_off = calc_off(io_mtx_A, dp_num);
	}

	/* free locals */
	free(mtx_P);

	return Success;
}


/* comparison function to sort eigenvalues for eigengap heuristic */
int f_and_idx_compare (const void* a, const void* b) {
	f_and_idx* f_a = (f_and_idx*)a;
	f_and_idx* f_b = (f_and_idx*)b;

	/* a > b  => return > 0*/
	if (f_a->f > f_b->f) return 1; 
	/* a < b  => return < 0*/
	if (f_a->f < f_b->f) return (-1); 

	/* a = b  => return according to idx*/
	if (f_a->idx > f_b->idx) return 1; 
	if (f_a->idx < f_b->idx) return (-1);

	return 0;
}


/* the eigengap heuristic */
enum status eigengap_heuristic(
	f_and_idx* eigenvalues, int n,
	int* o_k) /* output - estimated k */
{
	F_TYPE* eigengaps = NULL;
	int status = Success;
	int i = 0;
	F_TYPE max_gap = 0;

	/* calculate eigengaps */
	eigengaps = calloc(sizeof(F_TYPE), (n-1));
	if (NULL == eigengaps) {
		status = Error;
		return status;
	}

	for (i = 0; i < (n-1); ++i)
		eigengaps[i] = F_TYPE_ABS(eigenvalues[i+1].f - eigenvalues[i].f);

	/* find argmax */
	o_k = 0;
	for (i = 0; i < (n/2); ++i) /* TODO: find out if n/2 or n instead */
	{
		if (eigengaps[i] > max_gap)
		{
			max_gap = eigengaps[i];
			*o_k = i;
		}
	}

	return status;
}


/* kmeans algorithm, from hw 1 */
enum status kmeans(
	F_TYPE** input_dps, int dp_num, int dim,
	F_TYPE** output_centrds, int* output_cluster_assign,
	int k, int max_iter)
{
	int status = Success;

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
		last_iter_centrds[i] = last_iter_centrds_mem + (i*dim);


	centrds_sum_mem = calloc(sizeof(F_TYPE), k*dim);
	assert(centrds_sum_mem != NULL);
	centrds_sum = calloc(sizeof(F_TYPE*), k);
	assert(centrds_sum != NULL);
	for (i = 0; i < k; ++i)
		centrds_sum[i] = centrds_sum_mem + (i*dim);
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
			memset(centrds_sum[i], 0, sizeof(F_TYPE)*dim);


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
			vec_div_by_scalar(output_centrds[i], centrds_sum[i], centrds_ref_cnt[i], dim);


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
			memcpy(last_iter_centrds[i], output_centrds[i], sizeof(F_TYPE)*dim);
		iter_num++;
	}


	free(last_iter_centrds);
	free(last_iter_centrds_mem);
	free(centrds_sum_mem);
	free(centrds_sum);
	free(centrds_ref_cnt);

	return status;
}


void print_matrix(F_TYPE** matrix, int n, int m)
{
	int i = 0;
	int j = 0;

	for (i = 0; i < n; ++i)
	{
		for (j = 0; j < m; ++j)
		{
			printf(F_TYPE_OUTPUT_FORMAT_SPEC, matrix[i][j]);
			if (j != m-1) printf(",");
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

	int status = Success;


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


	/* memory for weighted adjacency matrix */
	F_TYPE** wam = NULL;
	F_TYPE* wam_mem = NULL;

	/* memory for diagonal degree matrix */
	F_TYPE** ddg = NULL;
	F_TYPE* ddg_mem = NULL;

	/* memory for normalized graph Laplacian */
	F_TYPE** lnorm = NULL;
	F_TYPE* lnorm_mem = NULL;

	/* memory for normalized graph Laplacian */
	F_TYPE** eigenvectors = NULL;
	F_TYPE* eigenvectors_mem = NULL;
	f_and_idx* eigenvalues = NULL;

	F_TYPE** k_eigenvectors = NULL;
	F_TYPE* k_eigenvectors_mem = NULL;

	F_TYPE* eigenvalues_only = NULL;


	/* memory for returned centroids */
	F_TYPE** u_mtx = NULL;
	F_TYPE* u_mtx_mem = NULL;
	F_TYPE** t_mtx = NULL;
	F_TYPE* t_mtx_mem = NULL;
	F_TYPE* dummy_vector = NULL;

	/* memory for returned centroids */
	F_TYPE** output_centroids = NULL;
	F_TYPE* output_centroids_mem = NULL;

	/* memory for output cluster assignments */
	int d = 0;
	int n = 0;
	int* output_cluster_assign = NULL;

	/* temporary variables */
	int i = 0;
	int j = 0;
	F_TYPE scanned_num = 0;
	int curr_dim = 0;
	F_TYPE* curr_vector = NULL;
	int scan_status = 0;

	/*===================*/
	/* parsing arguments */
	/*===================*/

	if (3 != argc) {
		DEBUG_PRINT("main: bad num of args, expecting 3\n");
		PRINT_INVALID_INPUT();
		return Error;
	}

	k = atoi(argv[1]);
	goal = goal_enum(argv[2]);
	if (-1 == goal) {
		DEBUG_PRINT("bad goal value: expected one of {spk, wam, ddg, lnorm, jacobi}");
		PRINT_INVALID_INPUT();
		return Error;
	}

	/* TODO: should error here be invalid input or error? */
	finput = fopen(argv[3], "r");
	if (NULL == finput) {
		DEBUG_PRINT("fopen returned with error");
		PRINT_INVALID_INPUT();
		return Error;
	}

	/*=================*/
	/* scanning inputs */
	/*=================*/

	/* scan stdin and build the datapoints in a matrix */
	while (TRUE) {
		scan_status = scan_next_val(&scanned_num, finput);

		if (scan_status < 0) {
			DEBUG_PRINT("bad input file format\n");
			PRINT_INVALID_INPUT();
			status = Error;
			goto close;
		}

		/* when done reading the file */
		if (scan_status == 2) {
			if (dim == 0) {
				DEBUG_PRINT("not able to read any num from file.\n");
				PRINT_INVALID_INPUT();
				status = Error;
				goto close;

			}
			if (curr_dim != 0) {
				#ifdef DEBUG
				printf("file terminated in the middle of a vector (dp num %d, in idx %d)\n",
					dp_num, curr_dim);
				#endif
				PRINT_INVALID_INPUT();
				status = Error;
				goto close;
			}

			/* if scanned succesfully - allocate pointers array */
			input_dps = realloc(input_dps, (sizeof(F_TYPE*)*dp_num));
			assert(input_dps != NULL);
			for (i = 0; i < dp_num; ++i)
				input_dps[i] = input_dp_mem+(i*dim);

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
			PRINT_INVALID_INPUT();
			status = Error;
			goto close;
		}
		curr_vector[curr_dim-1] = scanned_num;


		/* 	when scanned full vector length: reallocate the
			matrix and memcpy the vector to the matrix */
		if (scan_status == 1) {
			if (curr_dim != dim) {
				printf("bad input vector length: first vector is of length %d while %d'th vector is of length %d\n", 
					dim, (dp_num+1), curr_dim);
				status = Error;
				goto close;
			}

			dp_num++;

			input_dp_mem = realloc(input_dp_mem, (sizeof(F_TYPE)*(dp_num*dim)));
			assert(input_dp_mem != NULL);
			memcpy(input_dp_mem+((dp_num-1)*dim), curr_vector, sizeof(F_TYPE)*dim);

			curr_dim = 0;
		}
	}
	free(curr_vector);
	fclose(finput);
	goto algorithm;

	close:
		free(input_dp_mem);
		free(input_dps);
		free(curr_vector);
		fclose(finput);
		return status;


	algorithm:

	/*===========================*/
	/* Weighted Adjacency Matrix */
	/*===========================*/

	/* allocate memory for WAM */
	wam_mem = calloc(sizeof(F_TYPE), dp_num*dp_num);
	assert(wam_mem != NULL);
	wam = calloc(sizeof(F_TYPE*), dp_num);
	assert(wam != NULL);
	for (i = 0; i < dp_num; ++i)
		wam[i] = wam_mem + (i*dp_num);

	/* calc Weighted Adjacency Matrix */
	status = calc_weighted_adjacency_matrix(dim, input_dps, dp_num, wam);

	/* free no-longed needed memory */
	free(input_dps); free(input_dp_mem);

	/* finish run if needed */
	if ((Success != status) || (1 == goal)) {
		free(wam); free(wam_mem);

		if (Success != status)
			PRINT_ERROR();
		else
			print_matrix(wam, dp_num, dp_num);

		return status;
	}


	/*========================*/
	/* Diagonal Degree Matrix */
	/*========================*/

	/* allocate memory for DDG */
	ddg_mem = calloc(sizeof(F_TYPE), dp_num*dp_num);
	assert(ddg_mem != NULL);
	ddg = calloc(sizeof(F_TYPE*), dp_num);
	assert(ddg != NULL);
	for (i = 0; i < dp_num; ++i)
		ddg[i] = ddg_mem + (i*dp_num);

	/* calc Diagonal Degree Matrix */
	status = calc_diagonal_degree_matrix(wam, dp_num, ddg);


	/* finish run if needed */
	if ((Success != status)  || (goal == 2)) {
		free(wam); free(wam_mem);
		free(ddg); free(ddg_mem);

		if (Success != status)
			PRINT_ERROR();
		else
			print_matrix(ddg, dp_num, dp_num);

		return status;
	}


	/*============================*/
	/* Normalized Graph Laplacian */
	/*============================*/

	/* allocate memory for lnorm */
	lnorm_mem = calloc(sizeof(F_TYPE), dp_num*dp_num);
	assert(lnorm_mem != NULL);
	lnorm = calloc(sizeof(F_TYPE*), dp_num);
	assert(lnorm != NULL);
	for (i = 0; i < dp_num; ++i)
		lnorm[i] = lnorm_mem + (i*dp_num);

	status = calc_normalized_graph_laplacian(wam, ddg, dp_num, lnorm);

	/* free no-longer needed memory */
	free(wam); free(wam_mem);
	free(ddg); free(ddg_mem);


	if ((Success != status) || (3 == goal)) {
		free(lnorm); free(lnorm_mem);

		if (Success != status)
			PRINT_ERROR();
		else 
			print_matrix(lnorm, dp_num, dp_num);

		return status;
	}


	/*==============================*/
	/* Jacobi Eigenvalues procedure */
	/*==============================*/

	/* allocate memory for eigenvectors matrix */
	eigenvectors_mem = calloc(sizeof(F_TYPE), dp_num*dp_num);
	assert(eigenvectors_mem != NULL);
	eigenvectors = calloc(sizeof(F_TYPE*), dp_num);
	assert(eigenvectors != NULL);
	for (i = 0; i < dp_num; ++i)
		eigenvectors[i] = eigenvectors_mem + (i*dp_num);

	if (Success != status) {
		free(lnorm); free(lnorm_mem);
		goto end_jacobi_stage;
	}

	/**
	* @details This method received a symmetric matrix, A. It returns a status value,
	* 			 an eignvalues matrix pointed by mtx_A (changes the received matrix) and
	* 			 an eignvectors matrix pointed by mtx_V.
	* @note
	* 		- This method perform in-place work on mtx_A. If you wish to save the values
	* 		  of that matrix, you should save it before using this function.
	* 		- This method does not check the validity of mtx_A, this method should be
	* 		  used cautiously only on symmetric matrices.
	*/ 
	status =  find_eigenvalues_jacobi(lnorm, dp_num, eigenvectors);


	/* allocate memory for eigenvalues */
	eigenvalues = calloc(sizeof(f_and_idx), dp_num);

	if (NULL == eigenvalues) {
		free(lnorm); free(lnorm_mem);
		status = Error;
		goto end_jacobi_stage;
	}

	for (i = 0; i < dp_num; ++i){
		eigenvalues[i].f = lnorm[i][i];
		eigenvalues[i].idx = i;
	}


	/* sort eigenvalues, with respect to original indices */
	/*----------------------------------------------------*/
	qsort(eigenvalues, sizeof(f_and_idx), dp_num, &f_and_idx_compare);

	/* run eigengap heuristic - if needed */
	/* TODO: make sure eigenvalues are greater than 0 */
	if (0 == k)
		status = eigengap_heuristic(eigenvalues, dp_num, &k);

	/* save only k eigenvectors */
	/*--------------------------*/
	/* allocate memory for eigenvectors matrix */
	k_eigenvectors_mem = calloc(sizeof(F_TYPE), k*dp_num);
	assert(k_eigenvectors_mem != NULL);
	k_eigenvectors = calloc(sizeof(F_TYPE*), k);
	assert(k_eigenvectors != NULL);
	for (i = 0; i < k; ++i)
		k_eigenvectors[i] = k_eigenvectors_mem + (i*dp_num);


	for (i = 0; i < k; ++i)
		memcpy(
			k_eigenvectors[i], /* dest */
			eigenvectors[eigenvalues[i].idx], /* source - eigenvector from corresponding idx */
			sizeof(F_TYPE)*dp_num /* size */
		);

	free(eigenvectors); free(eigenvectors_mem);

	if ((Success != status) || (4 == goal)) {

		end_jacobi_stage:
		if (Success != status)
			PRINT_ERROR();
		else { /* TODO: is this what i'm supposed to print? */
			eigenvalues_only = calloc(sizeof(F_TYPE*), dp_num);
			if (NULL == eigenvalues) {
				free(eigenvalues); free(k_eigenvectors); free(k_eigenvectors_mem);
				return Error;
			}
			for (i = 0; i < dp_num; ++i)
				eigenvalues_only[i] = eigenvalues[i].f;

			print_matrix(&eigenvalues_only, 1, dp_num);
			print_matrix(k_eigenvectors, k, dp_num);
		}

		free(eigenvalues); free(k_eigenvectors); free(k_eigenvectors_mem);

		return status;
	}

	free(eigenvalues);


	/*=========================================*/
	/* preparing matrix for kmeans (steps 4-5) */
	/*=========================================*/
	/* transpose eigenvectors matrix */
	/*-------------------------------*/
	u_mtx_mem = calloc(sizeof(F_TYPE), k*dp_num);
	assert(u_mtx_mem != NULL);
	u_mtx = calloc(sizeof(F_TYPE*), dp_num);
	assert(u_mtx != NULL);
	for (i = 0; i < dp_num; ++i)
		u_mtx[i] = u_mtx_mem + (i*k);

	transpose_matrix(k_eigenvectors, u_mtx, k, dp_num);


	/* normalize u_matrix by rows */
	/*----------------------------*/
	t_mtx_mem = calloc(sizeof(F_TYPE), k*dp_num);
	assert(t_mtx_mem != NULL);
	t_mtx = calloc(sizeof(F_TYPE*), dp_num);
	assert(t_mtx != NULL);
	for (i = 0; i < dp_num; ++i)
		t_mtx[i] = t_mtx_mem + (i*k);

	dummy_vector = calloc(sizeof(F_TYPE*), k);
	assert(dummy_vector != NULL);

	for (i = 0; i < k; ++i)
	{
		for (j = 0; j < dp_num; ++j)
			t_mtx[i][j] = u_mtx[i][j] / calc_l2_norm(u_mtx[i], dummy_vector, k);
	}
	

	/*==================*/
	/* kmeans algorithm */
	/*==================*/
	/* allocate memory for output centroids */
	d = k;
	n = dp_num;
	output_centroids_mem = calloc(sizeof(F_TYPE), k*d);
	assert(output_centroids_mem != NULL);
	output_centroids = calloc(sizeof(F_TYPE*), k);
	assert(output_centroids != NULL);
	for (i = 0; i < k; ++i)
		output_centroids[i] = output_centroids_mem + (i*d);

	/* allocate memory for datapoint assignment to cluster */
	output_cluster_assign = calloc(sizeof(int), n);
	assert(output_cluster_assign != NULL);

	
	/* THE KMEANS PROCEDURE CALL */
	status = kmeans(
		t_mtx, n, d, 
		output_centroids, output_cluster_assign,
		k, max_iter);
	/* TODO: handle kmeans status */

	/* print or return the output centroids */
	print_matrix(output_centroids, k, dim);

	free(u_mtx); free(u_mtx_mem);
	free(t_mtx); free(t_mtx_mem);

	free(output_centroids); free(output_centroids_mem);

	free(output_cluster_assign);

	return status;
}