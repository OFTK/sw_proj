/*============*/
/* float type */
/*============*/
#define F_TYPE double
#define F_TYPE_FORMAT_SPEC ("%lf")
#define F_TYPE_OUTPUT_FORMAT_SPEC ("%.4f")
#define F_TYPE_MAX (DBL_MAX)
#define F_TYPE_ABS(x) (fabs(x))


/*===========*/
/* Constants */
/*===========*/

#define JACOBI_CONVERGENCE_SIGMA (0.001)
#define JACOBI_MAX_ITERATIONS (100)

#define KMEANS_MAX_ITERATIONS (300)

#define TRUE (1)
#define FALSE (0)


/*========*/
/* macros */
/*========*/
#define PRINT_ERROR() (printf("An Error Has Occured"))
#define PRINT_INVALID_INPUT() (printf("Invalid Input!"))


/*====================*/
/* enums and typedefs */
/*====================*/

enum status {
	Error=(-1),
	Success=0,
	Finish=1
};

enum goal {
	BAD=(-1),
	SPK=0,
	WAM=1,
	DDG=2,
	LNORM=3,
	JACOBI=4
};

typedef struct f_and_idx {
	F_TYPE f;
	int idx;
} f_and_idx;



/*==================*/
/* function headers */
/*==================*/

void print_matrix(F_TYPE** matrix, int n, int m);

enum status spkmeans_preperations(
	F_TYPE** input_dps, int dp_num, int dim,
	int* k, enum goal goal, 
	F_TYPE* o_mtx_mem, F_TYPE* o_eigenvalues);

enum status kmeans(
	F_TYPE** input_dps, int dp_num, int dim,
	F_TYPE** centroids, int* output_cluster_assign,
	int k, int max_iter);
