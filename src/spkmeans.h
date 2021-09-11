#define F_TYPE double

/*======*/
/* misc */
/*======*/
#define PRINT_ERROR() (printf("An Error Has Occured"))
#define PRINT_INVALID_INPUT() (printf("Invalid Input!"))

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

void print_matrix(F_TYPE** matrix, int n, int m);

enum status spkmeans_preperations(
	F_TYPE** input_dps, int dp_num, int dim,
	int* k, enum goal goal, 
	F_TYPE* o_mtx_mem, F_TYPE* o_eigenvalues);

enum status kmeans(
	F_TYPE** input_dps, int dp_num, int dim,
	F_TYPE** output_centrds, int* output_cluster_assign,
	int k, int max_iter);