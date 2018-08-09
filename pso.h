
/*===========Defining function pointers =========*/
typedef double (*fptr)(double*,int);

/*=========== Declaring external global variables =========*/
extern double * GBEST_MPSO_CCD;
extern int ND_MPSO_CCD;
extern int INDEX_MPSO_CCD;
extern fptr FUN_MPSO_CCD;

/*=========== Objective function selections ================*/
fptr functionLookup(char* fName);
double fun_1D(double t);

/*=========== Supportive functions for PSO =====================*/
void update_best(double * g, double * x, int nd);
int lower_bound_check(double g, double l, double u, double epsilon);
int upper_bound_check(double g, double l, double u, double epsilon);

/*=========== Optimization algorithms ========================*/
void pso(fptr fn, int nd, int np, int ni, double* lb, double* ub, double* value, double* gbest);
double brent_algorithm(double a, double b, double (*f) (double t), double tol);
void  ccd_local_optimizer(fptr fn, int nd, double* lb, double* ub, double* value, double* x);
void mpso_ccd(fptr fn, int nd, int np, int ni, int maxRun, double* lb, double* ub, double* value, double* gbest);
