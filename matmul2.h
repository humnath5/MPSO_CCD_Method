
struct matrix {

    double **  mat;
    int        row;
    int	       col;
};

struct matrix * get_matrix(int,int);
struct matrix * matrix_mult(struct matrix*,struct matrix*);

void copy(struct matrix**,double*,double);
struct matrix * transpose(struct matrix*);

double norm(struct matrix*,int);
void print_matrix(struct matrix *);
double ** get_mat(int,int);
void free_matrix(struct matrix*);


struct matrix* read_matrix(char* filename,int nrow,int ncol);
void write_matrix(char* filename,struct matrix *m);

