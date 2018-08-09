
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "matmul2.h"

 struct  matrix * get_matrix(int row, int col){

      struct  matrix *a;

      a = (struct matrix *) calloc(1,sizeof(struct matrix));
     if( row <1 ) {
     	printf("row index less than 1\n");
     	return NULL;}
     if( col <1 ) {
     	printf("column index less than 1\n");
     	return NULL;}
     a->row = row;
     a->col = col;
     a->mat = get_mat(row, col);
     return a;
     }




 struct matrix * matrix_mult(struct matrix *m1, struct matrix *m2)
 {
 int i, j, k;
 struct matrix  * ans;
 /* error checking!!*/

 //printf("m1:(%d x %d)\tm2:(%d x %d)\n",m1->row,m1->col,m2->row,m2->col);

 ans = get_matrix(m1->row, m2->col);

 for(i=0 ; i< m1->row; i++){
 	for(j=0; j<m2->col; j++){
 		for(k=0; k<m1->col; k++)
 		  (ans->mat)[i][j]+=((m1->mat)[i][k])*((m2->mat)[k][j]);}}



 return ans;
 }


struct matrix * transpose(struct matrix * m)
{
	struct matrix * ans;
	ans = get_matrix(m->col,m->row);

	int i, j;
	for(i=0;i<m->row;i++)
	{
		for(j=0;j<m->col;j++)
		{
			ans->mat[j][i] = m -> mat[i][j];
		}
	}
	return ans;
}




double norm(struct matrix * m, int p)
{
	int i,j;
	double sum = 0;
	for(i=0;i<m->row;i++)
	{
		for(j=0;j<m->col;j++)
		{
			sum += pow((m->mat[i][j]),p);
		}
	}

	double ans = pow(sum,1.0/(double)p);
	return ans;
}
struct matrix* read_matrix(char* filename,int nrow,int ncol){
    int i, j;
    FILE* fh;
    struct matrix* data;

    data = get_matrix(nrow, ncol);

    if ((fh = fopen(filename, "r")) == NULL){
        fprintf(fh,"Error: Cannot open %s\n", filename);
        exit(1);
    }


    for( i=0;i< data->row;i++)

	{
		for(j=0;j< data->col;j++)
		{
		  fscanf(fh, "%lf,", &(data->mat)[i][j]);
		}
	}
    fclose(fh);
    return data;
}


double ** get_mat(int n1, int n2)
{
/*  Purpose: This function allocates space for a n1 by n2 doubly
             indexed array of doubles.

    Return Value:  n1 by n2 matrix implemented as pointer to pointer
                   to double.  The indexing starts at [0][0] to
                   [n1-1][n2-1].  And the storage is allocated contiguously.



    Arguments:

        n1              row index (input)
        n2              column index (input)


        There are no Fatal Errors.

    Related Functions:  free_mat

 */
    int i;
    double ** mat, * temp_ptr;

            /*  Allocate space for the array  */
    temp_ptr   = (double *) calloc(n1*n2, sizeof(double));
    if((void *)temp_ptr == NULL){
       /* *inform = 4; */
        return NULL;
    }

    mat = (double **) calloc(n1, sizeof(double *));
    if((void *)(mat) == NULL){
        /**inform = 4;*/
        return NULL;
    }

    for(i=0; i< n1; i++)
        mat[i] = &(temp_ptr[i*n2]);

    return mat;
}


void free_matrix(struct matrix * m)
{
	free(*(m->mat));
	free(m->mat);
	free(m);
}

