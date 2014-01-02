//#define MODEL_FUNCTION
//#define BASE_TYPE
//#define FUNCTION_SPARSE_REVERSE_TWO
//#define FUNCTION_REVERSE_TWO_SPARSITY
//#define M_SIZE
//#define N_SIZE
#include <stdlib.h>

void FUNCTION_SPARSE_REVERSE_TWO(unsigned long pos, BASE_TYPE const *const * in, BASE_TYPE * const * out);
void FUNCTION_REVERSE_TWO_SPARSITY(unsigned long pos, unsigned long const** elements, unsigned long* nnz);

int MODEL_FUNCTION(BASE_TYPE const tx[], BASE_TYPE const ty[], BASE_TYPE px[], BASE_TYPE const py[]) {
    unsigned long e, i, j, jj, nnz;
    unsigned long const* pos;
    BASE_TYPE const * in[2];
    BASE_TYPE * out[1];
    BASE_TYPE x[N_SIZE];
    BASE_TYPE* compressed;
    int found;

    found = 0;
    for (jj = 0; jj < N_SIZE; jj++) {
        if (tx[jj * 2 + 1] != 0.0) {
            if (found || tx[jj * 2 + 1] != 1.0) {
                return 1; // error 
            }
            j = jj;
            found = 1;
        }
    }
    for (jj = 0; jj < N_SISE; jj++) {
        px[jj * 2 + 1] = 0;
    }
    if (!found) {
        return 0; //nothing to do
    }

    for (jj = 0; jj < N_SIZE; jj++)
        x[jj] = tx[jj * 2];

    FUNCTION_REVERSE_TWO_SPARSITY(j, &pos, &nnz);

    compressed = (BASE_TYPE*) malloc(nnz * sizeof (BASE_TYPE));
    in[0] = x;
    in[1] = py; // expected size is (k+1)*m
    out[0] = compressed;
    FUNCTION_SPARSE_REVERSE_TWO(i, in, out);

    for (e = 0; e < nnz; e++) {
        px[pos[e] * 2 + 1] = compressed[e];
    }

    free(compressed);
    return 0;
};
