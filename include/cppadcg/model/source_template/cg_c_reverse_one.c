//#define MODEL_FUNCTION
//#define BASE_TYPE
//#define FUNCTION_SPARSE_REVERSE_ONE 
//#define FUNCTION_REVERSE_ONE_SPARSITY
//#define M_SIZE
//#define N_SIZE
#include <stdlib.h>

void FUNCTION_SPARSE_REVERSE_ONE(unsigned long pos, BASE_TYPE const *const * in, BASE_TYPE * const * out);
void FUNCTION_REVERSE_ONE_SPARSITY(unsigned long pos, unsigned long const** elements, unsigned long* nnz);

int MODEL_FUNCTION(BASE_TYPE const x[], BASE_TYPE const ty[], BASE_TYPE px[], BASE_TYPE const py[]) {
    unsigned long e, i, ii, j, nnz;
    unsigned long const* pos;
    BASE_TYPE const * in[1];
    BASE_TYPE * out[1];
    BASE_TYPE* compressed;
    int found;

    found = 0;
    for (ii = 0; ii < M_SIZE; ii++) {
        if (py[ii] != 0.0) {
            if (found || py[ii] != 1.0) {
                return 1; // error 
            }
            i = ii;
            found = 1;
        }
    }
    for (j = 0; j < N_SISE; j++) {
        px[j] = 0;
    }
    if (!found) {
        return 0; //nothing to do
    }

    FUNCTION_REVERSE_ONE_SPARSITY(j, &pos, &nnz);

    compressed = (BASE_TYPE*) malloc(nnz * sizeof (BASE_TYPE));
    in[0] = x;
    out[0] = compressed;
    FUNCTION_SPARSE_REVERSE_ONE(i, in, out);

    for (e = 0; e < nnz; e++) {
        px[pos[e]] = compressed[e];
    }

    free(compressed);
    return 0;
};
