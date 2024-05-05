#include "stencil/solve.h"

#include <assert.h>
#include <math.h>
#include <omp.h>


#define BLOCK_SIZE_X 64
#define BLOCK_SIZE_Y 64
#define BLOCK_SIZE_Z 64

#define min(a, b) ((a) < (b) ? (a) : (b))


double power_of_17(usz exponent) {
    if (exponent == 0) return 1.0;
    if (exponent == 1) return 17.0;

    double result = 1.0;
    double base = 17.0;

    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result *= base;
        }
        base *= base;  
        exponent /= 2; 
    }

    return result;
}


void solve_jacobi(mesh_t* A, mesh_t const* B, mesh_t* C) {
    assert(A->dim_x == B->dim_x && B->dim_x == C->dim_x);
    assert(A->dim_y == B->dim_y && B->dim_y == C->dim_y);
    assert(A->dim_z == B->dim_z && B->dim_z == C->dim_z);

    usz const dim_x = A->dim_x;
    usz const dim_y = A->dim_y;
    usz const dim_z = A->dim_z;



    // Precompute
    f64 powers[STENCIL_ORDER];
    for (usz o = 1; o <= STENCIL_ORDER; ++o) {
        powers[o-1] = 1.0 / power_of_17(o);
    }
    usz dim_yz = dim_y * dim_z;
    usz dim_z_minus_ST = dim_z - STENCIL_ORDER;
    usz dim_x_minus_ST = dim_x - STENCIL_ORDER;
    usz dim_y_minus_ST = dim_y - STENCIL_ORDER;

    #pragma omp parallel for schedule(static,3)
    for (usz kk = STENCIL_ORDER; kk < dim_z_minus_ST; kk += BLOCK_SIZE_Z) {
        for (usz jj = STENCIL_ORDER; jj < dim_y_minus_ST; jj += BLOCK_SIZE_Y) {
            for (usz ii = STENCIL_ORDER; ii < dim_x_minus_ST; ii += BLOCK_SIZE_X) {
                
                
                for (usz i = ii; i < min(dim_x_minus_ST, ii + BLOCK_SIZE_X); ++i) {
                    for (usz j = jj; j < min(dim_y_minus_ST, jj + BLOCK_SIZE_Y); ++j) {
                        for (usz k = kk; k < min(dim_z_minus_ST, kk + BLOCK_SIZE_Z); ++k) {
                            usz idx = i * dim_yz + j * dim_z + k;
                            f64 sum = A->cells.value[idx] * B->cells.value[idx];
                            
                            for (usz o = 1; o <= STENCIL_ORDER; ++o) {
                                usz idx_plus_o = (i + o) * dim_yz + j * dim_z + k;
                                usz idx_minus_o = (i - o) * dim_yz + j * dim_z + k;
                                sum += (A->cells.value[idx_plus_o] * B->cells.value[idx_plus_o] +
                                        A->cells.value[idx_minus_o] * B->cells.value[idx_minus_o]) * powers[o - 1];

                                usz idx_plus_o_dimz = i * dim_yz + (j + o) * dim_z + k;
                                usz idx_minus_o_dimz = i * dim_yz + (j - o) * dim_z + k;
                                sum += (A->cells.value[idx_plus_o_dimz] * B->cells.value[idx_plus_o_dimz] +
                                        A->cells.value[idx_minus_o_dimz] * B->cells.value[idx_minus_o_dimz]) * powers[o - 1];

                                usz idx_plus_o_dimyz = i * dim_yz + j * dim_z + (k + o);
                                usz idx_minus_o_dimyz = i * dim_yz + j * dim_z + (k - o);
                                sum += (A->cells.value[idx_plus_o_dimyz] * B->cells.value[idx_plus_o_dimyz] +
                                        A->cells.value[idx_minus_o_dimyz] * B->cells.value[idx_minus_o_dimyz]) * powers[o - 1];
                            }
                            C->cells.value[idx] = sum;
                        }
                    }
                }
            }
        }
    }
    mesh_copy_core(A, C);
}
