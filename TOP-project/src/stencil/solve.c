#include "stencil/solve.h"

#include <assert.h>
#include <math.h>


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
        base *= base;  // Square the base
        exponent /= 2; // Divide the exponent by 2
    }

    return result;
}




/*
void solve_jacobi(mesh_t* A, mesh_t const* B, mesh_t* C) {
    assert(A->dim_x == B->dim_x && B->dim_x == C->dim_x);
    assert(A->dim_y == B->dim_y && B->dim_y == C->dim_y);
    assert(A->dim_z == B->dim_z && B->dim_z == C->dim_z);

    usz const dim_x = A->dim_x;
    usz const dim_y = A->dim_y;
    usz const dim_z = A->dim_z;

    for(usz kk = STENCIL_ORDER; kk < dim_z - STENCIL_ORDER ; kk+=BLOCK_SIZE_Z ){
        for(usz jj = STENCIL_ORDER; jj < dim_y - STENCIL_ORDER ; jj += BLOCK_SIZE_Y){
            for(usz ii = STENCIL_ORDER; ii < dim_x - STENCIL_ORDER ; ii += BLOCK_SIZE_X){

                for (usz k = kk; k < min(dim_z - STENCIL_ORDER, kk + BLOCK_SIZE_Z); ++k) {
                    for (usz j =  jj ; j < min(dim_y - STENCIL_ORDER,jj + BLOCK_SIZE_Y); ++j) {
                        for (usz i = ii ; i < min(dim_x - STENCIL_ORDER, ii + BLOCK_SIZE_X); ++i) {
                            usz idx = i * dim_y * dim_z + j * dim_z + k;
                            C->cells.value[idx] = A->cells.value[idx] * B->cells.value[idx];

                            for (usz o = 1; o <= STENCIL_ORDER; ++o) {
                                C->cells.value[idx] += A->cells.value[(i + o) * dim_y * dim_z + j * dim_z + k] *
                                                        B->cells.value[(i + o) * dim_y * dim_z + j * dim_z + k] / pow(17.0, (f64)o);
                                C->cells.value[idx] += A->cells.value[(i - o) * dim_y * dim_z + j * dim_z + k] *
                                                        B->cells.value[(i - o) * dim_y * dim_z + j * dim_z + k] / pow(17.0, (f64)o);
                                C->cells.value[idx] += A->cells.value[i * dim_y * dim_z + (j + o) * dim_z + k] *
                                                        B->cells.value[i * dim_y * dim_z + (j + o) * dim_z + k] / pow(17.0, (f64)o);
                                C->cells.value[idx] += A->cells.value[i * dim_y * dim_z + (j - o) * dim_z + k] *
                                                        B->cells.value[i * dim_y * dim_z + (j - o) * dim_z + k] / pow(17.0, (f64)o);
                                C->cells.value[idx] += A->cells.value[i * dim_y * dim_z + j * dim_z + (k + o)] *
                                                        B->cells.value[i * dim_y * dim_z + j * dim_z + (k + o)] / pow(17.0, (f64)o);
                                C->cells.value[idx] += A->cells.value[i * dim_y * dim_z + j * dim_z + (k - o)] *
                                                        B->cells.value[i * dim_y * dim_z + j * dim_z + (k - o)] / pow(17.0, (f64)o);
                            }
                        }
                    } 
                }
            }    
        }        
    }
    mesh_copy_core(A, C);
}



*/

void solve_jacobi(mesh_t* A, mesh_t const* B, mesh_t* C) {
    assert(A->dim_x == B->dim_x && B->dim_x == C->dim_x);
    assert(A->dim_y == B->dim_y && B->dim_y == C->dim_y);
    assert(A->dim_z == B->dim_z && B->dim_z == C->dim_z);

    usz const dim_x = A->dim_x;
    usz const dim_y = A->dim_y;
    usz const dim_z = A->dim_z;

    for(usz kk = STENCIL_ORDER; kk < dim_z - STENCIL_ORDER ; kk+=BLOCK_SIZE_Z ){
        for(usz jj = STENCIL_ORDER; jj < dim_y - STENCIL_ORDER ; jj += BLOCK_SIZE_Y){
            for(usz ii = STENCIL_ORDER; ii < dim_x - STENCIL_ORDER ; ii += BLOCK_SIZE_X){

                for (usz k = kk; k < min(dim_z - STENCIL_ORDER, kk + BLOCK_SIZE_Z); ++k) {
                    for (usz j =  jj ; j < min(dim_y - STENCIL_ORDER,jj + BLOCK_SIZE_Y); ++j) {
                        for (usz i = ii ; i < min(dim_x - STENCIL_ORDER, ii + BLOCK_SIZE_X); ++i) {
                            usz idx = i * dim_y * dim_z + j * dim_z + k;
                            C->cells.value[idx] = A->cells.value[idx] * B->cells.value[idx];

                            for (usz o = 1; o <= STENCIL_ORDER; ++o) {
                                C->cells.value[idx] += A->cells.value[(i + o) * dim_y * dim_z + j * dim_z + k] *
                                                        B->cells.value[(i + o) * dim_y * dim_z + j * dim_z + k] / power_of_17(o);
                                C->cells.value[idx] += A->cells.value[(i - o) * dim_y * dim_z + j * dim_z + k] *
                                                        B->cells.value[(i - o) * dim_y * dim_z + j * dim_z + k] /  power_of_17(o);
                                C->cells.value[idx] += A->cells.value[i * dim_y * dim_z + (j + o) * dim_z + k] *
                                                        B->cells.value[i * dim_y * dim_z + (j + o) * dim_z + k] /  power_of_17(o);
                                C->cells.value[idx] += A->cells.value[i * dim_y * dim_z + (j - o) * dim_z + k] *
                                                        B->cells.value[i * dim_y * dim_z + (j - o) * dim_z + k] /  power_of_17(o);
                                C->cells.value[idx] += A->cells.value[i * dim_y * dim_z + j * dim_z + (k + o)] *
                                                        B->cells.value[i * dim_y * dim_z + j * dim_z + (k + o)] /  power_of_17(o);
                                C->cells.value[idx] += A->cells.value[i * dim_y * dim_z + j * dim_z + (k - o)] *
                                                        B->cells.value[i * dim_y * dim_z + j * dim_z + (k - o)] /  power_of_17(o);
                            }
                        }
                    } 
                }
            }    
        }        
    }
    mesh_copy_core(A, C);
}
*/






/*
void solve_jacobi(mesh_t* A, mesh_t const* B, mesh_t* C) {
    assert(A->dim_x == B->dim_x && B->dim_x == C->dim_x);
    assert(A->dim_y == B->dim_y && B->dim_y == C->dim_y);
    assert(A->dim_z == B->dim_z && B->dim_z == C->dim_z);

    usz const dim_x = A->dim_x;
    usz const dim_y = A->dim_y;
    usz const dim_z = A->dim_z;
    for (usz k = STENCIL_ORDER; k < dim_z - STENCIL_ORDER; ++k) {
        for (usz j = STENCIL_ORDER; j < dim_y - STENCIL_ORDER; ++j) {
            for (usz i = STENCIL_ORDER; i < dim_x - STENCIL_ORDER; ++i) {
                C->cells[i][j][k].value = A->cells[i][j][k].value * B->cells[i][j][k].value;

                for (usz o = 1; o <= STENCIL_ORDER; ++o) {
                    C->cells[i][j][k].value += A->cells[i + o][j][k].value *
                                               B->cells[i + o][j][k].value / pow(17.0, (f64)o);
                    C->cells[i][j][k].value += A->cells[i - o][j][k].value *
                                               B->cells[i - o][j][k].value / pow(17.0, (f64)o);
                    C->cells[i][j][k].value += A->cells[i][j + o][k].value *
                                               B->cells[i][j + o][k].value / pow(17.0, (f64)o);
                    C->cells[i][j][k].value += A->cells[i][j - o][k].value *
                                               B->cells[i][j - o][k].value / pow(17.0, (f64)o);
                    C->cells[i][j][k].value += A->cells[i][j][k + o].value *
                                               B->cells[i][j][k + o].value / pow(17.0, (f64)o);
                    C->cells[i][j][k].value += A->cells[i][j][k - o].value *
                                               B->cells[i][j][k - o].value / pow(17.0, (f64)o);
                }
            }
        }
    }

    mesh_copy_core(A, C);
}

*/