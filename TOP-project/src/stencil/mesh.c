#include "stencil/mesh.h"

#include "logging.h"

#include <assert.h>
#include <stdlib.h>

mesh_t mesh_new(usz dim_x, usz dim_y, usz dim_z, mesh_kind_t kind) {
    usz const ghost_size = 2 * STENCIL_ORDER;

    f64* values = malloc((dim_x + ghost_size) * (dim_y + ghost_size) * (dim_z + ghost_size) * sizeof(f64));
    if(NULL == values ){
        error("failed to allocate memory for values of mesh %zu bytes",(dim_x + ghost_size) * (dim_y + ghost_size) * (dim_z + ghost_size) * sizeof(f64) );
    }

    cell_kind_t* kinds = malloc((dim_x + ghost_size) * (dim_y + ghost_size) * (dim_z + ghost_size) * sizeof(f64));
    if(kinds == NULL){
        free(values) ;
        error("failed to allocate memory for kinds of mesh %zu bytes ",(dim_x + ghost_size) * (dim_y + ghost_size) * (dim_z + ghost_size) * sizeof(f64));
    }

    mesh_t mesh ;
    mesh.dim_x = dim_x + ghost_size ;
    mesh.dim_y = dim_y + ghost_size ;
    mesh.dim_z = dim_z + ghost_size ;
    mesh.kind = kind ;
    mesh.cells.value = values ;
    mesh.cells.kind = kinds ;
    return mesh ;

}

void mesh_drop(mesh_t* self) {
    if (NULL != self->cells.value) {
        free(self->cells.value);
    }
    if (NULL != self->cells.kind) {
        free(self->cells.kind) ;
    }
}

/*
static char const* mesh_kind_as_str(mesh_t const* self) {
    static char const* MESH_KINDS_STR[] = {
        "CONSTANT",
        "INPUT",
        "OUTPUT",
    };
    return MESH_KINDS_STR[(usz)self->kind];
}
*/



cell_kind_t mesh_set_cell_kind(mesh_t const* self, usz i, usz j, usz k) {
    if ((i >= STENCIL_ORDER && i < self->dim_x - STENCIL_ORDER) &&
        (j >= STENCIL_ORDER && j < self->dim_y - STENCIL_ORDER) &&
        (k >= STENCIL_ORDER && k < self->dim_z - STENCIL_ORDER))
    {
        return CELL_KIND_CORE;
    } else {
        return CELL_KIND_PHANTOM;
    }
}

void mesh_copy_core(mesh_t* dst, mesh_t const* src) {
    assert(dst->dim_x == src->dim_x);
    assert(dst->dim_y == src->dim_y);
    assert(dst->dim_z == src->dim_z);

    for (usz k = STENCIL_ORDER; k < dst->dim_z - STENCIL_ORDER; ++k) {
        for (usz j = STENCIL_ORDER; j < dst->dim_y - STENCIL_ORDER; ++j) {
            for (usz i = STENCIL_ORDER; i < dst->dim_x - STENCIL_ORDER; ++i) {
                usz idx = i * dst->dim_y * dst->dim_z + j * dst->dim_z + k;
                assert(dst->cells.kind[idx] == CELL_KIND_CORE);
                assert(src->cells.kind[idx] == CELL_KIND_CORE);
                dst->cells.value[idx] = src->cells.value[idx];
            }
        }
    }
}

