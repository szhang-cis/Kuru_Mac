#ifndef SPARSEASSEMBLYNATIVE_H
#define SPARSEASSEMBLYNATIVE_H

#include <algorithm>
#include <cstdint>
#include <cstdlib>

#ifndef LL_TYPES
#define LL_TYPES
using Real = double;
using Integer = std::int64_t;
using UInteger = std::uint64_t;
#endif


#ifndef BINARY_LOCATE_DEF_
#define BINARY_LOCATE_DEF_
template<class RandomIt, class T>
inline RandomIt binary_locate(RandomIt first, RandomIt last, const T& val) {
  if(val == *first) return first;
  auto d = std::distance(first, last);
  if(d==1) return first;
  auto center = (first + (d/2));
  if(val < *center) return binary_locate(first, center, val);
  return binary_locate(center, last, val);
}
#endif




inline void SparseAssemblyNativeCSR_(
    const double *coeff,
    const int *data_local_indices,
    const int *data_global_indices,
    int elem,
    int local_capacity,
    double *data
    ) {

    // FILL DATA VALUES FOR CSR
    for (int i=0; i<local_capacity; ++i) {
        data[data_global_indices[local_capacity*elem+i]] += coeff[data_local_indices[local_capacity*elem+i]];
    }
}



inline void SparseAssemblyNativeCSR_RecomputeDataIndex_(
    const double *coeff,
    int *indices,
    int *indptr,
    double *data,
    int elem,
    int nvar,
    int nodeperelem,
    const UInteger *elements,
    const Integer *sorter) {

    int ndof = nvar*nodeperelem;
    int local_capacity = ndof*ndof;

    int *current_row_column = (int*)malloc(sizeof(int)*ndof);
    int *full_current_column = (int*)malloc(sizeof(int)*local_capacity);

    for (int counter=0; counter<nodeperelem; ++ counter) {
        const int const_elem_retriever = nvar*elements[elem*nodeperelem+counter];
        for (int ncounter=0; ncounter<nvar; ++ncounter) {
            current_row_column[nvar*counter+ncounter] = const_elem_retriever+ncounter;
        }
    }

    int *current_row_column_local = (int*)malloc(sizeof(int)*ndof);
    int *full_current_column_local = (int*)malloc(sizeof(int)*local_capacity);
    for (int counter=0; counter<nodeperelem; ++ counter) {
        const int node = sorter[elem*nodeperelem+counter];
        for (int ncounter=0; ncounter<nvar; ++ncounter) {
            current_row_column_local[nvar*counter+ncounter] = node*nvar+ncounter;
        }
    }

    for (int i=0; i<ndof; ++i) {
        const int current_local_ndof = current_row_column_local[i]*ndof;
        const int current_global_ndof = current_row_column[i];
        const int current_global_row = indptr[current_global_ndof];
        const int nnz = indptr[current_global_ndof+1] - current_global_row;
        int *search_space = (int*)malloc(sizeof(int)*nnz);
        for (int k=0; k<nnz; ++k) {
            search_space[k] = indices[current_global_row+k];
        }

        for (int j=0; j<ndof; ++j) {
            // int Iterr = std::find(search_space,search_space+nnz,current_row_column[j]) - search_space;
            int Iterr = binary_locate(search_space,search_space+nnz,current_row_column[j]) - search_space;
            full_current_column[i*ndof+j] = current_global_row + Iterr;
            full_current_column_local[i*ndof+j] = current_local_ndof+current_row_column_local[j];
        }
        free(search_space);
    }

    // FILL DATA VALUES FOR CSR
    for (int i=0; i<local_capacity; ++i) {
        data[full_current_column[i]] += coeff[full_current_column_local[i]];
    }

    free(current_row_column);
    free(full_current_column);

    free(current_row_column_local);
    free(full_current_column_local);
}


#endif
