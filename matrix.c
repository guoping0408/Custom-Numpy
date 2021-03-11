#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/

/*
 * Generates a random double between `low` and `high`.
 */
double rand_double(double low, double high) {
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/*
 * Generates a random matrix with `seed`.
 */
void rand_matrix(matrix *result, unsigned int seed, double low, double high) {
    srand(seed);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Allocate space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. Remember to set all fieds of the matrix struct.
 * `parent` should be set to NULL to indicate that this matrix is not a slice.
 * You should return -1 if either `rows` or `cols` or both have invalid values, or if any
 * call to allocate memory in this function fails. If you don't set python error messages here upon
 * failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix(matrix **mat, int rows, int cols) {
    /* TODO: YOUR CODE HERE */

    // return -1 if either `rows` or `cols` or both have invalid values
    if (rows <= 0 || cols <= 0) {
        PyErr_SetString(PyExc_ValueError, "Invalid rows or cols");
        return -1;
    }

    *mat = malloc(sizeof(matrix));
    // return -1 if any call to allocate memory fails
    if (*mat == NULL){
        PyErr_SetString(PyExc_RuntimeError, "malloc for matrix* fails");
        return -1;
    } 

    matrix* toAlloc = *mat;
    toAlloc->rows = rows;
    toAlloc->cols = cols;
    toAlloc->ref_cnt = 1;
    toAlloc->ref_cnt_ptr = &toAlloc->ref_cnt;
    toAlloc->parent = NULL;
    if (cols == 1 || rows == 1) {
        toAlloc->is_1d = 1;
    } else {
        toAlloc->is_1d = 0;
    }

    double* whole_matrix = calloc(cols*rows, sizeof(double));
    if (whole_matrix == NULL) {
        free(*mat);
        PyErr_SetString(PyExc_RuntimeError, "malloc for matrix*'s data fails");
        return -1;
    }
    
    // This field keeps track of the start address of the entire matrix
    // that we calloc all at once
    // later in deallocate_matrix, we just have to call free(mat->whole_matrix) once :)
    toAlloc->whole_matrix = whole_matrix;

    toAlloc->data = malloc(sizeof(double*) * rows);
    if (toAlloc->data == NULL) {
        free(*mat);
        free(whole_matrix);
        PyErr_SetString(PyExc_RuntimeError, "malloc for matrix*'s data fails");
        return -1;
    }
    for(int r = 0; r < rows; r++) {
        toAlloc->data[r] = whole_matrix + r*cols;
        // TODO: check whole_matrix + 1 would add sizeof(double) bytes to the address of whole_matrix?
    }
    return 0;
}

/*
 * Allocate space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * This is equivalent to setting the new matrix to be
 * from[row_offset:row_offset + rows, col_offset:col_offset + cols]
 * If you don't set python error messages here upon failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int row_offset, int col_offset,
                        int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    // printf("Doing allocate_matrix_ref\n");
    if (row_offset < 0 || col_offset < 0|| rows < 0 || cols < 0) {
        PyErr_SetString(PyExc_ValueError, "From allocate_matrix_ref: Invalid rows or cols");
        return -1;
    }

    // printf("Checking...\n");
    // TODO: check row_offset(&row_offset+rows) and col_offset(&col_offset+cols)
    // are less than from->rows and from->cols, respectively
    
    if ((row_offset + rows) > from->rows || (col_offset + cols) > from->cols) {
        PyErr_SetString(PyExc_RuntimeError, "From allocate_matrix_ref: Invalid rows or cols");
        return -1; 
    }
    

    *mat = malloc(sizeof(matrix));
    // return -1 if any call to allocate memory fails
    if (*mat == NULL){
        PyErr_SetString(PyExc_RuntimeError, "From allocate_matrix_ref: malloc for matrix* fails");
        return -1;
    } 

    // printf("Allocating toAlloc->data\n");
    matrix* toAlloc = *mat;
    toAlloc->data = malloc(sizeof(double*) * rows);
    if (toAlloc->data == NULL) {
        free(*mat);
        PyErr_SetString(PyExc_RuntimeError, "From allocate_matrix_ref: malloc for matrix*'s data fails");
        return -1;
    }
    toAlloc->rows = rows;
    toAlloc->cols = cols;
    toAlloc->whole_matrix = from->whole_matrix;
    toAlloc->ref_cnt_ptr = from->ref_cnt_ptr;
    
    toAlloc->parent = from;
    if (cols == 1 || rows == 1) {
        toAlloc->is_1d = 1;
    } else {
        toAlloc->is_1d = 0;
    }
    // printf("col_offset is %d\n", col_offset);
    // printf("row_offset is %d\n", row_offset);
    // printf("The from matrix is \n");
    
    // for(int i = 0; i < from->rows; i++) {
    //     for(int j = 0; j < from->cols; j++) {
    //         printf("%f  ", from->data[i][j]);
    //     }
    //     printf("\n");
    // }

    for(int r = 0; r < rows; r++) {
        // from->data[row_offset+r][col_offset];
        // toAlloc->data[r] = from->whole_matrix + (r+row_offset)*from->cols + col_offset;
        toAlloc->data[r] = &from->data[row_offset+r][col_offset];
    }

    // printf("ref_allocated matrix is \n");
    // for(int r = 0; r < rows; r++) {
    //     for(int c = 0; c < cols; c++) {
    //         printf("%f ", toAlloc->data[r][c]);
    //     }
    //     printf("\n");
    // }
    *(from->ref_cnt_ptr) += 1;
    toAlloc->ref_cnt = *(from->ref_cnt_ptr);
    return 0;

    
}

/*
 * This function will be called automatically by Python when a numc matrix loses all of its
 * reference pointers.
 * You need to make sure that you only free `mat->data` if no other existing matrices are also
 * referring this data array.
 * See the spec for more information.
 */
void deallocate_matrix(matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if (mat == NULL) {
        return;
    }
    
    *(mat->ref_cnt_ptr) -= 1;

    if (*(mat->ref_cnt_ptr) > 0) {
        free(mat->data); 
        free(mat);
        return;
    }

    
    // TODO : We calloc a big chuch of whole_matrix of size rows * cols at once
    // so when freeing mat->data[0], we essentially free the entire whole_matrix
    // free(mat->data[0]);
    // printf("WE ARE DEALLOCATING!!!!\n");
    free(mat->whole_matrix);
    free(mat->data);
    free(mat);
}

/*
 * Return the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid.
 */
double get(matrix *mat, int row, int col) {
    /* TODO: YOUR CODE HERE */
    // to Incorporate slicing
    if (mat->parent != NULL) {
        return mat->data[row][col];
    }
    else {
        return mat->whole_matrix[row*mat->cols+col];
    }
    
}

/*
 * Set the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val) {
    /* TODO: YOUR CODE HERE */
    if (mat->parent != NULL) {
        mat->data[row][col] = val;
    } else {
        mat->whole_matrix[row*mat->cols+col] = val;
    }

    // printf("Set mat row %d col %d with value %d\n", mat->rows, mat->cols, mat->data[row][col]);
}

/*
 * Set all entries in mat to val
 */
void fill_matrix(matrix *mat, double val) {
    /* TODO: YOUR CODE HERE */
    
    if (mat->parent == NULL) {
        // For the root matrix, parallelize using whole_matrix (1d chuck of memory)
        int size = mat->rows * mat->cols;
        __m256d v = _mm256_set1_pd (val);
        // memset
        #pragma omp parallel for
        for(int i = 0; i < size/16*16; i+=16) {
            _mm256_storeu_pd (mat->whole_matrix+i, v);
            _mm256_storeu_pd (mat->whole_matrix+i+4, v);
            _mm256_storeu_pd (mat->whole_matrix+i+8, v);
            _mm256_storeu_pd (mat->whole_matrix+i+12, v);
            // mat->whole_matrix[i] = val;
        }
        

        // Tailing case 
        for(int i = size/16*16; i<size; i++) {
             mat->whole_matrix[i] = val;
        }
        return;
    }

    #pragma omp parallel for
    for(int r = 0; r < mat->rows; r++) {
        for(int c = 0; c < mat->cols; c++) {
            mat->data[r][c] = val;
        }
    }
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */

    if (!(mat1->rows == mat2->rows && mat1->cols == mat2->cols )) {
        PyErr_SetString(PyExc_ValueError, "Error when add: shape not matched");
        return -1;
    }
    
    // Allocate result matrix if it's not allocated
    if (result == NULL) {
        int success = allocate_matrix(&result, mat1->rows, mat1->cols);
        if (success != 0) {
            return success;
        }
    } else {
        // check if result matrix has the same shape as mat1&mat2
        if (!(result->rows == mat1->rows && result->cols == mat1->cols)) {
            PyErr_SetString(PyExc_ValueError, "Error when add: shape not matched");
            return -1;
        }
    }

    // Speed up for non-slicing matrix
    if (mat1->parent == NULL && mat2->parent == NULL) {
        // Hive machines have 4 cores with two hyperthreads each
        // Can do up to 8 software threads

        int size = mat1->rows * mat2->cols;
        int unrolling = 4;
        int SIMD = 4;
        int step = unrolling * SIMD;
        int tail = size / step * step;
        #pragma omp parallel for
        for(int i = 0; i < tail; i += step) {
            // UNROLL 1
            
            __m256d add1 =  _mm256_add_pd(_mm256_loadu_pd(mat1->whole_matrix + i), 
                                          _mm256_loadu_pd(mat2->whole_matrix + i));
            
            // Store 256-bits (composed of 4 packed double-precision (64-bit) floating-point elements) 
            // from add1 into memory pointed at result->whole_matrix+i.
            _mm256_storeu_pd(result->whole_matrix+i, add1);
            
            // UNROLL 2 
            
            // starting from i+4 becase we deal with 4 doubles at once in unroll1 already
            // So, we jump to i + 4 then do add again
            __m256d add2 =  _mm256_add_pd(_mm256_loadu_pd(mat1->whole_matrix + i + 4), 
                                          _mm256_loadu_pd(mat2->whole_matrix + i + 4));
            _mm256_storeu_pd(result->whole_matrix+i+4, add2);


            // UNROLL 3
            __m256d add3 =  _mm256_add_pd(_mm256_loadu_pd(mat1->whole_matrix + i + 8), 
                                          _mm256_loadu_pd(mat2->whole_matrix + i + 8));
            _mm256_storeu_pd(result->whole_matrix+i+8, add3);


            // UNROLL 4
            __m256d add4 =  _mm256_add_pd(_mm256_loadu_pd(mat1->whole_matrix + i + 12), 
                                          _mm256_loadu_pd(mat2->whole_matrix + i + 12));
            _mm256_storeu_pd(result->whole_matrix+i+12, add4);
        }

        // deal with the tailing
        // do ONE double at a time since we don't wanna worry about
        // another tailing case for SIMD
        for(int i = tail; i < size; i++) {
            result->whole_matrix[i] = mat1->whole_matrix[i] + mat2->whole_matrix[i];
        }
        return 0;
    }

    #pragma omp parallel for
    for(int r = 0; r < result->rows; r++) {
        for(int c = 0; c < result->cols; c++) {
            result->data[r][c] = mat1->data[r][c] + mat2->data[r][c];
        }
    }
    return 0;
}

/*
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if (!(mat1->rows == mat2->rows && mat1->cols == mat2->cols )) {
        PyErr_SetString(PyExc_ValueError, "Error when add: shape not matched");
        return -1;
    }

    // Allocate result matrix if it's not allocated
    if (result == NULL) {
        int success = allocate_matrix(&result, mat1->rows, mat1->cols);
        if (success != 0) {
            return success;
        }
    } else {
        // check if result matrix has the same shape as mat1&mat2
        if (!(result->rows == mat1->rows && result->cols == mat1->cols)) {
            PyErr_SetString(PyExc_ValueError, "Error when add: shape not matched");
            return -1;
        }
    }
    // Speed up for non-slicing matrix
    if (mat1->parent == NULL && mat2->parent == NULL) {
        // Hive machines have 4 cores with two hyperthreads each
        // Can do up to 8 software threads

        int size = mat1->rows * mat2->cols;
        int unrolling = 4;
        int SIMD = 4;
        int step = unrolling * SIMD;
        int tail = size / step * step;
        #pragma omp parallel for
        for(int i = 0; i < tail; i += step) {
            // UNROLL 1
            __m256d add1 =  _mm256_sub_pd(_mm256_loadu_pd(mat1->whole_matrix + i), 
                                          _mm256_loadu_pd(mat2->whole_matrix + i));
            
            // Store 256-bits (composed of 4 packed double-precision (64-bit) floating-point elements) 
            // from add1 into memory pointed at result->whole_matrix+i.
            _mm256_storeu_pd(result->whole_matrix+i, add1);
            
            // UNROLL 2 
            // starting from i+4 becase we deal with 4 doubles at once in unroll1 already
            // So, we jump to i + 4 then do add again
            __m256d add2 =  _mm256_sub_pd(_mm256_loadu_pd(mat1->whole_matrix + i + 4), 
                                          _mm256_loadu_pd(mat2->whole_matrix + i + 4));
            _mm256_storeu_pd(result->whole_matrix+i+4, add2);


            // UNROLL 3
            __m256d add3 =  _mm256_sub_pd(_mm256_loadu_pd(mat1->whole_matrix + i + 8), 
                                          _mm256_loadu_pd(mat2->whole_matrix + i + 8));
            _mm256_storeu_pd(result->whole_matrix+i+8, add3);


            // UNROLL 4
            __m256d add4 =  _mm256_sub_pd(_mm256_loadu_pd(mat1->whole_matrix + i + 12), 
                                          _mm256_loadu_pd(mat2->whole_matrix + i + 12));
            _mm256_storeu_pd(result->whole_matrix+i+12, add4);
        }

        // deal with the tailing
        // do ONE double at a time since we don't wanna worry about
        // another tailing case for SIMD
        for(int i = tail; i < size; i++) {
            result->whole_matrix[i] = mat1->whole_matrix[i] - mat2->whole_matrix[i];
        }
        return 0;
    }

    #pragma omp parallel for
    for(int r = 0; r < result->rows; r++) {
        for(int c = 0; c < result->cols; c++) {
            result->data[r][c] = mat1->data[r][c] - mat2->data[r][c];
        }
    }
    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if (mat1->cols != mat2->rows) {
        PyErr_SetString(PyExc_ValueError, "Error when mul: unmatched shape");
        return -1;
    }
    // Allocate result matrix if it's not allocated
    if (result == NULL) {
        int success = allocate_matrix(&result, mat1->rows, mat2->cols);
        if (success != 0) {
            return success;
        }
    } 
    // check if result matrix has the correct shape
    if (!(result->rows == mat1->rows && result->cols == mat2->cols)) {
        PyErr_SetString(PyExc_ValueError, "Error when mul: shape not matched");
        return -1;
    }

    // Speed up for non-slicing matrix
    if (mat1->parent == NULL && mat2->parent == NULL) {
        // Since sometimes result and (mat1 || mat2) might be the same matrix
        // Create tmp matrix for holding the temporary result without affecting the operand
        // during the process of multiplication
        // matrix* tmp = NULL;
        // int success = allocate_matrix(&tmp, mat1->rows, mat2->cols);
        // if (success != 0) {
        //     return success;
        // }
        int unrolling = 4;
        int SIMD = 4;
        int step = unrolling * SIMD;

        int tail = mat2->cols / step * step;

        register __m256d m1_i_j_x4;

        register __m256d mat2rowj_1;
        register __m256d mat2rowj_2;
        register __m256d mat2rowj_3;
        register __m256d mat2rowj_4;

        register __m256d before1;
        register __m256d before2;
        register __m256d before3;
        register __m256d before4;

        register __m256d after1;
        register __m256d after2;
        register __m256d after3;
        register __m256d after4;


        #pragma omp parallel for
        for(int i = 0; i < mat1->rows; i++) {
            // #pragma omp parallel for
            // 3 - 2
            for(int j = 0; j < mat1->cols/4 * 4; j+=4) {
                // ************For every element of matrix1[i][j]**********
                // mat1[i][j]
                __m256d m1_i_j_x4 = _mm256_set1_pd(mat1->whole_matrix[i*mat1->cols+j]);

                // mat1[i][j+1]
                __m256d m1_i_j1_x4 = _mm256_set1_pd(mat1->whole_matrix[i*mat1->cols+j+1]);

                // mat1[i][j+2]
                __m256d m1_i_j2_x4 = _mm256_set1_pd(mat1->whole_matrix[i*mat1->cols+j+2]);

                // mat1[i][j+3]
                __m256d m1_i_j3_x4 = _mm256_set1_pd(mat1->whole_matrix[i*mat1->cols+j+3]);

                // // mat1[i][j+4]
                // __m256d m1_i_j4_x4 = _mm256_set1_pd(mat1->whole_matrix[i*mat1->cols+j+4]);

                // // mat1[i][j+5]
                // __m256d m1_i_j5_x4 = _mm256_set1_pd(mat1->whole_matrix[i*mat1->cols+j+5]);

                // // mat1[i][j+6]
                // __m256d m1_i_j6_x4 = _mm256_set1_pd(mat1->whole_matrix[i*mat1->cols+j+6]);

                // // mat1[i][j+7]
                // __m256d m1_i_j7_x4 = _mm256_set1_pd(mat1->whole_matrix[i*mat1->cols+j+7]);
                
                // We extract mat1[i][j], then times it with all the elements in row mat[j]
                // then store them in every element in the row result[i]
                // mat2->cols = tail
                
                for(int k = 0; k < tail; k += step) {
                    /* In the end of this loop
                     * result(tmp)[i][k:k+4] += m1[i][j]_x4 * m2[j][k:k+4]
                     * 
                     * So, result[i][j] will get updated as in the next iteration of j
                     * since k will always starts at 0 in every j
                     */
                    // mat2->whole_matrix+ j * result->cols + k == mat2->data[j][k: k+4]
                    mat2rowj_1 = _mm256_loadu_pd(mat2->whole_matrix+ j * result->cols + k);
                    mat2rowj_2 = _mm256_loadu_pd(mat2->whole_matrix+ j * result->cols + k + 4);
                    mat2rowj_3 = _mm256_loadu_pd(mat2->whole_matrix+ j * result->cols + k + 8);
                    mat2rowj_4 = _mm256_loadu_pd(mat2->whole_matrix+ j * result->cols + k + 12);

                    // j+1 ********************
                    __m256d mat2rowj1_1 = _mm256_loadu_pd(mat2->whole_matrix+ (j+1) * result->cols + k);
                    __m256d mat2rowj1_2 = _mm256_loadu_pd(mat2->whole_matrix+ (j+1) * result->cols + k + 4);
                    __m256d mat2rowj1_3 = _mm256_loadu_pd(mat2->whole_matrix+ (j+1) * result->cols + k + 8);
                    __m256d mat2rowj1_4 = _mm256_loadu_pd(mat2->whole_matrix+ (j+1) * result->cols + k + 12);

                    // j+2 ********************
                    __m256d mat2rowj2_1 = _mm256_loadu_pd(mat2->whole_matrix+ (j+2) * result->cols + k);
                    __m256d mat2rowj2_2 = _mm256_loadu_pd(mat2->whole_matrix+ (j+2) * result->cols + k + 4);
                    __m256d mat2rowj2_3 = _mm256_loadu_pd(mat2->whole_matrix+ (j+2) * result->cols + k + 8);
                    __m256d mat2rowj2_4 = _mm256_loadu_pd(mat2->whole_matrix+ (j+2) * result->cols + k + 12);

                    // j+3 ********************
                    __m256d mat2rowj3_1 = _mm256_loadu_pd(mat2->whole_matrix+ (j+3) * result->cols + k);
                    __m256d mat2rowj3_2 = _mm256_loadu_pd(mat2->whole_matrix+ (j+3) * result->cols + k + 4);
                    __m256d mat2rowj3_3 = _mm256_loadu_pd(mat2->whole_matrix+ (j+3) * result->cols + k + 8);
                    __m256d mat2rowj3_4 = _mm256_loadu_pd(mat2->whole_matrix+ (j+3) * result->cols + k + 12);

                    // // j+4 ********************
                    // __m256d mat2rowj4_1 = _mm256_loadu_pd(mat2->whole_matrix+ (j+4) * result->cols + k);
                    // __m256d mat2rowj4_2 = _mm256_loadu_pd(mat2->whole_matrix+ (j+4) * result->cols + k + 4);
                    // __m256d mat2rowj4_3 = _mm256_loadu_pd(mat2->whole_matrix+ (j+4) * result->cols + k + 8);
                    // __m256d mat2rowj4_4 = _mm256_loadu_pd(mat2->whole_matrix+ (j+4) * result->cols + k + 12);

                    // // j+5 ********************
                    // __m256d mat2rowj5_1 = _mm256_loadu_pd(mat2->whole_matrix+ (j+5) * result->cols + k);
                    // __m256d mat2rowj5_2 = _mm256_loadu_pd(mat2->whole_matrix+ (j+5) * result->cols + k + 4);
                    // __m256d mat2rowj5_3 = _mm256_loadu_pd(mat2->whole_matrix+ (j+5) * result->cols + k + 8);
                    // __m256d mat2rowj5_4 = _mm256_loadu_pd(mat2->whole_matrix+ (j+5) * result->cols + k + 12);

                    // // j+6 ********************
                    // __m256d mat2rowj6_1 = _mm256_loadu_pd(mat2->whole_matrix+ (j+6) * result->cols + k);
                    // __m256d mat2rowj6_2 = _mm256_loadu_pd(mat2->whole_matrix+ (j+6) * result->cols + k + 4);
                    // __m256d mat2rowj6_3 = _mm256_loadu_pd(mat2->whole_matrix+ (j+6) * result->cols + k + 8);
                    // __m256d mat2rowj6_4 = _mm256_loadu_pd(mat2->whole_matrix+ (j+6) * result->cols + k + 12);

                    // // j+7 ********************
                    // __m256d mat2rowj7_1 = _mm256_loadu_pd(mat2->whole_matrix+ (j+7) * result->cols + k);
                    // __m256d mat2rowj7_2 = _mm256_loadu_pd(mat2->whole_matrix+ (j+7) * result->cols + k + 4);
                    // __m256d mat2rowj7_3 = _mm256_loadu_pd(mat2->whole_matrix+ (j+7) * result->cols + k + 8);
                    // __m256d mat2rowj7_4 = _mm256_loadu_pd(mat2->whole_matrix+ (j+7) * result->cols + k + 12);
                    

                    
                    // tmp->whole_matrix+ i * result->cols + k == tmp->data[i][k:k+4]
                    before1 = _mm256_loadu_pd(result->whole_matrix+ i * result->cols + k);
                    before2 = _mm256_loadu_pd(result->whole_matrix+ i * result->cols + k + 4);
                    before3 = _mm256_loadu_pd(result->whole_matrix+ i * result->cols + k + 8);
                    before4 = _mm256_loadu_pd(result->whole_matrix+ i * result->cols + k + 12);

                    after1 = _mm256_fmadd_pd(m1_i_j_x4, mat2rowj_1, before1);
                    after2 = _mm256_fmadd_pd(m1_i_j_x4, mat2rowj_2, before2);
                    after3 = _mm256_fmadd_pd(m1_i_j_x4, mat2rowj_3, before3);
                    after4 = _mm256_fmadd_pd(m1_i_j_x4, mat2rowj_4, before4);

                    // j+1 ****************************
                    after1 = _mm256_fmadd_pd(m1_i_j1_x4, mat2rowj1_1, after1);
                    after2 = _mm256_fmadd_pd(m1_i_j1_x4, mat2rowj1_2, after2);
                    after3 = _mm256_fmadd_pd(m1_i_j1_x4, mat2rowj1_3, after3);
                    after4 = _mm256_fmadd_pd(m1_i_j1_x4, mat2rowj1_4, after4);

                    // j+2 ****************************
                    after1 = _mm256_fmadd_pd(m1_i_j2_x4, mat2rowj2_1, after1);
                    after2 = _mm256_fmadd_pd(m1_i_j2_x4, mat2rowj2_2, after2);
                    after3 = _mm256_fmadd_pd(m1_i_j2_x4, mat2rowj2_3, after3);
                    after4 = _mm256_fmadd_pd(m1_i_j2_x4, mat2rowj2_4, after4);

                    // j+3 ****************************
                    after1 = _mm256_fmadd_pd(m1_i_j3_x4, mat2rowj3_1, after1);
                    after2 = _mm256_fmadd_pd(m1_i_j3_x4, mat2rowj3_2, after2);
                    after3 = _mm256_fmadd_pd(m1_i_j3_x4, mat2rowj3_3, after3);
                    after4 = _mm256_fmadd_pd(m1_i_j3_x4, mat2rowj3_4, after4);

                    // // j+4 ****************************
                    // after1 = _mm256_fmadd_pd(m1_i_j4_x4, mat2rowj4_1, after1);
                    // after2 = _mm256_fmadd_pd(m1_i_j4_x4, mat2rowj4_2, after2);
                    // after3 = _mm256_fmadd_pd(m1_i_j4_x4, mat2rowj4_3, after3);
                    // after4 = _mm256_fmadd_pd(m1_i_j4_x4, mat2rowj4_4, after4);

                    // // j+5 ****************************
                    // after1 = _mm256_fmadd_pd(m1_i_j5_x4, mat2rowj5_1, after1);
                    // after2 = _mm256_fmadd_pd(m1_i_j5_x4, mat2rowj5_2, after2);
                    // after3 = _mm256_fmadd_pd(m1_i_j5_x4, mat2rowj5_3, after3);
                    // after4 = _mm256_fmadd_pd(m1_i_j5_x4, mat2rowj5_4, after4);

                    // // j+6 ****************************
                    // after1 = _mm256_fmadd_pd(m1_i_j6_x4, mat2rowj6_1, after1);
                    // after2 = _mm256_fmadd_pd(m1_i_j6_x4, mat2rowj6_2, after2);
                    // after3 = _mm256_fmadd_pd(m1_i_j6_x4, mat2rowj6_3, after3);
                    // after4 = _mm256_fmadd_pd(m1_i_j6_x4, mat2rowj6_4, after4);

                    // // j+7 ****************************
                    // after1 = _mm256_fmadd_pd(m1_i_j7_x4, mat2rowj7_1, after1);
                    // after2 = _mm256_fmadd_pd(m1_i_j7_x4, mat2rowj7_2, after2);
                    // after3 = _mm256_fmadd_pd(m1_i_j7_x4, mat2rowj7_3, after3);
                    // after4 = _mm256_fmadd_pd(m1_i_j7_x4, mat2rowj7_4, after4);
                    
            
                    _mm256_storeu_pd(result->whole_matrix+ i * result->cols + k, after1);
                    _mm256_storeu_pd(result->whole_matrix+ i * result->cols + k + 4, after2);
                    _mm256_storeu_pd(result->whole_matrix+ i * result->cols + k + 8, after3);
                    _mm256_storeu_pd(result->whole_matrix+ i * result->cols + k + 12, after4);

                }

                // The tailing case
                for(int k = tail; k < mat2->cols; k++) {
                    // printf("K is %d\n", k);
                    double matrix1_i_j = mat1->whole_matrix[i*mat1->cols+j];
                    double matrix2_j_k = mat2->whole_matrix[j*mat2->cols+k];

                    // j+1 ************************
                    double matrix1_i_j1 = mat1->whole_matrix[i*mat1->cols+j+1];
                    double matrix2_j1_k = mat2->whole_matrix[(j+1)*mat2->cols+k];

                    // j+2 ************************
                    double matrix1_i_j2 = mat1->whole_matrix[i*mat1->cols+j+2];
                    double matrix2_j2_k = mat2->whole_matrix[(j+2)*mat2->cols+k];

                    // j+3 ************************
                    double matrix1_i_j3 = mat1->whole_matrix[i*mat1->cols+j+3];
                    double matrix2_j3_k = mat2->whole_matrix[(j+3)*mat2->cols+k];

                    // // j+4 ************************
                    // double matrix1_i_j4 = mat1->whole_matrix[i*mat1->cols+j+4];
                    // double matrix2_j4_k = mat2->whole_matrix[(j+4)*mat2->cols+k];

                    // // j+5 ************************
                    // double matrix1_i_j5 = mat1->whole_matrix[i*mat1->cols+j+5];
                    // double matrix2_j5_k = mat2->whole_matrix[(j+5)*mat2->cols+k];

                    // // j+6 ************************
                    // double matrix1_i_j6 = mat1->whole_matrix[i*mat1->cols+j+6];
                    // double matrix2_j6_k = mat2->whole_matrix[(j+6)*mat2->cols+k];

                    // // j+7 ************************
                    // double matrix1_i_j7 = mat1->whole_matrix[i*mat1->cols+j+7];
                    // double matrix2_j7_k = mat2->whole_matrix[(j+7)*mat2->cols+k];


                    // printf("Adding matrix1_i_j %f * matrix2_j_k %f , matrix1_i_j1 %f matrix2_j1_k %f to %d %d\n",
                    //                 matrix1_i_j, matrix2_j_k, matrix1_i_j1, matrix2_j1_k, i, k);
                    // For every interation j, we sweep through result(tmp)[i][0:mat2->cols]
                    // tmp[i][k] += m1[i][j](fixed) * m2[j][k]
                    result->whole_matrix[i*result->cols + k] += matrix1_i_j * matrix2_j_k;

                    // j+1 ************************
                    result->whole_matrix[i*result->cols + k] += matrix1_i_j1 * matrix2_j1_k;

                    // j+2 ************************
                    result->whole_matrix[i*result->cols + k] += matrix1_i_j2 * matrix2_j2_k;

                    // j+3 ************************
                    result->whole_matrix[i*result->cols + k] += matrix1_i_j3 * matrix2_j3_k;

                    // // j+4 ************************
                    // result->whole_matrix[i*result->cols + k] += matrix1_i_j4 * matrix2_j4_k;

                    // // j+5 ************************
                    // result->whole_matrix[i*result->cols + k] += matrix1_i_j5 * matrix2_j5_k;

                    // // j+6 ************************
                    // result->whole_matrix[i*result->cols + k] += matrix1_i_j6 * matrix2_j6_k;

                    // // j+7 ************************
                    // result->whole_matrix[i*result->cols + k] += matrix1_i_j7 * matrix2_j7_k;
                }
            }

            // Tailing case
            
            for(int j = mat1->cols/4*4; j < mat1->cols; j += 1) {
                
                // ************For every element of matrix1[i][j]**********
                // m1_i_j_x4 == mat1->data[i][j] x 4 == mat1->whole_matrix[i*mat->cols+j]
                __m256d m1_i_j_x4 = _mm256_set1_pd(mat1->whole_matrix[i*mat1->cols+j]);
                
                // We extract mat1[i][j], then times it with all the elements in row mat[j]
                // then store them in every element in the row result[i]
                // mat2->cols = tail
                for(int k = 0; k < tail; k += step) {
                    /* In the end of this loop
                     * result(tmp)[i][k:k+4] += m1[i][j]_x4 * m2[j][k:k+4]
                     * 
                     * So, result[i][j] will get updated as in the next iteration of j
                     * since k will always starts at 0 in every j
                     */
                    // mat2->whole_matrix+ j * result->cols + k == mat2->data[j][k: k+4]
                    mat2rowj_1 = _mm256_loadu_pd(mat2->whole_matrix+ j * result->cols + k);
                    mat2rowj_2 = _mm256_loadu_pd(mat2->whole_matrix+ j * result->cols + k + 4);
                    mat2rowj_3 = _mm256_loadu_pd(mat2->whole_matrix+ j * result->cols + k + 8);
                    mat2rowj_4 = _mm256_loadu_pd(mat2->whole_matrix+ j * result->cols + k + 12);

                    
                    // tmp->whole_matrix+ i * result->cols + k == tmp->data[i][k:k+4]
                    before1 = _mm256_loadu_pd(result->whole_matrix+ i * result->cols + k);
                    before2 = _mm256_loadu_pd(result->whole_matrix+ i * result->cols + k + 4);
                    before3 = _mm256_loadu_pd(result->whole_matrix+ i * result->cols + k + 8);
                    before4 = _mm256_loadu_pd(result->whole_matrix+ i * result->cols + k + 12);

                    after1 = _mm256_fmadd_pd(m1_i_j_x4, mat2rowj_1, before1);
                    after2 = _mm256_fmadd_pd(m1_i_j_x4, mat2rowj_2, before2);
                    after3 = _mm256_fmadd_pd(m1_i_j_x4, mat2rowj_3, before3);
                    after4 = _mm256_fmadd_pd(m1_i_j_x4, mat2rowj_4, before4);

                    _mm256_storeu_pd(result->whole_matrix+ i * result->cols + k, after1);
                    _mm256_storeu_pd(result->whole_matrix+ i * result->cols + k + 4, after2);
                    _mm256_storeu_pd(result->whole_matrix+ i * result->cols + k + 8, after3);
                    _mm256_storeu_pd(result->whole_matrix+ i * result->cols + k + 12, after4);
                }

                // The tailing case
                for(int k = tail; k < mat2->cols; k++) {
                    double matrix1_i_j = mat1->whole_matrix[i*mat1->cols+j];
                    double matrix2_j_k = mat2->whole_matrix[j*mat2->cols+k];
                    
                    // TODO
                    // For every interation j, we sweep through result(tmp)[i][0:mat2->cols]
                    // tmp[i][k] += m1[i][j](fixed) * m2[j][k]
                    result->whole_matrix[i*result->cols + k] += matrix1_i_j * matrix2_j_k;
                }
            }

        }

        // #pragma omp parallel for
        // for (int i = 0; i < (result->rows * result->cols); i+=1) {
        //     result->whole_matrix[i] = tmp->whole_matrix[i];
        // }
        // deallocate_matrix(tmp);
        return 0;
    }


    // Actually doing the hard work, matrix multiplication
    for(int r = 0; r < mat1->rows; r++) {
        for(int c = 0; c < mat2->cols; c++) {
            double curr_sum = 0;
            for(int i = 0; i < mat1->cols; i++) {
                curr_sum += mat1->data[r][i] * mat2->data[i][c]; 
            }
            result->data[r][c] = curr_sum;
        }
    }
    return 0;

}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {
    /* TODO: YOUR CODE HERE */
    if (pow < 0) {
        PyErr_SetString(PyExc_ValueError, "Error when pow: negative pow");
        return -1;
    }

    if (mat->rows != mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Error when pow: can't pow a non-square matrix");
        return -1;
    }

    matrix* temp_result = NULL;
    int success;
    success = allocate_matrix(&temp_result, result->rows, result->cols);
    if (success != 0) {
        return -1;
    }

    matrix* exponential = NULL;
    if (allocate_matrix(&exponential, mat->rows, mat->cols) != 0 ||
        add_matrix(exponential, mat, exponential) != 0) {
        return -1;
    }

    matrix * temp_exp = NULL;
    success = allocate_matrix(&temp_exp, exponential->rows, exponential->cols);
    if (success != 0){
        return -1;
    }

    // Result = Identity
    for(int i = 0; i < result->rows; i++) {
        result->whole_matrix[i*result->cols+i] = 1.0;
    }

    register __m256d zero = _mm256_set1_pd (0);
    
    while(pow >= 1) {
        
        if (pow & 1 == 1) {
            success = mul_matrix(temp_result, result, exponential);
            if (success != 0) {
                return success;
            }

            // for(int i = 0; i < result->rows * result->cols; i++) {
            //     result->whole_matrix[i] = temp_result->whole_matrix[i];
            //     temp_result->whole_matrix[i] = 0;
            // }
            #pragma omp parallel for
            for(int i = 0; i < result->rows * result->cols/32*32; i+=32) {
                _mm256_storeu_pd (result->whole_matrix+i, _mm256_loadu_pd(temp_result->whole_matrix + i));
                _mm256_storeu_pd (result->whole_matrix+i+4, _mm256_loadu_pd(temp_result->whole_matrix + i+4));
                _mm256_storeu_pd (result->whole_matrix+i+8, _mm256_loadu_pd(temp_result->whole_matrix + i+8));
                _mm256_storeu_pd (result->whole_matrix+i+12, _mm256_loadu_pd(temp_result->whole_matrix + i+12));
                _mm256_storeu_pd (result->whole_matrix+i+16, _mm256_loadu_pd(temp_result->whole_matrix + i+16));
                _mm256_storeu_pd (result->whole_matrix+i+20, _mm256_loadu_pd(temp_result->whole_matrix + i+20));
                _mm256_storeu_pd (result->whole_matrix+i+24, _mm256_loadu_pd(temp_result->whole_matrix + i+24));
                _mm256_storeu_pd (result->whole_matrix+i+28, _mm256_loadu_pd(temp_result->whole_matrix + i+28));


                _mm256_storeu_pd (temp_result->whole_matrix+i, zero);
                _mm256_storeu_pd (temp_result->whole_matrix+i+4, zero);
                _mm256_storeu_pd (temp_result->whole_matrix+i+8, zero);
                _mm256_storeu_pd (temp_result->whole_matrix+i+12, zero);
                _mm256_storeu_pd (temp_result->whole_matrix+i+16, zero);
                _mm256_storeu_pd (temp_result->whole_matrix+i+20, zero);
                _mm256_storeu_pd (temp_result->whole_matrix+i+24, zero);
                _mm256_storeu_pd (temp_result->whole_matrix+i+28, zero);
            }
            for(int i = result->cols * result->rows / 32 * 32; i < result->cols * result->rows; i++) {
                result->whole_matrix[i] = temp_result->whole_matrix[i];
                temp_result->whole_matrix[i] = 0;
            }

        }
        
        success = mul_matrix(temp_exp, exponential, exponential);
        if (success != 0) {
            return success;
        }
        pow /= 2;
        // // exponential = temp_exp
        // for(int i = 0; i < exponential->rows * exponential->cols; i++) {
        //     exponential->whole_matrix[i] = temp_exp->whole_matrix[i];
        //     temp_exp->whole_matrix[i] = 0;
        // }
        

        #pragma omp parallel for
        for(int i = 0; i < exponential->cols * exponential->rows / 32 * 32; i += 32) {
            _mm256_storeu_pd (exponential->whole_matrix+i, _mm256_loadu_pd(temp_exp->whole_matrix + i));
            _mm256_storeu_pd (exponential->whole_matrix+i+4, _mm256_loadu_pd(temp_exp->whole_matrix + i+4));
            _mm256_storeu_pd (exponential->whole_matrix+i+8, _mm256_loadu_pd(temp_exp->whole_matrix + i+8));
            _mm256_storeu_pd (exponential->whole_matrix+i+12, _mm256_loadu_pd(temp_exp->whole_matrix + i+12));
            _mm256_storeu_pd (exponential->whole_matrix+i+16, _mm256_loadu_pd(temp_exp->whole_matrix + i+16));
            _mm256_storeu_pd (exponential->whole_matrix+i+20, _mm256_loadu_pd(temp_exp->whole_matrix + i+20));
            _mm256_storeu_pd (exponential->whole_matrix+i+24, _mm256_loadu_pd(temp_exp->whole_matrix + i+24));
            _mm256_storeu_pd (exponential->whole_matrix+i+28, _mm256_loadu_pd(temp_exp->whole_matrix + i+28));


            _mm256_storeu_pd (temp_exp->whole_matrix+i, zero);
            _mm256_storeu_pd (temp_exp->whole_matrix+i+4, zero);
            _mm256_storeu_pd (temp_exp->whole_matrix+i+8, zero);
            _mm256_storeu_pd (temp_exp->whole_matrix+i+12, zero);
            _mm256_storeu_pd (temp_exp->whole_matrix+i+16, zero);
            _mm256_storeu_pd (temp_exp->whole_matrix+i+20, zero);
            _mm256_storeu_pd (temp_exp->whole_matrix+i+24, zero);
            _mm256_storeu_pd (temp_exp->whole_matrix+i+28, zero);
        }

        for(int i = exponential->cols * exponential->rows / 32 * 32; i < exponential->cols * exponential->rows; i++) {
            exponential->whole_matrix[i] = temp_exp->whole_matrix[i];
            temp_exp->whole_matrix[i] = 0;
        }
        
    }
    deallocate_matrix(temp_result);
    deallocate_matrix(temp_exp);
    return 0;







    ///********Old code****

    // // if pow = 0:
    // // case 1: result == NULL
    // // do allocate result and do not do fill matrix
    // // case 2: result != NULL
    // // do not do allocate result and do fill matrix
    // if (pow==0){
    //     // (1) allocate
        
    //     if (result == NULL){
    //         int success = allocate_matrix(&result, mat->rows, mat->cols);
    //         if (success != 0) {
    //             return success;
    //         }
    //         // (2) Set the result to be identity matrix
    //         double identity = 1.0;
    //         #pragma omp parallel for
    //         for(int i = 0; i < result->cols; i++) {
    //             // result[i][i] = identity
    //             result->whole_matrix[result->cols*i + i] = identity;
    //         }
    //         return 0;
    //     } else {
    //         // when result != null
    //         // (1) fill matrix
            
    //         fill_matrix(result, 0.0);
    //         double identity = 1.0;

    //         // (2) Set the result to be identity matrix
    //         #pragma omp parallel for
    //         for(int i = 0; i < result->cols; i++) {
    //             // result[i][i] = identity
    //             result->whole_matrix[result->cols*i + i] = identity;
    //         }
    //         return 0;
    //     }
    // }

    // // in the following code, we know that pow != 0

    // // Allocate result matrix if it's not allocated
    // // if (result == NULL) {
    // //     int success = allocate_matrix(&result, mat->rows, mat->cols);
    // //     if (success != 0) {
    // //         return success;
    // //     }
    // // } 

    // if (!(result->rows == mat->rows && result->cols == mat->cols)) {
    //     PyErr_SetString(PyExc_RuntimeError, "Error when pow: shape not matched");
    //     return -1;
    // }

    // // // Always make result to be identity 
    // // // O(N^2)
    // // // not necessary
    // // if (pow == 0) {
    // //     fill_matrix(result, 0.0);

    // //     double identity = 1.0;
    // //     // Set the result to be identity matrix
    // //     #pragma omp parallel for
    // //     for(int i = 0; i < result->cols; i++) {
    // //         // result[i][i] = identity
    // //         result->whole_matrix[result->cols*i + i] = identity;
    // //     }
    // //     return 0;
    // // }
    // int first_time = 0;
    // if (mat->parent == NULL) {
    //     matrix* exponential = NULL;
    //     if (allocate_matrix(&exponential, mat->rows, mat->cols) != 0 ||
    //         add_matrix(exponential, mat, exponential) != 0) {
    //         return -1;
    //     }
        
       
    //     matrix* tmp_expo = NULL;
    //     if (allocate_matrix(&tmp_expo, exponential->rows, exponential->cols) != 0) {
    //         return -1;
    //     }

    //     matrix* tmp_result = NULL;
    //     if (allocate_matrix(&tmp_result, result->rows, result->cols) != 0) {
    //         return -1;
    //     }

    //     register __m256d zero = _mm256_set1_pd (0);

    //     int first_time = 0;
    //     while(pow >= 1) {
    //         if (pow & 1 == 1) {
    //             // temp = exp*result
    //             // tmp = result * exponential
    //             // if (this is the 1st iteration)
    //             // tmp-result = expon
    //             if (first_time == 0) {
    //                 // Deep copy from exponential to result
    //                 // result = exponential

    //                 #pragma omp parallel for
    //                 for(int i = 0; i < result->cols * result->rows / 32 * 32; i+= 32) {
    //                     _mm256_storeu_pd (result->whole_matrix+i, _mm256_loadu_pd(exponential->whole_matrix + i));
    //                     _mm256_storeu_pd (result->whole_matrix+i+4, _mm256_loadu_pd(exponential->whole_matrix + i+4));
    //                     _mm256_storeu_pd (result->whole_matrix+i+8, _mm256_loadu_pd(exponential->whole_matrix + i+8));
    //                     _mm256_storeu_pd (result->whole_matrix+i+12, _mm256_loadu_pd(exponential->whole_matrix + i+12));
    //                     _mm256_storeu_pd (result->whole_matrix+i+16, _mm256_loadu_pd(exponential->whole_matrix + i+16));
    //                     _mm256_storeu_pd (result->whole_matrix+i+20, _mm256_loadu_pd(exponential->whole_matrix + i+20));
    //                     _mm256_storeu_pd (result->whole_matrix+i+24, _mm256_loadu_pd(exponential->whole_matrix + i+24));
    //                     _mm256_storeu_pd (result->whole_matrix+i+28, _mm256_loadu_pd(exponential->whole_matrix + i+28));
    //                 }
    //                 for(int i = result->cols * result->rows / 32 * 32; i < result->cols * result->rows; i++) {
    //                     result->whole_matrix[i] = exponential->whole_matrix[i];
    //                 }
    //                 first_time = 1;
    //             } 
                
    //             else{
    //                 int success = mul_matrix(tmp_result, result, exponential);
    //                 if (success != 0) {
    //                     return success;
    //                 }
                    
                    

    //                 // Deep copy from tmp_result to result
    //                 #pragma omp parallel for
    //                 for(int i = 0; i < result->cols * result->rows / 32 * 32; i+= 32) {
    //                     _mm256_storeu_pd (result->whole_matrix+i, _mm256_loadu_pd(tmp_result->whole_matrix + i));
    //                     _mm256_storeu_pd (result->whole_matrix+i+4, _mm256_loadu_pd(tmp_result->whole_matrix + i+4));
    //                     _mm256_storeu_pd (result->whole_matrix+i+8, _mm256_loadu_pd(tmp_result->whole_matrix + i+8));
    //                     _mm256_storeu_pd (result->whole_matrix+i+12, _mm256_loadu_pd(tmp_result->whole_matrix + i+12));
    //                     _mm256_storeu_pd (result->whole_matrix+i+16, _mm256_loadu_pd(tmp_result->whole_matrix + i+16));
    //                     _mm256_storeu_pd (result->whole_matrix+i+20, _mm256_loadu_pd(tmp_result->whole_matrix + i+20));
    //                     _mm256_storeu_pd (result->whole_matrix+i+24, _mm256_loadu_pd(tmp_result->whole_matrix + i+24));
    //                     _mm256_storeu_pd (result->whole_matrix+i+28, _mm256_loadu_pd(tmp_result->whole_matrix + i+28));


    //                     _mm256_storeu_pd (tmp_result->whole_matrix+i, zero);
    //                     _mm256_storeu_pd (tmp_result->whole_matrix+i+4, zero);
    //                     _mm256_storeu_pd (tmp_result->whole_matrix+i+8, zero);
    //                     _mm256_storeu_pd (tmp_result->whole_matrix+i+12, zero);
    //                     _mm256_storeu_pd (tmp_result->whole_matrix+i+16, zero);
    //                     _mm256_storeu_pd (tmp_result->whole_matrix+i+20, zero);
    //                     _mm256_storeu_pd (tmp_result->whole_matrix+i+24, zero);
    //                     _mm256_storeu_pd (tmp_result->whole_matrix+i+28, zero);
    //                 }

    //                 for(int i = result->cols * result->rows / 32 * 32; i < result->cols * result->rows; i++) {
    //                     result->whole_matrix[i] = tmp_result->whole_matrix[i];
    //                     tmp_result->whole_matrix[i] = 0;
    //                 }
                    
                    
    //             }
    //         }
            
    //         int success = mul_matrix(tmp_expo, exponential, exponential);
            
    //         // Deep copy from tmp to expoential
    //         #pragma omp parallel for
    //         for(int i = 0; i < exponential->cols * exponential->rows / 32 * 32; i += 32) {
    //             _mm256_storeu_pd (exponential->whole_matrix+i, _mm256_loadu_pd(tmp_expo->whole_matrix + i));
    //             _mm256_storeu_pd (exponential->whole_matrix+i+4, _mm256_loadu_pd(tmp_expo->whole_matrix + i+4));
    //             _mm256_storeu_pd (exponential->whole_matrix+i+8, _mm256_loadu_pd(tmp_expo->whole_matrix + i+8));
    //             _mm256_storeu_pd (exponential->whole_matrix+i+12, _mm256_loadu_pd(tmp_expo->whole_matrix + i+12));
    //             _mm256_storeu_pd (exponential->whole_matrix+i+16, _mm256_loadu_pd(tmp_expo->whole_matrix + i+16));
    //             _mm256_storeu_pd (exponential->whole_matrix+i+20, _mm256_loadu_pd(tmp_expo->whole_matrix + i+20));
    //             _mm256_storeu_pd (exponential->whole_matrix+i+24, _mm256_loadu_pd(tmp_expo->whole_matrix + i+24));
    //             _mm256_storeu_pd (exponential->whole_matrix+i+28, _mm256_loadu_pd(tmp_expo->whole_matrix + i+28));
                

    //             _mm256_storeu_pd (tmp_expo->whole_matrix+i, zero);
    //             _mm256_storeu_pd (tmp_expo->whole_matrix+i+4, zero);
    //             _mm256_storeu_pd (tmp_expo->whole_matrix+i+8, zero);
    //             _mm256_storeu_pd (tmp_expo->whole_matrix+i+12, zero);
    //             _mm256_storeu_pd (tmp_expo->whole_matrix+i+16, zero);
    //             _mm256_storeu_pd (tmp_expo->whole_matrix+i+20, zero);
    //             _mm256_storeu_pd (tmp_expo->whole_matrix+i+24, zero);
    //             _mm256_storeu_pd (tmp_expo->whole_matrix+i+28, zero);
    //         }

    //         for(int i = exponential->cols * exponential->rows / 32 * 32; i < exponential->cols * exponential->rows; i++) {
    //             exponential->whole_matrix[i] = tmp_expo->whole_matrix[i];
    //             tmp_expo->whole_matrix[i] = 0;
    //         }
    //         if (success != 0) {
    //             return success;
    //         }
    //         pow /= 2;
    //     }
    //     deallocate_matrix(tmp_result);
    //     deallocate_matrix(tmp_expo);
    //     return 0;
    // }
    
    // // for(int i = 0; i < pow; i++) {
    // //     mul_matrix(result, result, mat);
    // // }

    // return 0;
}

/*
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    // Allocate result matrix if it's not allocated
    if (result == NULL) {
        int success = allocate_matrix(&result, mat->rows, mat->cols);
        if (success != 0) {
            return success;
        }
    } else {
        // check if result matrix has the same shape as mat1&mat2
        if (!(result->rows == mat->rows && result->cols == mat->cols)) {
            PyErr_SetString(PyExc_RuntimeError, "Error when add: shape not matched");
            return -1;
        }
    }



    if (mat->parent == NULL) {
        int size = mat->rows * mat->cols;
        #pragma omp parallel for
        for(int i = 0; i < size; i++) {
            result->whole_matrix[i] = -mat->whole_matrix[i];
        }
        return 0;
    }
    

    #pragma omp parallel for
    for(int r = 0; r < mat->rows; r++) {
        for(int c = 0; c < mat->cols; c++) {
            result->data[r][c] = -mat->data[r][c];
        }
    }

    return 0;
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int abs_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */

    // Allocate result matrix if it's not allocated
    if (result == NULL) {
        int success = allocate_matrix(&result, mat->rows, mat->cols);
        if (success != 0) {
            return success;
        }
    } else {
        // check if result matrix has the same shape as mat1&mat2
        if (!(result->rows == mat->rows && result->cols == mat->cols)) {
            PyErr_SetString(PyExc_RuntimeError, "Error when add: shape not matched");
            return -1;
        }
    }

    // Only for root matrix for speed up
    if (mat->parent == NULL) {
        int size = mat->rows * mat->cols;
        #pragma omp parallel for
        for(int i = 0; i < size; i++) {
            result->whole_matrix[i] = mat->whole_matrix[i] < 0 ? -mat->whole_matrix[i] : mat->whole_matrix[i];
        }
        return 0;
    }

    #pragma omp parallel for
    for(int r = 0; r < mat->rows; r++) {
        for(int c = 0; c < mat->cols; c++) {
            result->data[r][c] = mat->data[r][c] < 0 ? -mat->data[r][c] : mat->data[r][c];
        }
    }
    return 0;

}