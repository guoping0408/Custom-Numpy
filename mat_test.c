#include <stdio.h>

#include "CUnit/Basic.h"
#include "CUnit/CUnit.h"
#include "matrix.h"

/* Test Suite setup and cleanup functions: */
int init_suite(void) { return 0; }
int clean_suite(void) { return 0; }

/************* Test case functions ****************/
void add_test(void) {
    matrix *result = NULL;
    matrix *mat1 = NULL;
    matrix *mat2 = NULL;
    CU_ASSERT_EQUAL(allocate_matrix(&result, 2, 2), 0);
    CU_ASSERT_EQUAL(allocate_matrix(&mat1, 2, 2), 0);
    CU_ASSERT_EQUAL(allocate_matrix(&mat2, 2, 2), 0);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            set(mat1, i, j, i * 2 + j);
            set(mat2, i, j, i * 2 + j);
        }
    }
    add_matrix(result, mat1, mat2);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            CU_ASSERT_EQUAL(get(result, i, j), 2 * (i * 2 + j));
        }
    }
    printf("Deallocating... in add_test\n");
    deallocate_matrix(result);
    deallocate_matrix(mat1);
    deallocate_matrix(mat2);
}

void sub_test(void) {
    matrix *result = NULL;
    matrix *mat1 = NULL;
    matrix *mat2 = NULL;
    CU_ASSERT_EQUAL(allocate_matrix(&result, 2, 2), 0);
    CU_ASSERT_EQUAL(allocate_matrix(&mat1, 2, 2), 0);
    CU_ASSERT_EQUAL(allocate_matrix(&mat2, 2, 2), 0);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            set(mat1, i, j, i * 2 + j);
            set(mat2, i, j, (i * 2 + j) * 3);
        }
    }
    sub_matrix(result, mat1, mat2);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            CU_ASSERT_EQUAL(get(result, i, j), (-2) * (i * 2 + j));
        }
    }
    printf("Deallocating... in sub_test\n");
    deallocate_matrix(result);
    deallocate_matrix(mat1);
    deallocate_matrix(mat2);
}

void mul_test(void) {
    matrix *result = NULL;
    matrix *mat1 = NULL;
    matrix *mat2 = NULL;
    CU_ASSERT_EQUAL(allocate_matrix(&result, 3, 3), 0);
    CU_ASSERT_EQUAL(allocate_matrix(&mat1, 3, 3), 0);
    CU_ASSERT_EQUAL(allocate_matrix(&mat2, 3, 3), 0);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            set(mat1, i, j, i * 3 + j + 1);
            set(mat2, i, j, i * 3 + j + 1);
        }
    }
    printf("The matrix1 is \n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f ", get(mat1, i, j));
        }
        printf("\n");
    }

    printf("The matrix2 is \n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f ", get(mat2, i, j));
        }
        printf("\n");
    }
    
    mul_matrix(result, mat1, mat2);

    printf("The result is \n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            printf("%f ", get(result, i, j));
        }
        printf("\n");
    }

    CU_ASSERT_EQUAL(get(result, 0, 0), 30);
    CU_ASSERT_EQUAL(get(result, 0, 1), 36);
    CU_ASSERT_EQUAL(get(result, 0, 2), 42);
    CU_ASSERT_EQUAL(get(result, 1, 0), 66);
    CU_ASSERT_EQUAL(get(result, 1, 1), 81);
    CU_ASSERT_EQUAL(get(result, 1, 2), 96);
    CU_ASSERT_EQUAL(get(result, 2, 0), 102);
    CU_ASSERT_EQUAL(get(result, 2, 1), 126);
    CU_ASSERT_EQUAL(get(result, 2, 2), 150);
    printf("Deallocating... in mul_test\n");
    deallocate_matrix(result);
    deallocate_matrix(mat1);
    deallocate_matrix(mat2);
}

void neg_test(void) {
    matrix *result = NULL;
    matrix *mat = NULL;
    CU_ASSERT_EQUAL(allocate_matrix(&result, 2, 2), 0);
    CU_ASSERT_EQUAL(allocate_matrix(&mat, 2, 2), 0);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            set(mat, i, j, i * 2 + j);
        }
    }
    neg_matrix(result, mat);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            CU_ASSERT_EQUAL(get(result, i, j), -(i * 2 + j));
        }
    }
    printf("Deallocating... in neg_test\n");
    deallocate_matrix(result);
    deallocate_matrix(mat);
}

void abs_test(void) {
    matrix *result = NULL;
    matrix *mat = NULL;
    CU_ASSERT_EQUAL(allocate_matrix(&result, 2, 2), 0);
    CU_ASSERT_EQUAL(allocate_matrix(&mat, 2, 2), 0);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            if (j % 2 == 0)
                set(mat, i, j, i * 2 + j);
            else
                set(mat, i, j, -(i * 2 + j));
        }
    }
    abs_matrix(result, mat);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            CU_ASSERT_EQUAL(get(result, i, j), i * 2 + j);
        }
    }
    printf("Deallocating... in abs_test\n");
    deallocate_matrix(result);
    deallocate_matrix(mat);
}

void pow_test(void) {
    matrix *result = NULL;
    matrix *mat = NULL;
    CU_ASSERT_EQUAL(allocate_matrix(&result, 2, 2

                                   ),
                    0);
    CU_ASSERT_EQUAL(allocate_matrix(&mat, 2, 2), 0);
    set(mat, 0, 0, 1);
    set(mat, 0, 1, 1);
    set(mat, 1, 0, 1);
    set(mat, 1, 1, 0);
    pow_matrix(result, mat, 0);
    CU_ASSERT_EQUAL(get(result, 0, 0), 1);
    CU_ASSERT_EQUAL(get(result, 0, 1), 0);
    CU_ASSERT_EQUAL(get(result, 1, 0), 0);
    CU_ASSERT_EQUAL(get(result, 1, 1), 1);

    pow_matrix(result, mat, 3);
    CU_ASSERT_EQUAL(get(result, 0, 0), 3);
    CU_ASSERT_EQUAL(get(result, 0, 1), 2);
    CU_ASSERT_EQUAL(get(result, 1, 0), 2);
    CU_ASSERT_EQUAL(get(result, 1, 1), 1);
    
    for(int i = 0 ; i < mat->rows; i++) {
        for(int j = 0; j < mat->cols; j++) {
            printf("%f ", get(mat, i, j));
        }
        printf("\n");
    }
    pow_matrix(result, mat, 10);
    CU_ASSERT_EQUAL(get(result, 0, 0), 89);
    CU_ASSERT_EQUAL(get(result, 0, 1), 55);
    CU_ASSERT_EQUAL(get(result, 1, 0), 55);
    CU_ASSERT_EQUAL(get(result, 1, 1), 34);
    printf("get(result, 0, 0) %f \n", get(result, 0, 0));
    printf("get(result, 0, 1) %f \n", get(result, 0, 1));
    printf("get(result, 1, 0) %f \n", get(result, 1, 0));
    printf("get(result, 1, 1) %f \n", get(result, 1, 1));

    printf("Deallocating... in pow_test\n");
    deallocate_matrix(result);
    deallocate_matrix(mat);
}

void alloc_fail_test(void) {
    matrix *mat = NULL;
    CU_ASSERT_EQUAL(allocate_matrix(&mat, 0, 0), -1);
    CU_ASSERT_EQUAL(allocate_matrix(&mat, 0, 1), -1);
    CU_ASSERT_EQUAL(allocate_matrix(&mat, 1, 0), -1);
}

void alloc_success_test(void) {
    matrix *mat = NULL;
    CU_ASSERT_EQUAL(allocate_matrix(&mat, 3, 2), 0);
    CU_ASSERT_EQUAL(mat->parent, NULL);
    CU_ASSERT_EQUAL(mat->ref_cnt, 1);
    CU_ASSERT_EQUAL(mat->rows, 3);
    CU_ASSERT_EQUAL(mat->cols, 2);
    CU_ASSERT_NOT_EQUAL(mat->data, NULL);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            CU_ASSERT_EQUAL(get(mat, i, j), 0);
        }
    }
    printf("Deallocating... in alloc_success_test\n");
    deallocate_matrix(mat);
}

void alloc_ref_test(void) {
    matrix *mat1 = NULL;
    matrix *mat2 = NULL;
    matrix *from = NULL;
    allocate_matrix(&from, 3, 2);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            set(from, i, j, i * 2 + j);
        }
    }

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            printf("%f ", get(from, i, j));
        }
        printf("\n");
    }
    
    /* 2D slice */
    CU_ASSERT_EQUAL(allocate_matrix_ref(&mat1, from, 1, 0, 2, 2), 0);
    CU_ASSERT_PTR_EQUAL(mat1->parent, from);
    CU_ASSERT_EQUAL(mat1->parent->ref_cnt, 2);
    CU_ASSERT_EQUAL(mat1->rows, 2);
    CU_ASSERT_EQUAL(mat1->cols, 2);
    printf("2D slice\n");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            CU_ASSERT_EQUAL(get(mat1, i, j), get(from, i + 1, j));
            printf("get(mat1, i, j) %f    ", get(mat1, i, j));
            printf("get(from, i + 1, j) %f\n", get(from, i + 1, j));
        }
    }
    /* 1D slice */
    CU_ASSERT_EQUAL(allocate_matrix_ref(&mat2, from, 1, 0, 2, 1), 0);
    CU_ASSERT_PTR_EQUAL(mat2->parent, from);
    CU_ASSERT_EQUAL(mat2->parent->ref_cnt, 3);
    CU_ASSERT_EQUAL(mat2->rows, 2);
    CU_ASSERT_EQUAL(mat2->cols, 1);
    CU_ASSERT_NOT_EQUAL(mat2->is_1d, 0);
    printf("1D slice\n");
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 1; j++) {
            CU_ASSERT_EQUAL(get(mat2, i, j), get(from, i + 1, j));
            printf("get(mat2, i, j) %f    ", get(mat2, i, j));
            printf("get(from, i + 1, j) %f\n", get(from, i + 1, j));
        }
    }
    /* Now we compare the data in the reference matrix */
    printf("Deallocating... in alloc_ref_test\n");
    deallocate_matrix(from);
    deallocate_matrix(mat1);
    deallocate_matrix(mat2);
    printf("Finished deallocating... in alloc_ref_test\n");
}
void alloc_ref_test_recursion(void) {
    matrix *mat1 = NULL;
    matrix *mat2 = NULL;
    matrix *from = NULL;
    allocate_matrix(&from, 5, 5);
    for (int i = 0; i < 5; i++) {
        for (int j = 0; j < 5; j++) {
            set(from, i, j, i * 2 + j);
        }
    }

    CU_ASSERT_EQUAL(allocate_matrix_ref(&mat1, from, 1, 0, 3, 3), 0);
    CU_ASSERT_EQUAL(allocate_matrix_ref(&mat2, mat1, 0, 0, 2, 2), 0);

    printf("Deallocating... in alloc_ref_test_recursion\n");
    deallocate_matrix(from);
    deallocate_matrix(mat1);
    printf("We should only deallocate here\n");
    deallocate_matrix(mat2);
    printf("Finished deallocating... in alloc_ref_test_recursion\n");
}


/* Test the null case doesn't crash */
void dealloc_null_test(void) {
    matrix *mat = NULL;
    deallocate_matrix(mat);
}

void get_test(void) {
    matrix *mat = NULL;
    allocate_matrix(&mat, 2, 2);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            set(mat, i, j, i * 2 + j);
        }
    }
    CU_ASSERT_EQUAL(get(mat, 0, 0), 0);
    CU_ASSERT_EQUAL(get(mat, 0, 1), 1);
    CU_ASSERT_EQUAL(get(mat, 1, 0), 2);
    CU_ASSERT_EQUAL(get(mat, 1, 1), 3);
    printf("Deallocating... in get_test\n");
    deallocate_matrix(mat);
}

void set_test(void) {
    matrix *mat = NULL;
    allocate_matrix(&mat, 2, 2);
    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            set(mat, i, j, i * 2 + j);
        }
    }
    CU_ASSERT_EQUAL(get(mat, 0, 0), 0);
    CU_ASSERT_EQUAL(get(mat, 0, 1), 1);
    CU_ASSERT_EQUAL(get(mat, 1, 0), 2);
    CU_ASSERT_EQUAL(get(mat, 1, 1), 3);
    deallocate_matrix(mat);
}

/************* Test Runner Code goes here **************/

int main(void) {
    Py_Initialize();  // Need to call this so that Python.h functions won't
    // segfault
    CU_pSuite pSuite = NULL;

    /* initialize the CUnit test registry */
    if (CU_initialize_registry() != CUE_SUCCESS) return CU_get_error();

    /* add a suite to the registry */
    pSuite = CU_add_suite("mat_test_suite", init_suite, clean_suite);
    if (pSuite == NULL) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    /* add the tests to the suite */
    if ((CU_add_test(pSuite, "add_test", add_test) == NULL) ||
            (CU_add_test(pSuite, "sub_test", sub_test) == NULL) ||
            (CU_add_test(pSuite, "mul_test", mul_test) == NULL) ||
            (CU_add_test(pSuite, "neg_test", neg_test) == NULL) ||
            (CU_add_test(pSuite, "abs_test", abs_test) == NULL) ||
            (CU_add_test(pSuite, "pow_test", pow_test) == NULL) ||
            (CU_add_test(pSuite, "alloc_fail_test", alloc_fail_test) == NULL) ||
            (CU_add_test(pSuite, "alloc_success_test", alloc_success_test) == NULL) ||
            (CU_add_test(pSuite, "alloc_ref_test", alloc_ref_test) == NULL) ||
            (CU_add_test(pSuite, "dealloc_null_test", dealloc_null_test) == NULL) ||
            (CU_add_test(pSuite, "get_test", get_test) == NULL) ||
            (CU_add_test(pSuite, "set_test", set_test) == NULL) || 
            (CU_add_test(pSuite, "alloc_ref_test", set_test) == NULL)) {
        CU_cleanup_registry();
        return CU_get_error();
    }

    // Run all tests using the basic interface
    CU_basic_set_mode(CU_BRM_NORMAL);
    // CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    printf("\n");
    CU_basic_show_failures(CU_get_failure_list());
    printf("\n\n");

    /* Clean up registry and return */
    CU_cleanup_registry();
    return CU_get_error();
}