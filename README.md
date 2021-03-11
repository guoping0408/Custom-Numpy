# Custom-Numpy
This project is both a C project as well as a performance project—building a slower version of numpy called numc. numc is much faster than the naive implementations of matrix operations (but slower than the actual numpy, still). The Python-C interface is overloaded with some operators and defined some instance methods for numc.Matrix objects. Finally, the speed up the naive solution makes numc.Matrix operations faster.

Here's what I did in this project:
-
In matrix.h:
We declare a field called double* whole_matrix and we use it to do store our underlying 1-D matrix. We use the original double** data as a way to implement specifically for indexing of both slicing and non-slicing matrices. However, when it's a non-slicing matrix, we just directly use a whole_matrix.

In addition, we also create an int *ref_cnt_ptr. And every slicing matrix / parent matrix points to the same int *ref_cnt_ptr, so when we are calling deallcate_matrix, we are decrementing the value pointed by the same ref_cnt_ptr and actually free the underlying data when *ref_cnt_ptr == 0.


In matrix.c:

(1) allocate_matrix:
We calloc a big chunk of memory of size rows*cols for the whole_matrix, and we malloc memory of size ows and assign it to double** data. The root matrix has slice index 0, 0

(2) allocate_matrix_ref:
Each element in double** data points to the underlying whole_matrix based on the slicing index 

(3) deallocate_matrix:
We always use *(mat->ref_cnt_ptr) -= 1. And, we free the whole_matrix only when *(mat->ref_cnt_ptr) = 0;

(4) get / set:
We use naive 2d double** data to access the data when it's slicing matrix
We directly access the underlying whole_matrix 1-D when it's non-slicing matrix (root matrix)

(5) fill_matrix:
We use SIMD and openMP to parallelize storing data using the 1D whole_matrix when the matrix itself is a root matrix. We use the original naive 2D whole_matrix to fill the matrix when it's a slicing matrix

(6) add_matrix:
We use the original naive 2D whole_matrix to add two matrices together when it's a slicing matrix. When it's a both operand matrices are root matrices (parent == NULL)


(7) add matrix:
In this function we supported the addition between two matrices, which is simply adding all entries in the matrix together. As a result, we have to make sure that the dimensions of the two matrices have to equate, else we throw an error. We only consider speedup when the matrices are root matrices instead of slicing, that is why we do (if parent == NULL) before our speedup portion. For speedup, we employed two techniques. The first one SIMD, which is to store four doubles inside a register. We notice that we can use the 256 bit register, and one double only takes place of 64 double, which means that a 256 bit register can store up to 4 double. We utilized these characteristics and inside each loop, we added 4 doubles together for each 256 bit register. Aside from this, we also employed unrolling to minimize the number of times that we have to get inside and out of a loop. We choose our unrolling = 4, so in our while loop, each time we increase by steps and each step is equal to #unrolling * #SIMDS = 4 * 4 = 16. For matrices that are children matrix (slicing result), we do not employ any speedup techniques. We just access each of their entries by using their data attribute (2d array), do naive addition, and store the sum to the result matrices. 

(8) sub matrix:
In this function we supported the subtraction between two matrices, which is simply subtracting the value of all entries in the second matrix from the first matrix. The following is very similar to the add matrix. We have to make sure that the dimensions of the two matrices have to equate, else we throw an error. We only consider speedup when the matrices are root matrices instead of slicing, that is why we do (if parent == NULL) before our speedup portion. For speedup, we employed two techniques. The first one SIMD, which is to store four doubles inside a register. We notice that we can use the 256 bit register, and one double only takes place of 64 double, which means that a 256 bit register can store up to 4 double. We utilized this characteristics and inside each loop, we add 4 doubles together for each 256 bit register. Aside from this, we also employed unrolling to minimize the number of times that we have to get inside and out of a loop. We choose our unrolling = 4, so in our while loop, each time we increase by steps and each step is equal to #unrolling * #SIMDS = 4 * 4 = 16. For matrices that are children matrix (slicing result), we do not employ any speedup techniques. We just access each of their entries by using their data attribute (2d array), do naive subtraction, and store the subtraction result to the result matrices. 

(9) mul matrix:
In mul_matrix, we only speed up the the operation when both of the matrices are root matrices.
Specifically, in the outer loop, we ITERATE through the MAT1 and we multiply its element with  the corresponding ROW of the mat2 at once.

Moreover, we apply SIMD when multiplying one single element (in fact we take one element from mat1 and have 4 copy of it in one single 256-bit register) with the corresponding rows of the mat2.

In addition, we use unrolling of 4 so we get 4 elements at once while multiplying the mat1 elements with the entire rows of the mat2. Therefore after each second inner loop, we have 
result[i][:] += mat1[i][j(fixed)] * mat2[j(fixed)][:], so every elements int row i of result will get accumulated eventually  

By doing this, we can achieve great performance when accessing elements in the inner most loop, as we have 1D sequential access on the mat2 (as we are sweeping through the rows of the mat2).

In the end, since we have unrolling of 4 and SIMD step size of 4, we will do essentially 16 instruction in a time. But we also need to do one step at the time when the size of row of mat2 is not divisible by 16.

pow matrix:
We treat the power number as a binary number. For example, we consider pow = 5 as pow = 101. We want the result to be multiplied only when the binary bit is equal to one. That is why we do (if pow & 1 == 1) we do mul (result, result, exp), and we do mul (exp, exp, exp) all the time. 
We have several things else to take care of though. For example, our mul does not support matrix multiplication when the result pointer = one of the operand. In order to achieve result pointer = one of the operands, we must allocate a new result matrix inside mul, which slows us down. That is why we choose to not do allocate in mul and instead do allocate in pow. In pow function, we just have to allocate two temporary matrix outside of the while loop ONCE (one is for storing temporary result and the other is used for storing temporary exponential, which would then be deeply copied to result and exponential), which is significant less than doing allocation everytime when going inside mul. Another edge case we paid attention to is when pow = 0. In this case, we just want to return the identity matrix. We supported this edge case by assigning the value of result to the identity first before going inside the while loop. This way, when pow = 0, we would not go inside the while loop (while pow >= 1) at all, and the result would be what we have initialized it to, which would be the identity matrix. 

neg matrix:


In numc.c
Task 3 is all written in this doc.
(1) Number methods:
Those methods do arithmetic operation to 61c matrices. There is no speedup presented here as we just call the underlying arithmetic operation function in matrix.c. We, however, did take care of the error checking as well as allocating & freeing here.
(2) Instance methods:
    - Matrix61c_set_value: we appropriate extract the index row, col and value from the PyObject* args and apply appropriate error checking function before actually using them for the set fuction.

    - Matrix61c_get_value: we appropriate extract the index row, col from the PyObject* args and apply appropriate error checking function before actually using them for the get fuction.
    


(3) Indexing:
Subscript: we slice a matrix by manipulating the data attribute (2d array) in the 61c matrix. The slicing result and the original matrix share the same whole_matrix (1d array), but our function gets the slicing result by changing the pointers of data, the 2d array, so the data points to different locations of the whole_matrix array now.
Set_subscript: we first get the slicing result of the original matrix by calling subscript. Then, we divide the slicing result into 3 cases. 
CASE 1: resulting = 1 by 1
val must = int/float, else type error
CASE 2: resulting = 1d 
val must = 1d list with same length and all long inside the list (else value error), else type error
CASE 3: resulting = 2d
val must = 2d list with same length and all long inside the list (else value error), else type error

In setup.py
We just included all the required ARGS and DIRS to the extension and call setup
