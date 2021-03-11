#include "numc.h"
#include <structmember.h>

PyTypeObject Matrix61cType;

/* Helper functions for initalization of matrices and vectors */

/*
 * Return a tuple given rows and cols
 */
PyObject *get_shape(int rows, int cols) {
  if (rows == 1 || cols == 1) {
    return PyTuple_Pack(1, PyLong_FromLong(rows * cols));
  } else {
    return PyTuple_Pack(2, PyLong_FromLong(rows), PyLong_FromLong(cols));
  }
}
/*
 * Matrix(rows, cols, low, high). Fill a matrix random double values
 */
int init_rand(PyObject *self, int rows, int cols, unsigned int seed, double low,
              double high) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    rand_matrix(new_mat, seed, low, high);
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(rows, cols, val). Fill a matrix of dimension rows * cols with val
 */
int init_fill(PyObject *self, int rows, int cols, double val) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed)
        return alloc_failed;
    else {
        fill_matrix(new_mat, val);
        ((Matrix61c *)self)->mat = new_mat;
        ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    }
    return 0;
}

/*
 * Matrix(rows, cols, 1d_list). Fill a matrix with dimension rows * cols with 1d_list values
 */
int init_1d(PyObject *self, int rows, int cols, PyObject *lst) {
    if (rows * cols != PyList_Size(lst)) {
        PyErr_SetString(PyExc_ValueError, "Incorrect number of elements in list");
        return -1;
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    int count = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j, PyFloat_AsDouble(PyList_GetItem(lst, count)));
            count++;
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(2d_list). Fill a matrix with dimension len(2d_list) * len(2d_list[0])
 */
int init_2d(PyObject *self, PyObject *lst) {
    int rows = PyList_Size(lst);
    if (rows == 0) {
        PyErr_SetString(PyExc_ValueError,
                        "Cannot initialize numc.Matrix with an empty list");
        return -1;
    }
    int cols;
    if (!PyList_Check(PyList_GetItem(lst, 0))) {
        PyErr_SetString(PyExc_ValueError, "List values not valid");
        return -1;
    } else {
        cols = PyList_Size(PyList_GetItem(lst, 0));
    }
    for (int i = 0; i < rows; i++) {
        if (!PyList_Check(PyList_GetItem(lst, i)) ||
                PyList_Size(PyList_GetItem(lst, i)) != cols) {
            PyErr_SetString(PyExc_ValueError, "List values not valid");
            return -1;
        }
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j,
                PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(lst, i), j)));
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * This deallocation function is called when reference count is 0
 */
void Matrix61c_dealloc(Matrix61c *self) {
    deallocate_matrix(self->mat);
    Py_TYPE(self)->tp_free(self);
}

/* For immutable types all initializations should take place in tp_new */
PyObject *Matrix61c_new(PyTypeObject *type, PyObject *args,
                        PyObject *kwds) {
    /* size of allocated memory is tp_basicsize + nitems*tp_itemsize*/
    Matrix61c *self = (Matrix61c *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

/*
 * This matrix61c type is mutable, so needs init function. Return 0 on success otherwise -1
 */
int Matrix61c_init(PyObject *self, PyObject *args, PyObject *kwds) {
    /* Generate random matrices */
    if (kwds != NULL) {
        PyObject *rand = PyDict_GetItemString(kwds, "rand");
        if (!rand) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (!PyBool_Check(rand)) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (rand != Py_True) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        PyObject *low = PyDict_GetItemString(kwds, "low");
        PyObject *high = PyDict_GetItemString(kwds, "high");
        PyObject *seed = PyDict_GetItemString(kwds, "seed");
        double double_low = 0;
        double double_high = 1;
        unsigned int unsigned_seed = 0;

        if (low) {
            if (PyFloat_Check(low)) {
                double_low = PyFloat_AsDouble(low);
            } else if (PyLong_Check(low)) {
                double_low = PyLong_AsLong(low);
            }
        }

        if (high) {
            if (PyFloat_Check(high)) {
                double_high = PyFloat_AsDouble(high);
            } else if (PyLong_Check(high)) {
                double_high = PyLong_AsLong(high);
            }
        }

        if (double_low >= double_high) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        // Set seed if argument exists
        if (seed) {
            if (PyLong_Check(seed)) {
                unsigned_seed = PyLong_AsUnsignedLong(seed);
            }
        }

        PyObject *rows = NULL;
        PyObject *cols = NULL;
        if (PyArg_UnpackTuple(args, "args", 2, 2, &rows, &cols)) {
            if (rows && cols && PyLong_Check(rows) && PyLong_Check(cols)) {
                return init_rand(self, PyLong_AsLong(rows), PyLong_AsLong(cols), unsigned_seed, double_low,
                                 double_high);
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    }
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyObject *arg3 = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 3, &arg1, &arg2, &arg3)) {
        /* arguments are (rows, cols, val) */
        if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && (PyLong_Check(arg3)
                || PyFloat_Check(arg3))) {
            if (PyLong_Check(arg3)) {
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyLong_AsLong(arg3));
            } else
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyFloat_AsDouble(arg3));
        } else if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && PyList_Check(arg3)) {
            /* Matrix(rows, cols, 1D list) */
            return init_1d(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), arg3);
        } else if (arg1 && PyList_Check(arg1) && arg2 == NULL && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_2d(self, arg1);
        } else if (arg1 && arg2 && PyLong_Check(arg1) && PyLong_Check(arg2) && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), 0);
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return -1;
    }
}

/*
 * List of lists representations for matrices
 */
PyObject *Matrix61c_to_list(Matrix61c *self) {
    int rows = self->mat->rows;
    int cols = self->mat->cols;
    PyObject *py_lst = NULL;
    if (self->mat->is_1d) {  // If 1D matrix, print as a single list
        py_lst = PyList_New(rows * cols);
        int count = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(py_lst, count, PyFloat_FromDouble(get(self->mat, i, j)));
                count++;
            }
        }
    } else {  // if 2D, print as nested list
        py_lst = PyList_New(rows);
        for (int i = 0; i < rows; i++) {
            PyList_SetItem(py_lst, i, PyList_New(cols));
            PyObject *curr_row = PyList_GetItem(py_lst, i);
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(curr_row, j, PyFloat_FromDouble(get(self->mat, i, j)));
            }
        }
    }
    return py_lst;
}

PyObject *Matrix61c_class_to_list(Matrix61c *self, PyObject *args) {
    PyObject *mat = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 1, &mat)) {
        if (!PyObject_TypeCheck(mat, &Matrix61cType)) {
            PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
            return NULL;
        }
        Matrix61c* mat61c = (Matrix61c*)mat;
        return Matrix61c_to_list(mat61c);
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }
}

/*
 * Add class methods
 */
PyMethodDef Matrix61c_class_methods[] = {
    {"to_list", (PyCFunction)Matrix61c_class_to_list, METH_VARARGS, "Returns a list representation of numc.Matrix"},
    {NULL, NULL, 0, NULL}
};

/*
 * Matrix61c string representation. For printing purposes.
 */
PyObject *Matrix61c_repr(PyObject *self) {
    PyObject *py_lst = Matrix61c_to_list((Matrix61c *)self);
    return PyObject_Repr(py_lst);
}

/* NUMBER METHODS */

/*
 * Add the second numc.Matrix (Matrix61c) object to the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_add(Matrix61c* self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    if (self == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "self is null");
        return NULL;
    }
    
    if (!PyObject_TypeCheck(args, &Matrix61cType)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
    
    // Have confirmed that args is indeed a matrix. Cast its type.
    Matrix61c* b = (Matrix61c*) args;
    
    // Check dimension
    if (b->mat->rows != self->mat->rows){
        PyErr_SetString(PyExc_ValueError,
                        "add operation between two matrices with different dimensions");
        return NULL;
    }
    
    if (b->mat->cols != self->mat->cols){
        PyErr_SetString(PyExc_ValueError,
                        "add operation between two matrices with different dimensions");
        return NULL;
    }
    
    // create a new matrix for storing the result
    matrix* c;
    int alloc_failed = allocate_matrix(&c, self->mat->rows, self->mat->cols);
    
    // check allocation failure
    if (alloc_failed != 0) {
        PyErr_SetString(PyExc_RuntimeError, "allocation failed");
        return NULL;
    }
    
    int add_failed = add_matrix(c, self->mat, b->mat);
    // chekc add failure
    if (add_failed != 0) {
        PyErr_SetString(PyExc_RuntimeError, "add failed.");
        return NULL;
    }
    
    Matrix61c *final = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    final->mat = c;
    final->shape = get_shape(c->rows, c->cols);
    return (PyObject*)final;
}

/*
 * Substract the second numc.Matrix (Matrix61c) object from the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_sub(Matrix61c* self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
        if (self == NULL) {
            PyErr_SetString(PyExc_RuntimeError, "self is null");
            return NULL;
        }
        
        if (!PyObject_TypeCheck(args, &Matrix61cType)) {
            PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
            return NULL;
        }
        
        // Have confirmed that args is indeed a matrix. Cast its type.
        Matrix61c* b = (Matrix61c*) args;
        
        // Check dimension
        if (b->mat->rows != self->mat->rows){
            PyErr_SetString(PyExc_ValueError,
                            "sub operation between two matrices with different dimensions");
            return NULL;
        }
        
        if (b->mat->cols != self->mat->cols){
            PyErr_SetString(PyExc_ValueError,
                            "sub operation between two matrices with different dimensions");
            return NULL;
        }
        
        // create a new matrix for storing the result
        matrix* c;
        int alloc_failed = allocate_matrix(&c, self->mat->rows, self->mat->cols);
        
        // check allocation failure
        if (alloc_failed != 0) {
            PyErr_SetString(PyExc_RuntimeError, "allocation failed");
            return NULL;
        }
        
        int sub_failed = sub_matrix(c, self->mat, b->mat);
        // chekc add failure
        if (sub_failed != 0) {
            PyErr_SetString(PyExc_RuntimeError, "sub failed.");
            return NULL;
        }
        
        Matrix61c *final = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
        final->mat = c;
        final->shape = get_shape(c->rows, c->cols);
        return (PyObject*)final;
    }


/*
 * NOT element-wise multiplication. The first operand is self, and the second operand
 * can be obtained by casting `args`.
 */
PyObject *Matrix61c_multiply(Matrix61c* self, PyObject *args) {
    /* TODO: YOUR CODE HERE */
        if (self == NULL) {
            PyErr_SetString(PyExc_RuntimeError, "self is null");
            return NULL;
        }
        
        if (!PyObject_TypeCheck(args, &Matrix61cType)) {
            PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
            return NULL;
        }
        
        // Have confirmed that args is indeed a matrix. Cast its type.
        Matrix61c* b = (Matrix61c*) args;
        
        // Check dimension        
        if (self->mat->cols != b->mat->rows){
            PyErr_SetString(PyExc_ValueError,
                            "mul: b's rows != a's cols");
            return NULL;
        }
        
        // create a new matrix for storing the result
        matrix* c;
        int alloc_failed = allocate_matrix(&c, self->mat->rows, b->mat->cols);
        
        // check allocation failure
        if (alloc_failed != 0) {
            PyErr_SetString(PyExc_RuntimeError, "allocation failed");
            return NULL;
        }
        
        int mul_failed = mul_matrix(c, self->mat, b->mat);
        // chekc add failure
        if (mul_failed != 0) {
            PyErr_SetString(PyExc_RuntimeError, "mul failed.");
            return NULL;
        }
        
        Matrix61c *final = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
        final->mat = c;
        final->shape = get_shape(c->rows, c->cols);
        return (PyObject*)final;
    }


/*
 * Negates the given numc.Matrix.
 */
PyObject *Matrix61c_neg(Matrix61c* self) {
    /* TODO: YOUR CODE HERE */
        if (self == NULL) {
            PyErr_SetString(PyExc_RuntimeError, "self is null");
            return NULL;
        }
        
        // create a new matrix for storing the result
        matrix* c;
        int alloc_failed = allocate_matrix(&c, self->mat->rows, self->mat->cols);
        
        // check allocation failure
        if (alloc_failed != 0) {
            PyErr_SetString(PyExc_RuntimeError, "allocation failed");
            return NULL;
        }
        
        int neg_failed = neg_matrix(c, self->mat);

        // check add failure
        if (neg_failed != 0) {
            PyErr_SetString(PyExc_RuntimeError, "neg failed.");
            return NULL;
        }
        
        Matrix61c *final = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
        final->mat = c;
        final->shape = get_shape(c->rows, c->cols);
        return (PyObject*)final;
}

/*
 * Take the element-wise absolute value of this numc.Matrix.
 */
PyObject *Matrix61c_abs(Matrix61c *self) {
    /* TODO: YOUR CODE HERE */
    if (self == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "self is null");
        return NULL;
    }
    
    // create a new matrix for storing the result
    matrix* c;
    int alloc_failed = allocate_matrix(&c, self->mat->rows, self->mat->cols);
    
    // check allocation failure
    if (alloc_failed != 0) {
        PyErr_SetString(PyExc_RuntimeError, "allocation failed");
        return NULL;
    }
    
    int abs_failed = abs_matrix(c, self->mat);
    // check add failure
    if (abs_failed != 0) {
        PyErr_SetString(PyExc_RuntimeError, "abs failed.");
        return NULL;
    }
    
    Matrix61c *final = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    final->mat = c;
    final->shape = get_shape(c->rows, c->cols);
    return (PyObject*)final;
    }

/*
 * Raise numc.Matrix (Matrix61c) to the `pow`th power. You can ignore the argument `optional`.
 */
PyObject *Matrix61c_pow(Matrix61c *self, PyObject *pow, PyObject *optional) {
    /* TODO: YOUR CODE HERE */
    if (self == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "self is null");
        return NULL;
    }

    // Throws error when pow is not an integer.
    if (!PyObject_TypeCheck(pow, &PyLong_Type)) {
        PyErr_SetString(PyExc_TypeError, "Pow must of type PyLong_Type!");
        return NULL;
    }
    int power = PyLong_AsLong(pow);

    // ValueError if self is not a square matrix or if pow is negative.
    if (power < 0 || self->mat->rows != self->mat->cols) {
        PyErr_SetString(PyExc_ValueError, " Not a square matrix or if pow is negative.");
        return NULL;
    }
    
    // create a new matrix for storing the result
    matrix* c;
    int alloc_failed = allocate_matrix(&c, self->mat->rows, self->mat->cols);
    
    // check allocation failure
    if (alloc_failed) {
        PyErr_SetString(PyExc_RuntimeError, "allocation failed");
        return NULL;
    }
    
    int pow_failed = pow_matrix(c, self->mat, power);
    // chekc add failure
    if (pow_failed) {
        PyErr_SetString(PyExc_RuntimeError, "pow failed.");
        return NULL;
    }
    
    Matrix61c *final = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
    final->mat = c;
    final->shape = get_shape(c->rows, c->cols);
    return (PyObject*)final;
}

/*
 * Create a PyNumberMethods struct for overloading operators with all the number methods you have
 * define. You might find this link helpful: https://docs.python.org/3.6/c-api/typeobj.html
 */
PyNumberMethods Matrix61c_as_number = {
    /* TODO: YOUR CODE HERE */
    .nb_add = (binaryfunc) Matrix61c_add,
    .nb_subtract = (binaryfunc) Matrix61c_sub,
    .nb_multiply = (binaryfunc) Matrix61c_multiply,
    .nb_remainder = (binaryfunc) Py_NotImplemented,
    .nb_divmod = (binaryfunc) Py_NotImplemented,
    .nb_power = (ternaryfunc) Matrix61c_pow,
    .nb_negative = (unaryfunc) Matrix61c_neg,
    .nb_positive = (unaryfunc) Py_NotImplemented,
    .nb_absolute = (unaryfunc) Matrix61c_abs
};


/* INSTANCE METHODS */

/*
 * Given a numc.Matrix self, parse `args` to (int) row, (int) col, and (double/int) val.
 * Return None in Python (this is different from returning null).
 */
PyObject *Matrix61c_set_value(Matrix61c *self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    if (self == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "self is null");
        return NULL;
    }

    if (!PyObject_TypeCheck(args, &PyTuple_Type)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type PyTupleObject!");
        return NULL;
    }

    // Have confirmed that args is indeed a PyTupleObjecte. Cast its type.
    PyTupleObject* pyargs = (PyTupleObject*) args;

    // TypeError if the number of arguments parsed from args is not 3
    PyObject* Pysize = PyLong_FromSsize_t(PyTuple_GET_SIZE(pyargs));
    int size = PyLong_AsLong(Pysize);
    if (size != 3) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
    
    PyObject *Pyarg_i = PyTuple_GetItem(pyargs, 0);
    PyObject *Pyarg_j = PyTuple_GetItem(pyargs, 1);
    PyObject *Pyarg_val = PyTuple_GetItem(pyargs, 2);

    // if i and j are not integers, or if val is not a float or int.
    if (!PyObject_TypeCheck(Pyarg_i, &PyLong_Type)) {
        PyErr_SetString(PyExc_TypeError, 
        "Argument must of type PyLongObject!");
        return NULL;
    }

    if (!PyObject_TypeCheck(Pyarg_j, &PyLong_Type)) {
        PyErr_SetString(PyExc_TypeError, 
        "Argument must of type PyLongObject!");
        return NULL;
    }

    if (!PyObject_TypeCheck(Pyarg_val, &PyLong_Type) 
        && !PyObject_TypeCheck(Pyarg_val, &PyFloat_Type)) {
        PyErr_SetString(PyExc_TypeError, 
        "Argument must of type PyLongObject or PyFloatObject!");
        return NULL;
    }

    int i = PyLong_AsLong(Pyarg_i);
    int j = PyLong_AsLong(Pyarg_j);
    double val;
    if (PyObject_TypeCheck(Pyarg_val, &PyLong_Type)) {
        val = (double) PyLong_AsLong(Pyarg_val);
    } else {
        val = PyFloat_AsDouble(Pyarg_val);
    }
    // Check if out-of-bound
    if (i >= self->mat->rows || j >= self->mat->cols){
        PyErr_SetString(PyExc_IndexError,
                        "add operation between two matrices with different dimensions");
        return NULL;
    }

    // Finally, setting the value
    set(self->mat, i, j, val);

    return Py_None;
}

/*
 * Given a numc.Matrix `self`, parse `args` to (int) row and (int) col.
 * Return the value at the `row`th row and `col`th column, which is a Python
 * float/int.
 */
PyObject *Matrix61c_get_value(Matrix61c *self, PyObject* args) {
   /* TODO: YOUR CODE HERE */
    if (self == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "self is null");
        return NULL;
    }

    if (!PyObject_TypeCheck(args, &PyTuple_Type)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type PyTupleObject!");
        return NULL;
    }

    // Have confirmed that args is indeed a PyTupleObjecte. Cast its type.
    PyTupleObject* pyargs = (PyTupleObject*) args;


    // TypeError if the number of arguments parsed from args is not 2
    PyObject* Pysize = PyLong_FromSsize_t(PyTuple_GET_SIZE(pyargs));
    int size = PyLong_AsLong(Pysize);
    if (size != 2) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
    
    
    PyObject *Pyarg_i = PyTuple_GetItem(pyargs, 0);
    PyObject *Pyarg_j = PyTuple_GetItem(pyargs, 1);
    

    // if i and j are not integers
    if (!PyObject_TypeCheck(Pyarg_i, &PyLong_Type)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type PyLongObject!");
        return NULL;
    }

    if (!PyObject_TypeCheck(Pyarg_j, &PyLong_Type)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type PyLongObject!");
        return NULL;
    }


    int i = PyLong_AsLong(Pyarg_i);
    int j = PyLong_AsLong(Pyarg_j);
    


    // Check if out-of-bound
    if (i >= self->mat->rows || j >= self->mat->cols){
        PyErr_SetString(PyExc_IndexError,
                        "add operation between two matrices with different dimensions");
        return NULL;
    }

    // Finally, getting the value
    double ret = get(self->mat, i, j);
    PyObject* res = PyFloat_FromDouble(ret);
    return res;
}

/*
 * Create an array of PyMethodDef structs to hold the instance methods.
 * Name the python function corresponding to Matrix61c_get_value as "get" and Matrix61c_set_value
 * as "set"
 * You might find this link helpful: https://docs.python.org/3.6/c-api/structures.html
 */
PyMethodDef Matrix61c_methods[] = {
    /* TODO: YOUR CODE HERE */
    {"get", Matrix61c_get_value, METH_VARARGS, "function Matrix61c_get_value"},
    {"set", Matrix61c_set_value, METH_VARARGS, "function Matrix61c_set_value"},
    {NULL, NULL, 0, NULL}
};

/* INDEXING */

/*
 * Given a numc.Matrix `self`, index into it with `key`. Return the indexed result.
 */
PyObject *Matrix61c_subscript(Matrix61c* self, PyObject* key) {
    /* TODO: YOUR CODE HERE */

    // Input parameter type check
    // CASE 1: 1D matrix
    if (self->mat->is_1d == 1) {
        // For a 1D matrix, the key could either be an integer or a single slice
        if (!PyObject_TypeCheck(key, &PyLong_Type) &&
            !PyObject_TypeCheck(key, &PySlice_Type)) {
            PyErr_SetString(PyExc_TypeError, "Subscript key input for 1D matrix must be type PyLong_Type or PySlice_Type!");
            return NULL;
        }

        // If 1D matrix input is an integer
        if (PyObject_TypeCheck(key, &PyLong_Type)) {
            // printf("1D case and the key is an integer\n");
            int index = PyLong_AsLong(key);
            // IndexError if key is an integer but is out of range, 
            
            
            if (self->mat->rows == 1) {
                //printf("The rows is 1\n");
                if (index >= self->mat->cols) {
                    PyErr_SetString(PyExc_IndexError, "key is out of bound!");
                    return NULL;
                }
                return PyFloat_FromDouble(get(self->mat, 0, index));
            }
            if (self->mat->cols == 1) {
                //printf("The cols is 1\n");
                if (index >= self->mat->rows) {
                    PyErr_SetString(PyExc_IndexError, "key is out of bound!");
                    return NULL;
                }
                return PyFloat_FromDouble(get(self->mat, index, 0));
            }

        } else if (PyObject_TypeCheck(key, &PySlice_Type)) {
            // If 1D matrix input is a slice
            // eg: a[1:3]: 1 is slice0 & 3 is slice1
            Py_ssize_t *slice0 = malloc(sizeof(Py_ssize_t));
            Py_ssize_t *slice1 = malloc(sizeof(Py_ssize_t));
            Py_ssize_t *step = malloc(sizeof(Py_ssize_t));
            Py_ssize_t length;
            //printf("1D case, and the input is slice\n");
            if (self->mat->cols == 1) {
                length = (Py_ssize_t)self->mat->rows;
            } else {
                length = (Py_ssize_t)self->mat->cols;
            }
            //printf("The length is %d\n", length);
            if ((int)length < 1){
                free(slice0);
                free(slice1);
                free(step);
                PyErr_SetString(PyExc_ValueError, "length < 1");
                return NULL;
            }
            //printf("Doing Pyslice_getindices\n");
            int ind_failed = PySlice_GetIndices
            (key, length, slice0, slice1, step);
            //("Pyslice_getindices is good\n");
            if (ind_failed == -1){
                free(slice0);
                free(slice1);
                free(step);
                PyErr_SetString(PyExc_RuntimeError, "indices failed");
                return NULL;
            }
            if ((int)*step != 1){
                free(slice0);
                free(slice1);
                free(step);
                PyErr_SetString(PyExc_ValueError, "step != 1");
                return NULL;
            }
            //printf("Finished checking\n");

            
            matrix *c;
            if (self->mat->rows == 1){
                //("1D case allocate_ref with slice0 %d, slice1 %d\n", (int)*slice0, (int)*slice1);
                int ref_failed = allocate_matrix_ref
                    (&c, self->mat, 0, (int)*slice0, 1, ((int)*slice1)-((int)*slice0));
                if (ref_failed){
                    PyErr_SetString(PyExc_RuntimeError, "allocate_ref failed");
                    free(slice0);
                    free(slice1);
                    free(step);
                    return NULL;
                }
                Matrix61c *final = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
                final->mat = c;
                final->shape = get_shape(c->rows, c->cols);
                free(slice0);
                free(slice1);
                free(step);
                return (PyObject*)final;
            }
            // eg: b[1:3]: 1 is slice0 & 3 is slice1
            if (self->mat->cols == 1){
                //printf("1D case allocate_ref with slice0 %d, slice1 %d\n", (int)*slice0, (int)*slice1);
                // int allocate_matrix_ref(matrix **mat, matrix *from, int row_offset, int col_offset,
                //        int rows, int cols)
                int ref_failed = allocate_matrix_ref
                    (&c, self->mat, (int)*slice0, 0, ((int)*slice1)-((int)*slice0), 1);
                if (ref_failed){
                    PyErr_SetString(PyExc_RuntimeError, "allocate_ref failed");
                    free(slice0);
                    free(slice1);
                    free(step);
                    return NULL;
                }
                Matrix61c *final = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
                final->mat = c;
                final->shape = get_shape(c->rows, c->cols);
                free(slice0);
                free(slice1);
                free(step);
                return (PyObject*)final;
            }
        } else {
            // Since the key is neigher an integer nor a slice
            // And, matrix is 1D, so it's a TypeError
            PyErr_SetString(PyExc_TypeError, "Subscript key input for 1D matrix must be type PyLong_Type or PySlice_Type!");
            return NULL;
        }

// CASE 2: 2D matrix
    } else {
        // check if key is not an integer, a slice, or a length-2 tuple of slices/ints.
        if (!PyObject_TypeCheck(key, &PyLong_Type) && !PyObject_TypeCheck(key, &PySlice_Type) &&
            PyLong_AsLong(PyLong_FromSsize_t(PyTuple_GET_SIZE(key))) != 2) {
            PyErr_SetString(PyExc_TypeError, "Argument must of type PyLong_Type or PySlice_Type!");
            return NULL;
        }


        if (PyObject_TypeCheck(key, &PyLong_Type)) {
            // Case 0 : If integer
            
            int index = PyLong_AsLong(key);
            
            // IndexError if key is an integer but is out of range, 
            if (index >= self->mat->rows) {
                PyErr_SetString(PyExc_IndexError, "key is out of bound!");
                return NULL;
            }
            // now index is our row number
            matrix *c;
            int ref_failed = allocate_matrix_ref
                    (&c, self->mat, index, 0, 1, self->mat->cols);

            if (ref_failed){
                PyErr_SetString(PyExc_RuntimeError, "allocate_ref failed");
                return NULL;
            }
            Matrix61c *final = (Matrix61c *) Matrix61c_new(&Matrix61cType, NULL, NULL);
            final->mat = c;
            final->shape = get_shape(c->rows, c->cols);
            return (PyObject*)final;
            
        } else if (PyObject_TypeCheck(key, &PySlice_Type)) {
            // If slice
            // Case 1 : key is a single slice
            // eg: a[0:2], 0 is slice0 & 2 is slice1
            
            Py_ssize_t *slice0 = malloc(sizeof(Py_ssize_t));
            Py_ssize_t *slice1 = malloc(sizeof(Py_ssize_t));;
            Py_ssize_t *step = malloc(sizeof(Py_ssize_t));
            Py_ssize_t length = (Py_ssize_t)self->mat->rows;
            
            if ((int)*slice1 - (int)*slice0 < 1){
                PyErr_SetString(PyExc_ValueError, "length < 1");
                free(slice0);
                free(slice1);
                free(step);
                return NULL;
            }

            if ((int)length < 1){
                PyErr_SetString(PyExc_ValueError, "length < 1");
                free(slice0);
                free(slice1);
                free(step);
                return NULL;
            }
            int ind_failed = PySlice_GetIndices
            (key, length, slice0, slice1, step);
            if (ind_failed){
                PyErr_SetString(PyExc_RuntimeError, "indices failed");
                free(slice0);
                free(slice1);
                free(step);
                return NULL;
            }
            if ((int)*step != 1){
                PyErr_SetString(PyExc_ValueError, "step != 1");
                free(slice0);
                free(slice1);
                free(step);
                return NULL;
            }
            
            // now we got all parameter and checked for all errors
            // use allocate ref to allocate and set the result to be
            // [row_offset:row_offset + rows, col_offset:col_offset + cols]
            matrix * c;
            int ref_failed = allocate_matrix_ref
                    (& c, self->mat, (int)*slice0, 0, ((int)*slice1)-((int)*slice0), self->mat->cols);
            if (ref_failed){
                PyErr_SetString(PyExc_RuntimeError, "allocate_ref failed");
                free(slice0);
                free(slice1);
                free(step);
                return NULL;
            }
            Matrix61c* final = (Matrix61c*)
                            Matrix61c_new(&Matrix61cType, NULL, NULL);
            final->mat = c;
            final->shape = PyTuple_Pack(2, PyLong_FromLong(c->rows), PyLong_FromLong(c->cols));
            free(slice0);
            free(slice1);
            free(step);
            return (PyObject*)final;

        } else if (PyObject_TypeCheck(key, &PyTuple_Type)){
            // if a tuple of two slices/ints
            PyObject* first;
            PyObject* second;
            
            first = PyTuple_GetItem(key, (Py_ssize_t)0);
            second = PyTuple_GetItem(key, (Py_ssize_t)1);
             
            // case 2: key is a tuple of two slices
            // eg: a[0:2, 0:2], 0 is slice0 & 2 is slice1; 0 is slice2 & 2 is slice3
            if (PySlice_Check(first) && PySlice_Check(second)){
                //printf("This is two slice\n");
                Py_ssize_t *slice0 = malloc(sizeof(Py_ssize_t));
                Py_ssize_t *slice1 = malloc(sizeof(Py_ssize_t));;
                Py_ssize_t *step = malloc(sizeof(Py_ssize_t));
                Py_ssize_t length = (Py_ssize_t)self->mat->rows;

                if ((int)length < 1){
                    PyErr_SetString(PyExc_ValueError, "length < 1");
                    free(slice0);
                    free(slice1);
                    free(step);
                    return NULL;
                }
                int ind_failed = PySlice_GetIndices
                (first, length, slice0, slice1, step);
                if (ind_failed){
                    PyErr_SetString(PyExc_RuntimeError, "indices failed");
                    free(slice0);
                    free(slice1);
                    free(step);
                    return NULL;
                }
                if ((int)*step != 1){
                    PyErr_SetString(PyExc_ValueError, "step != 1");
                    free(slice0);
                    free(slice1);
                    free(step);
                    return NULL;
                }
                
                Py_ssize_t *slice2 = malloc(sizeof(Py_ssize_t));
                Py_ssize_t *slice3 = malloc(sizeof(Py_ssize_t));
                Py_ssize_t *step2 = malloc(sizeof(Py_ssize_t));
                int ind_failed2 = PySlice_GetIndices
                (second, length, slice2, slice3, step2);
                if (ind_failed2){
                    PyErr_SetString(PyExc_RuntimeError, "indices failed");
                    free(slice2);
                    free(slice3);
                    free(step2);
                    return NULL;
                }
                if ((int)*step2 != 1){
                    PyErr_SetString(PyExc_ValueError, "step != 1");
                    free(slice2);
                    free(slice3);
                    free(step2);
                    return NULL;
                }
                
                // if both slices are 1 by 1 : we want to only return a single valu
                if (((int)*slice1)-((int)*slice0)==1 && ((int)*slice3)-((int)*slice2)){
                    //printf("This is 1x1\n");
                    int r = (int)*slice0;
                    int c = (int)*slice2;
                    free(slice0);
                    free(slice1);
                    free(slice2);
                    free(slice3);
                    free(step);
                    free(step2);
                    return PyFloat_FromDouble(get(self->mat, r, c));
                }
                // now we got all parameter and checked for all errors
                matrix *c;
                int ref_failed = allocate_matrix_ref
                (&c, self->mat, (int)*slice0, (int)*slice2, ((int)*slice1)-((int)*slice0), ((int)*slice3)-((int)*slice2));
                
                if (ref_failed){
                    PyErr_SetString(PyExc_RuntimeError, "allocate_ref failed");
                    free(slice0);
                    free(slice1);
                    free(slice2);
                    free(slice3);
                    free(step);
                    free(step2);
                    return NULL;
                }
                Matrix61c* final = (Matrix61c*)
                                Matrix61c_new(&Matrix61cType, NULL, NULL);
                final->mat = c;
                final->shape = PyTuple_Pack(2, PyLong_FromLong(c->rows), PyLong_FromLong(c->cols));
                free(slice0);
                free(slice1);
                free(slice2);
                free(slice3);
                free(step);
                free(step2);
                return (PyObject*)final;
            }
            
            // case 3: key is a tuple of (slice, int)
            // eg: a[0:2, 0], 0 is slice0 & 2 is slice1; 0 is slice2
            else if (PySlice_Check(first) && PyObject_TypeCheck(second, &PyLong_Type)){
                Py_ssize_t *slice0 = malloc(sizeof(Py_ssize_t));
                Py_ssize_t *slice1 = malloc(sizeof(Py_ssize_t));;
                Py_ssize_t *step = malloc(sizeof(Py_ssize_t));
                Py_ssize_t length = (Py_ssize_t)self->mat->rows;
                if ((int)length < 1){
                    PyErr_SetString(PyExc_ValueError, "length < 1");
                    free(slice0);
                    free(slice1);
                    free(step);
                    return NULL;
                }
                int ind_failed = PySlice_GetIndices
                (first, length, slice0, slice1, step);
                if (ind_failed){
                    PyErr_SetString(PyExc_RuntimeError, "indices failed");
                    free(slice0);
                    free(slice1);
                    free(step);
                    return NULL;
                }
                if ((int)*step != 1){
                    PyErr_SetString(PyExc_ValueError, "step != 1");
                    free(slice0);
                    free(slice1);
                    free(step);
                    return NULL;
                }
                
                // now we got all parameter and checked for all errors
                matrix *c;
                
                // we did: second = PyTuple_GetItem(key, (Py_ssize_t)1);
                int slice2 = PyLong_AsLong(second);  
                
                if (slice2 >= self->mat->cols) {
                    PyErr_SetString(PyExc_IndexError,
                        "key is out of bound");
                    return NULL;
                }
                int ref_failed = allocate_matrix_ref
                (& c, self->mat, (int)*slice0, slice2, ((int)*slice1)-((int)*slice0), slice2 + 1);
                
                if (ref_failed){
                    PyErr_SetString(PyExc_RuntimeError, "allocate_ref failed");
                    free(slice0);
                    free(slice1);
                    free(step);
                    return NULL;
                }
                Matrix61c* final = (Matrix61c*)
                                Matrix61c_new(&Matrix61cType, NULL, NULL);
                final->mat = c;
                final->shape = PyTuple_Pack(2, PyLong_FromLong(c->rows), PyLong_FromLong(c->cols));
                free(slice0);
                free(slice1);
                free(step);
                return (PyObject*)final;
            }
            
            // case 4: key is a tuple of (int, slice)
            // eg: a[0, 0:2], 0 is slice0 ; 0 is slice1 & 2 is slice2
            else if (PyObject_TypeCheck(first, &PyLong_Type) && PySlice_Check(second)){
                Py_ssize_t *slice1 = malloc(sizeof(Py_ssize_t));
                Py_ssize_t *slice2 = malloc(sizeof(Py_ssize_t));;
                Py_ssize_t *step = malloc(sizeof(Py_ssize_t));
                Py_ssize_t length = (Py_ssize_t)self->mat->rows;
                if ((int)length < 1){
                    PyErr_SetString(PyExc_ValueError, "length < 1");
                    free(slice1);
                    free(slice2);
                    free(step);
                    return NULL;
                }
                int ind_failed = PySlice_GetIndices
                (second, length, slice1, slice2, step);
                if (ind_failed){
                    PyErr_SetString(PyExc_RuntimeError, "indices failed");
                    free(slice1);
                    free(slice2);
                    free(step);
                    return NULL;
                }
                if ((int)*step != 1){
                    PyErr_SetString(PyExc_ValueError, "step != 1");
                    free(slice1);
                    free(slice2);
                    free(step);
                    return NULL;
                }
                
                // now we got all parameter and checked for all errors
                matrix *c;
                
                int slice0 = PyLong_AsLong(first);

                if (slice0 >= self->mat->rows) {
                    PyErr_SetString(PyExc_IndexError, "key is out of bound");
                    return NULL;
                }

                int ref_failed = allocate_matrix_ref
                (& c, self->mat, slice0, (int)*slice1, slice0 + 1, ((int)*slice2)-((int)*slice1));
                
                if (ref_failed){
                    PyErr_SetString(PyExc_RuntimeError, "allocate_ref failed");
                    free(slice1);
                    free(slice2);
                    free(step);
                    return NULL;
                }
                Matrix61c* final = (Matrix61c*)
                                Matrix61c_new(&Matrix61cType, NULL, NULL);
                final->mat = c;
                final->shape = PyTuple_Pack(2, PyLong_FromLong(c->rows), PyLong_FromLong(c->cols));
                free(slice1);
                free(slice2);
                free(step);
                return (PyObject*)final;
            }
            
            // case 5: key is (int, int)
            else if (PyObject_TypeCheck(first, &PyLong_Type) && PyObject_TypeCheck(second, &PyLong_Type)){
                PyObject *Pyarg_i = PyTuple_GetItem(key, 0);
                PyObject *Pyarg_j = PyTuple_GetItem(key, 1);

                // Error check
                // if i and j are not integers
                if (!PyObject_TypeCheck(Pyarg_i, &PyLong_Type)) {
                    PyErr_SetString(PyExc_TypeError, "Argument must of type PyLongObject!");
                    return NULL;
                }

                if (!PyObject_TypeCheck(Pyarg_j, &PyLong_Type)) {
                    PyErr_SetString(PyExc_TypeError, "Argument must of type     PyLongObject!");
                    return NULL;
                }

                int i = PyLong_AsLong(Pyarg_i);
                int j = PyLong_AsLong(Pyarg_j);

                if (i >= self->mat->rows || j > self->mat->cols) {
                    PyErr_SetString(PyExc_IndexError, "key is out of bound");
                    return NULL;
                }
                return PyFloat_FromDouble(get(self->mat, i, j));
            }

            else{
                // throw error
                PyErr_SetString(PyExc_TypeError, "Subscript key input for 2D matrix is wrong!");
                return NULL;
            }
        }
    }
}

/*
 * Given a numc.Matrix `self`, index into it with `key`, and set the indexed result to `v`.
 */
int Matrix61c_set_subscript(Matrix61c* self, PyObject *key, PyObject *v) {
    /* TODO: YOUR CODE HERE */

    // get the slicing result
    // !! remember that if we encounter an error along the way, we need to free slice
    Matrix61c* slice = (Matrix61c *) Matrix61c_subscript(self, key);
    //printf("Finish subscript\n");
    //printf("The slice type is \n");
    // CASE 1 : resulting slice is 1 by 1
    if (PyObject_TypeCheck(slice, &PyFloat_Type)){
        // two sub cases
        // SUB 1 : self = 1D and v = int
        //printf("Case 1\n");
        if (!PyObject_TypeCheck(v, &PyFloat_Type) 
            && !PyObject_TypeCheck(v, &PyLong_Type)){
            // type error
            PyErr_SetString(PyExc_TypeError, 
            "set: slice = single value, but key not float/long type");
            return -1;
        }
 
        if (self->mat->is_1d == 1){
            // Then key is got to be an PyLong_Type
            if (self->mat->cols == 1) {
                set(self->mat, PyLong_AsLong(key), 0, PyLong_AsLong(v));
            } else {
                set(self->mat, 0, PyLong_AsLong(key), PyLong_AsLong(v));
            }
            // return 0 if success
            return 0;
        }

        // SUB 2 : self = 2D, key is a tuple of two slices, and v is an intint
        else { 
            // a[0:1, 0:1] will return a single number
            //printf("self is 2D and key is a tuple\n");
            PyObject* first;
            PyObject* second;
            first = PyTuple_GetItem(key, (Py_ssize_t)0);
            second = PyTuple_GetItem(key, (Py_ssize_t)1);
            //printf("Finished getting first and sec\n");
            Py_ssize_t length = (Py_ssize_t)self->mat->rows;
            Py_ssize_t *slice0 = malloc(sizeof(Py_ssize_t));
            Py_ssize_t *slice1 = malloc(sizeof(Py_ssize_t));;
            Py_ssize_t *step = malloc(sizeof(Py_ssize_t));
            PySlice_GetIndices(first, length, slice0, slice1, step);
            //printf("Finished getting pyslice slice0, slice1\n");
                
            Py_ssize_t *slice2 = malloc(sizeof(Py_ssize_t));
            Py_ssize_t *slice3 = malloc(sizeof(Py_ssize_t));
            Py_ssize_t *step2 = malloc(sizeof(Py_ssize_t));
            PySlice_GetIndices(second, length, slice2, slice3, step2);
            //printf("Finished getting pyslice slice2, slice3\n");
            
            //printf("slice0 %d  slice 2%d\n", *slice0, *slice2);
            //printf("v %ld\n", PyLong_AsLong(v));
            set(self->mat, (int)*slice0, (int)*slice2, PyLong_AsLong(v));

            free(slice0);
            free(slice1);
            free(slice2);
            free(slice3);
            free(step);
            free(step2);
            return 0;
        }
    }
    
    // CASE 2 : resulting slice is a 1D matrix
    else if (slice->mat->is_1d == 1){
        //printf("This case 2\n");
        if (PyList_Check(v)){
            //printf("List check passed\n");
            int len = PyList_Size(v);
            int size;
            if (slice->mat->cols == 1) {
                size = slice->mat->rows;
            } else {
                size = slice->mat->cols;
            }
            //printf("List length %d and size %d\n", len, size);
            if (len != size){
                // value error
                PyErr_SetString(PyExc_ValueError, 
                "set: case2, v len != size of the 1D slice");
                Matrix61c_dealloc(slice);
                return -1;
            }
            
            for (int i = 0; i < len; i ++){
                //printf("index %d\n", i);
                PyObject* temp = PyList_GetItem(v, i);
                if (!PyObject_TypeCheck(temp, &PyLong_Type) 
                && !PyObject_TypeCheck(temp, &PyLong_Type)){
                    // value error
                    PyErr_SetString(PyExc_ValueError, 
                    "set: case2, 1D slice, v's element not float/long type");
                    free (slice);
                    return -1;
                }
                if (slice->mat->cols == 1) {
                    //printf("Setting %d %d with value %d \n", i, 0, PyLong_AsLong(temp));
                    set(slice->mat, i, 0, PyLong_AsLong(temp));
                } else {
                    // (0, 0) & (0, 1) <= 2
                    //printf("Setting %d %d with value %d \n", 0, i, PyLong_AsLong(temp));
                    set(slice->mat, 0, i, PyLong_AsLong(temp));
                }
            }

            Matrix61c_dealloc(slice);
            return 0;
        }
        else{
            // type error
            PyErr_SetString(PyExc_TypeError, "set: 1D slice, v must be a PyList!");
            Matrix61c_dealloc(slice);
            return -1;
        }
    }
    
    // CASE 3 : resulting slice is a 2D matrix
    else if (slice->mat->is_1d != 1){
        //printf("This is case3\n");
        if (PyList_Check(v)){
            int len = PyList_Size(v);

            // v should be a 2D list where the ith element of this list is a 1D list of integers/floats
            if (len != slice->mat->rows){
                // value error
                PyErr_SetString(PyExc_ValueError, 
                "set: case3, 2D slice, v's length != # of 2D matrix rows");
                Matrix61c_dealloc(slice);
                return -1;
            }
            for (int i = 0; i < len; i ++){
                PyObject* temp = PyList_GetItem(v, i);
                
                if (!PyList_Check(temp)){
                    // value error
                    PyErr_SetString(PyExc_ValueError, 
                    "set: case3, 2D slice, v is not a 2D list");
                    Matrix61c_dealloc(slice);
                    return -1;
                }
                
                int subLen = PyList_Size(temp);
                if (subLen != slice->mat->cols){
                    // value error
                    PyErr_SetString(PyExc_ValueError, 
                    "set: case3, 2D slice, v's element's len != # of 2D matrix cols");
                    Matrix61c_dealloc(slice);
                    return -1;
                }
        
                for (int j = 0; j < subLen; j ++){
                    PyObject* temp2 = PyList_GetItem(temp, j);
                    if (!PyObject_TypeCheck(temp2, &PyLong_Type)
                     && !PyObject_TypeCheck(temp2, &PyFloat_Type)){
                        // value error
                        PyErr_SetString(PyExc_ValueError, 
                        "set: case3, 2D slice, v's element's element is not float/long type");
                        Matrix61c_dealloc(slice);
                        return -1;
                    }
                    set(slice->mat, i, j, PyLong_AsLong(temp2));
                }
            }
            Matrix61c_dealloc(slice);
            return 0;
        }
    } else{
            //printf("Nothing\n");
            // type error
            PyErr_SetString(PyExc_RuntimeError, 
            "slicing result is not an int or a matrix");
            Matrix61c_dealloc(slice);
            return -1;
        }
    }

PyMappingMethods Matrix61c_mapping = {
    NULL,
    (binaryfunc) Matrix61c_subscript,
    (objobjargproc) Matrix61c_set_subscript,
};

/* INSTANCE ATTRIBUTES*/
PyMemberDef Matrix61c_members[] = {
    {
        "shape", T_OBJECT_EX, offsetof(Matrix61c, shape), 0,
        "(rows, cols)"
    },
    {NULL}  /* Sentinel */
};

PyTypeObject Matrix61cType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "numc.Matrix",
    .tp_basicsize = sizeof(Matrix61c),
    .tp_dealloc = (destructor)Matrix61c_dealloc,
    .tp_repr = (reprfunc)Matrix61c_repr,
    .tp_as_number = &Matrix61c_as_number,
    .tp_flags = Py_TPFLAGS_DEFAULT |
    Py_TPFLAGS_BASETYPE,
    .tp_doc = "numc.Matrix objects",
    .tp_methods = Matrix61c_methods,
    .tp_members = Matrix61c_members,
    .tp_as_mapping = &Matrix61c_mapping,
    .tp_init = (initproc)Matrix61c_init,
    .tp_new = Matrix61c_new
};


struct PyModuleDef numcmodule = {
    PyModuleDef_HEAD_INIT,
    "numc",
    "Numc matrix operations",
    -1,
    Matrix61c_class_methods
};

/* Initialize the numc module */
PyMODINIT_FUNC PyInit_numc(void) {
    PyObject* m;

    if (PyType_Ready(&Matrix61cType) < 0)
        return NULL;

    m = PyModule_Create(&numcmodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&Matrix61cType);
    PyModule_AddObject(m, "Matrix", (PyObject *)&Matrix61cType);
    printf("CS61C Fall 2020 Project 4: numc imported!\n");
    fflush(stdout);
    return m;
}
