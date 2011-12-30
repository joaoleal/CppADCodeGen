/* --------------------------------------------------------------------------
CppAD: C++ Algorithmic Differentiation: Copyright (C) 2011 Ciengis

CppAD is distributed under multiple licenses. This distribution is under
the terms of the 
                    Common Public License Version 1.0.

A copy of this license is included in the COPYING file of this distribution.
Please visit http://www.coin-or.org/CppAD/ for information on other licenses.
-------------------------------------------------------------------------- */

#ifndef CPPAD_CODEGEN_MATRIX_INCLUDED
#define	CPPAD_CODEGEN_MATRIX_INCLUDED

#include <algorithm>

CPPAD_BEGIN_NAMESPACE

template<class Type>
class Row;

template<class Type>
class Matrix {
private:
    /// maximum number of elements that should ever be in this vector
    const static size_t MAX_CAPACITY_;
    // array containing the data
    Type* data_;
    std::vector<Row<Type> > rows_;
    /// maximum number of Type elements current allocation can hold
    size_t capacity_;
    /// current number of rows
    size_t m_;
    /// current number of columns
    size_t n_;
public:

    Matrix() {
        m_ = 0;
        n_ = 0;
        capacity_ = 0;
        data_ = NULL;
    }

    Matrix(size_t m, size_t n) {
        CPPAD_ASSERT_UNKNOWN(m > MAX_CAPACITY_ / n); // check for overflow
        m_ = m;
        n_ = n;
        capacity_ = m_ * n_;
        data_ = new Type[m * n];
        rows_.resize(m);

        for (size_t i = 0; i < m; i++) {
            rows_[i] = Row<Type > (this, i);
        }

        // set the default data;
        std::fill(data_, data_[capacity_], Type());
    }

    inline size_t rows() const {
        return m_;
    }

    inline size_t columns() const {
        return n_;
    }

    inline size_t capacity() const {
        return capacity_;
    }

    /// changes the current number of rows

    inline size_t resizeRows(size_t m) {
        if (m_ == m) {
            return;
        }

        if (m > m_) {
            CPPAD_ASSERT_UNKNOWN(m > MAX_CAPACITY_ / n_); // check for overflow

            size_t min_cap = m*n_;
            if (min_cap < capacity_) {
                // resize data

                // save old information
                size_t old_m = m_;
                Type* old_data = data_;
                size_t old_cap = capacity_;
                capacity_ = m * n_;

                data_ = new Type[capacity_];

                // move data
                std::copy(old_data, old_data[old_cap], data_);
                for (size_t i = 0; i < old_m; i++) {
                    rows_[i].resetRowIndex(i);
                }
                rows_.resize(m);
                for (size_t i = old_m; i < m; i++) {
                    rows_[i] = Row<Type > (this, i);
                }

                // delete old data
                delete [] old_data;
            }

            // set default data in the new region
            std::fill(data_[m_ * n_], data_[m * n_], Type());
        }

        m_ = m;
    }

    /*!
    Remove all the elements from this matrix but leave the capacity
    and data pointer as is.
     */
    void erase() {
        m_ = 0;
        n_ = 0;
    }

    // ----------------------------------------------------------------------
    /// non-constant element access; i.e., we can change this element value

    inline Row<Type>& operator[](
            /// row index, must be less than the number of rows
            size_t i
            ) {
        CPPAD_ASSERT_UNKNOWN(i < m_);
        return rows_[i];
    }
    // ----------------------------------------------------------------------
    /// constant element access; i.e., we cannot change this element value

    inline const Row<Type>& operator[](
            /// row index, must be less than the number of rows
            size_t i
            ) const {
        CPPAD_ASSERT_UNKNOWN(i < m_);
        return rows_[i];
    }

    virtual ~Matrix() {
        delete [] data_;
    }

private:

    explicit Matrix(const Matrix& other) {
        CPPAD_ASSERT_UNKNOWN(false); // not implemented yet
    }

    void operator=(const Matrix& rhs) {
        CPPAD_ASSERT_UNKNOWN(false); // not implemented yet
    }

    template<class T>
    friend class Row;
};

template<class Type>
const size_t Matrix<Type>::MAX_CAPACITY_ = std::numeric_limits<size_t>::max();

template<class Type>
class Row {
private:
    const Matrix<Type>* matrix_; // the matrix where this row is located
    Type* data_; // the row data

private:

    Row() {
        matrix_ = NULL;
        data_ = NULL;
    }

    Row(const Matrix<Type>* matrix, size_t row) {
        matrix_ = matrix;
        data_ = &matrix_->data_[row * matrix->columns()];
    }

    inline void resetRowIndex(size_t row) {
        data_ = &matrix_->data_[row * matrix_->columns()];
    }
public:
    // ----------------------------------------------------------------------
    /// non-constant element access; i.e., we can change this element value

    inline Type& operator[](
            /// element index, must be less than the number of columns
            size_t j
            ) {
        CPPAD_ASSERT_UNKNOWN(j < matrix_->n_);
        return data_[j];
    }
    // ----------------------------------------------------------------------
    /// constant element access; i.e., we cannot change this element value

    inline const Type& operator[](
            /// element index, must be less than the number of columns
            size_t j
            ) const {
        CPPAD_ASSERT_UNKNOWN(j < matrix_->n_);
        return data_[j];
    }

    virtual ~Row() {
        // nothing to do
    }
};


CPPAD_END_NAMESPACE

#endif	/* CPPAD_CODEGEN_MATRIX_INCLUDED */

