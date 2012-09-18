#ifndef MATRIX_OP_HPP
#define MATRIX_OP_HPP
// Copyright Fang-Cheng Yeh 2010
// Distributed under the BSD License
//
/*
Copyright (c) 2010, Fang-Cheng Yeh
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#include <iterator>
#include <cmath>
#include <vector>
#include <algorithm>

namespace math
{


template<typename lhs_type,typename rhs_type>
typename std::iterator_traits<lhs_type>::value_type
vector_op_dot(lhs_type v1,lhs_type v1_end,rhs_type v2)
{
    typedef typename std::iterator_traits<lhs_type>::value_type value_type;
    if (v1 == v1_end)
        return value_type(0);
    value_type sum((*v1)*(*v2));
    if (++v1 != v1_end)
        do
        {
            ++v2;
            sum += (*v1)*(*v2);
        }
        while (++v1 != v1_end);
    return sum;
}

/*
perform x <- -x
*/
template<typename lhs_type>
void vector_op_negate(lhs_type x_begin,lhs_type x_end)
{
    for (;x_begin != x_end;++x_begin)
        (*x_begin) = -(*x_begin);
}

/*
perform x <- ax
*/
template<typename lhs_type,typename scalar_type>
void vector_op_scale(lhs_type x_begin,lhs_type x_end,scalar_type a)
{
    for (;x_begin != x_end;++x_begin)
        (*x_begin) *= a;
}

/*
perform y <- ax
*/
template<typename lhs_type,typename rhs_type,typename scalar_type>
void vector_op_scale(lhs_type x_begin,lhs_type x_end,rhs_type y,scalar_type a)
{
    if (x_begin != x_end)
        do
        {
            (*y) = (*x_begin)*a;
            if (++x_begin == x_end)
                return;
            ++y;
        }
        while (1);
}

/*
calculate norm2
*/
template<typename lhs_type>
typename std::iterator_traits<lhs_type>::value_type
vector_op_norm2(lhs_type x_begin,lhs_type x_end)
{
    typedef typename std::iterator_traits<lhs_type>::value_type value_type;
    // dimension = 0
    if (x_begin == x_end)
        return value_type(0);

    // dimension = 1
    value_type sum(*x_begin);
    if (++x_begin == x_end)
        return std::abs(sum);

    // dimension >= 2
    sum *= sum;
    value_type x2(*x_begin);
    x2 *= x2;
    sum += x2;
    while (++x_begin != x_end)
    {
        x2 = *x_begin;
        x2 *= x2;
        sum += x2;
    }
    return std::sqrt(sum);
}

/*
swap vector x, y
*/
template<typename lhs_type,typename rhs_type>
void vector_op_swap(lhs_type x_begin,lhs_type x_end,rhs_type y)
{
    if (x_begin != x_end)
        do
        {
            std::swap(*x_begin,*y);
            if (++x_begin == x_end)
                return;
            ++y;
        }
        while (1);
}

/*
perform y <- ax+y
*/
template<typename lhs_type,typename rhs_type,typename scalar_type>
void vector_op_axpy(lhs_type y_begin,lhs_type y_end,scalar_type a,rhs_type x)
{
    if (y_begin != y_end)
        do
        {
            *y_begin += (*x)*a;
            if (++y_begin == y_end)
                return;
            ++x;
        }
        while (1);
}


/*
perform x <- ay+x
*/
template<typename lhs_type,typename rhs_type,typename scalar_type>
void vector_op_aypx(lhs_type y_begin,lhs_type y_end,scalar_type a,rhs_type x)
{
    if (y_begin != y_end)
        do
        {
            (*x) += (*y_begin)*a;
            if (++y_begin == y_end)
                return;
            ++x;
        }
        while (1);
}



/*
perform x <- c*x+s*y
perform y <- c*y-s*x
*/
template<typename lhs_type,typename rhs_type,typename scalar_type>
void vector_op_rot(lhs_type x_begin,lhs_type x_end,rhs_type y,scalar_type c,scalar_type s)
{
    typename std::iterator_traits<lhs_type>::value_type x_temp;
    if (x_begin != x_end)
        do
        {

            x_temp = (*x_begin)*c + (*y)*s;
            *y     = (*y)*c       - (*x_begin)*s;
            *x_begin = x_temp;

            if (++x_begin == x_end)
                return;
            ++y;
        }
        while (1);
}



/*
perform A=x*y'
*/

template<typename left_input_iterator,
typename right_input_iterator,
typename output_iterator>
void vector_op_gen(left_input_iterator x,left_input_iterator x_end,right_input_iterator y,output_iterator out)
{
    size_t dim = (x_end-x);
    right_input_iterator y_end = y + dim;
    if (x != x_end)
        do
        {
            vector_op_scale(y,y_end,out,*x);
            if (++x == x_end)
                return;
            out += dim;
        }
        while (1);
}

/*
perform A=x*y'
*/

template<typename left_input_iterator,
typename right_input_iterator,
typename output_iterator>
void vector_op_gen(left_input_iterator x,left_input_iterator x_end,right_input_iterator y,right_input_iterator y_end,output_iterator out)
{
    size_t dim = (y_end-y);
    if (x != x_end)
        do
        {
            vector_op_scale(y,y_end,out,*x);
            if (++x == x_end)
                return;
            out += dim;
        }
        while (1);
}

// data type for specifying the matrix dimension in compiler time
template<size_t row,size_t col>
struct dim
{
    size_t row_count(void)const
    {
        return row;
    }
    size_t col_count(void)const
    {
        return col;
    }
    size_t size(void)const
    {
        return row*col;
    }
};

// data type for specifying the matrix dimension in runtime

struct dyndim
{
    size_t row,col;
    dyndim(void) {}
    dyndim(size_t row_):row(row_),col(1) {}
    dyndim(size_t row_,size_t col_):row(row_),col(col_) {}
    size_t row_count(void)const
    {
        return row;
    }
    size_t col_count(void)const
    {
        return col;
    }
    size_t size(void)const
    {
        return row*col;
    }
};


/*
perform y = Ax
*/

template<typename left_input_iterator,
typename right_input_iterator,
typename output_iterator,
typename left_dim_type>
void matrix_vector_product(left_input_iterator lhs,right_input_iterator x,output_iterator y,const left_dim_type& ldim)
{
    left_input_iterator lhs_end = lhs + ldim.size();
    if (lhs == lhs_end)
        return;
    size_t common_col_count = ldim.col_count();
    left_input_iterator lhs_next;
    do
    {
        lhs_next = lhs + common_col_count;
        *y = vector_op_dot(lhs,lhs_next,x);
        if (lhs_next == lhs_end)
            return;
        lhs = lhs_next;
        ++y;
    }
    while (1);
}


/**
perform A*B

INPUT: must be random access iterator
OUTPUT: random access iterator or bidirectional iterator

*/


template<typename left_input_iterator,
typename right_input_iterator,
typename output_iterator,
typename left_dim_type,
typename right_dim_type>
void matrix_product(left_input_iterator lhs			    /*A*/,
                    right_input_iterator rhs			/*B*/,
                    output_iterator out					/*output*/,
                    const left_dim_type& ldim			/* the dimension of A*/,
                    const right_dim_type& rdim			/* the dimension of B*/)
{
    size_t common_col_count = ldim.col_count();
    size_t right_col_count = rdim.col_count();
    left_input_iterator lhs_end = lhs + ldim.size();
    right_input_iterator rhs_end = rhs + rdim.col_count();

    while (lhs != lhs_end)
    {
        left_input_iterator lhs_to = lhs + common_col_count;
        for (right_input_iterator rhs_col = rhs;rhs_col != rhs_end;++rhs_col,++out)
        {
            right_input_iterator rhs_from = rhs_col;
            left_input_iterator lhs_from = lhs;
            typename std::iterator_traits<left_input_iterator>::value_type sum((*lhs_from)*(*rhs_from));
            if (++lhs_from != lhs_to)
                do
                {
                    rhs_from += right_col_count;
                    sum += (*lhs_from)*(*rhs_from);
                }
                while (++lhs_from != lhs_to);
            *out = sum;
        }
        lhs = lhs_to;
    }
}

/**
perform A*Bt

INPUT: must be random access iterator
OUTPUT: random access iterator or bidirectional iterator

*/
template<typename left_input_iterator,
typename right_input_iterator,
typename output_iterator,
typename left_dim_type,
typename right_dim_type>
void matrix_product_transpose(
    left_input_iterator lhs			    /*A*/,
    right_input_iterator rhs			/*B*/,
    output_iterator out					/*output*/,
    const left_dim_type& ldim			/* the dimension of A*/,
    const right_dim_type& rdim			/* the dimension of B*/)
{
    size_t common_col_count = ldim.col_count();
    left_input_iterator lhs_end = lhs + ldim.size();
    right_input_iterator rhs_end = rhs + rdim.size();

    for (;lhs != lhs_end;lhs += common_col_count)
        for (right_input_iterator rhs_iter = rhs;rhs_iter != rhs_end;rhs_iter += common_col_count,++out)
            *out = vector_op_dot(lhs,lhs+common_col_count,rhs_iter);

}


/*
perform A*At
INPUT: must be random access iterator
OUTPUT: must be random access iterator
*/

template<typename input_iterator,
typename output_iterator,
typename dim_type>
void matrix_square(input_iterator lhs,output_iterator out,const dim_type& dim)
{
    output_iterator iter = out;

    size_t common_col_count = dim.col_count();
    input_iterator rhs = lhs;
    input_iterator lhs_end = lhs + dim.size();
    input_iterator rhs_end = rhs + dim.size();
    for (size_t row = 0;lhs != lhs_end;lhs += common_col_count,++row)
    {
        input_iterator rhs_iter = rhs;
        for (size_t col = 0;rhs_iter != rhs_end;rhs_iter += common_col_count,++out,++col)
            if (row <= col)// skip the symmetric part
                *out = vector_op_dot(lhs,lhs+common_col_count,rhs_iter);
    }

    size_t row_count = dim.row_count();
    if (row_count > 1)
    {
        input_iterator col_wise = iter + 1;
        input_iterator row_wise = iter + row_count;
        size_t shift = row_count + 1;
        for (size_t length = row_count - 1;1;col_wise += shift,row_wise += shift)
        {
            input_iterator col_from = col_wise;
            input_iterator col_to = col_wise + length;
            input_iterator row_from = row_wise;
            while (1)
            {
                *row_from = *col_from;
                if (++col_from == col_to)
                    break;
                row_from += row_count;
            }
            if (--length <= 0)
                break;
        }
    }
}







/**
example:
\code
double sym[]={12, 8, 3, 1,
               8, 4, 2, 5,
               3, 2,11, 4,
               1, 5, 4, 7};
is_symmetric(sym,dim<4,4>());
\endcode
*/
template<typename input_iterator,typename dim_type>
bool matrix_is_symmetric(input_iterator iter,const dim_type& dim)
{
    input_iterator col_wise = iter + 1;
    input_iterator row_wise = iter + dim.col_count();
    size_t row_count = dim.row_count();
    size_t shift = dim.col_count() + 1;
    for (size_t length = dim.col_count() - 1;length > 0;--length)
    {
        input_iterator col_from = col_wise;
        input_iterator col_to = col_wise + length;
        input_iterator row_from = row_wise;
        while (1)
        {
            if (*col_from != *row_from)
                return false;
            ++col_from;
            if (col_from == col_to)
                break;
            row_from += row_count;
        }
        col_wise += shift;
        row_wise += shift;
    }
    return true;
}

/**
example:
\code
double sym[]={12, 8, 3, 1,
       9, 4, 2, 5,
       5, 3,11, 4,
       2, 5, 6, 7};
matrix_transpose(sym,dim<4,4>());
\endcode
*/
template<typename input_iterator,size_t matrix_dim>
void matrix_inplace_transpose(input_iterator A,dim<matrix_dim,matrix_dim>)
{
    if (matrix_dim > 1)
    {
        input_iterator col_wise = A + 1;
        input_iterator row_wise = A + matrix_dim;
        size_t shift = matrix_dim + 1;
        for (size_t length = matrix_dim - 1;1;col_wise += shift,row_wise += shift)
        {
            input_iterator col_from = col_wise;
            input_iterator col_to = col_wise + length;
            input_iterator row_from = row_wise;
            while (1)
            {
                std::swap(*col_from,*row_from);
                if (++col_from == col_to)
                    break;
                row_from += matrix_dim;
            }
            if (--length <= 0)
                break;
        }
    }
}



template<typename input_iterator,typename output_iterator,typename dim_type>
void matrix_transpose(input_iterator in,output_iterator out,const dim_type& dim)
{
    size_t col = 0;
    size_t out_leap = dim.row_count();
    output_iterator out_col = out + col;
    output_iterator out_end = out + dim.size()-out_leap;// last leap position
    for (input_iterator end = in + dim.size();in != end;++in)
    {
        *out_col = *in;
        if (out_col >= out_end)
        {
            ++col;
            out_col = out + col;
        }
        else
            out_col += out_leap;
    }
}

template<typename io_iterator,typename dim_type>
void matrix_inplace_transpose(io_iterator io,const dim_type& dim)
{
    typedef typename std::iterator_traits<io_iterator>::value_type value_type;
    std::vector<value_type> temp(io,io+dim.size());
    matrix_transpose(temp.begin(),io,dim);
}


template<typename input_iterator,typename value_type>
void matrix_col_rotate_dyn(input_iterator col1,input_iterator col2,
                           value_type c,value_type s,size_t row_count,size_t col_count)
{
    value_type temp;
    size_t row = 0;
    do
    {
        temp = *col2;
        *col2 = s * (*col1) + c * temp;
        *col1 = c * (*col1) - s * temp;
        if (++row == row_count)
            break;
        col1 += col_count;
        col2 += col_count;
    }
    while (1);
}



template <typename iterator_type,typename dim_type>
void matrix_identity(iterator_type I,const dim_type& dim)
{
    size_t size = dim.size();
    typedef typename std::iterator_traits<iterator_type>::value_type value_type;
    std::fill(I,I+size,value_type(0));
    size_t leap_size = dim.col_count()+1;
    for (size_t index = 0;index < size;index += leap_size)
        I[index] = value_type(1);
}

/**

	double A[]={8, 1, 3,
				7, 0, 2,
                12, 3, 2};
	size_t pivot[3];
	math::matrix_lu_decomposition(A,pivot,math::dim<3,3>());

*/
template<typename io_iterator,typename pivot_iterator,typename dim_type>
bool matrix_lu_decomposition(io_iterator A,pivot_iterator pivot,const dim_type& dim)
{
    typedef typename std::iterator_traits<io_iterator>::value_type value_type;
    const size_t dimension = dim.row_count();
    const size_t size = dim.size();
    for (size_t k = 0;k < dimension;++k)
        pivot[k] = k;
    for (size_t k = 0,row_k = 0;k < dimension;++k,row_k+=dimension)
    {
        // pivoting
        {
            value_type max_value(0);
            size_t max_index = 0;
            size_t max_row = k;
            for (size_t i = k,index_ik = row_k + k;i < dimension;++i,index_ik += dimension)
            {
                value_type value = A[index_ik];
                if (value < 0)
                    value = -value;
                if (value > max_value)
                {
                    max_value = value;
                    max_index = index_ik;
                    max_row = i;
                }
            }
            if (max_value == 0)
                return false; // singularity
            if (max_row != k) // row swap is needed
            {
                vector_op_swap(A+row_k,A+row_k+dim.col_count(),A+max_index-k);
                std::swap(pivot[k],pivot[max_row]);
            }
        }

        //  reduce the matrix
        value_type bjj = A[row_k + k];
        for (size_t row_i = row_k + dimension;row_i < size;row_i += dimension)
        {
            value_type temp = A[row_i + k] /= bjj;
            for (size_t j = k+1;j < dimension;++j)
                A[row_i + j] -= temp*A[row_k + j];
        }
    }
    return true;
}


template<typename input_iterator,typename dim_type>
typename std::iterator_traits<input_iterator>::value_type
matrix_lu_determinant(input_iterator A,const dim_type& dim)
{
    typedef typename std::iterator_traits<input_iterator>::value_type value_type;
    value_type result = A[0];
    size_t leap_size = dim.col_count()+1;
    for (size_t index = leap_size;index < dim.size();index += leap_size)
        result *= A[index];
    return result;
}

template<typename input_iterator,size_t count>
typename std::iterator_traits<input_iterator>::value_type
matrix_determinant(input_iterator iter,dim<count,count>)
{
    //not_yet_implemented();
}

template<typename input_iterator>
typename std::iterator_traits<input_iterator>::value_type
matrix_determinant(input_iterator iter,dim<4,4>)
{
    typedef typename std::iterator_traits<input_iterator>::value_type value_type;
    value_type v_9_14_m_10_13  = (*(iter+ 9))*(*(iter+14))-(*(iter+10))*(*(iter+13));
    value_type v_10_15_m_11_14 = (*(iter+10))*(*(iter+15))-(*(iter+11))*(*(iter+14));
    value_type v_11_12_m_8_15 = (*(iter+11))*(*(iter+12))-(*(iter+ 8))*(*(iter+15));
    value_type v_8_13_m_9_12 = ((*(iter+ 8))*(*(iter+13))-(*(iter+ 9))*(*(iter+12)));
    value_type v_11_13_m_9_15 = ((*(iter+11))*(*(iter+13))-(*(iter+ 9))*(*(iter+15)));
    value_type v_8_14_m_10_12 = ((*(iter+ 8))*(*(iter+14))-(*(iter+10))*(*(iter+12)));
    return
        (*iter)*
        (
            (*(iter+ 5))*v_10_15_m_11_14+
            (*(iter+ 6))*v_11_13_m_9_15+
            (*(iter+ 7))*v_9_14_m_10_13
        )-
        *(iter+1)*
        (
            (*(iter+ 4))*v_10_15_m_11_14+
            (*(iter+ 6))*v_11_12_m_8_15+
            (*(iter+ 7))*v_8_14_m_10_12
        )+
        *(iter+2)*
        (
            (*(iter+ 4))*(-v_11_13_m_9_15)+
            (*(iter+ 5))*v_11_12_m_8_15+
            (*(iter+ 7))*v_8_13_m_9_12
        )-
        *(iter+3)*
        (
            (*(iter+ 4))*v_9_14_m_10_13+
            (*(iter+ 5))*(-v_8_14_m_10_12)+
            (*(iter+ 6))*v_8_13_m_9_12
        );
}


/**
example:
\code
double sym[]={12, 8, 3,
               8, 4, 2,
               3, 2,11};
std::cout << la::matrix_determinant(sym,la::dim<3,3>());
\endcode
*/
template<typename input_iterator>
typename std::iterator_traits<input_iterator>::value_type
matrix_determinant(input_iterator iter,dim<3,3>)
{
    return (*(iter  ))*((*(iter+4))*(*(iter+8))-(*(iter+5))*(*(iter+7)))+
           (*(iter+1))*((*(iter+5))*(*(iter+6))-(*(iter+3))*(*(iter+8)))+
           (*(iter+2))*((*(iter+3))*(*(iter+7))-(*(iter+4))*(*(iter+6)));
}

template<typename input_iterator>
typename std::iterator_traits<input_iterator>::value_type
matrix_determinant(input_iterator iter,dim<2,2>)
{
    return (*iter)*(*(iter+3)) - (*(iter+1))*(*(iter+2));
}
/**
    solve Ax=b

*/
template<typename input_iterator1,typename input_iterator2,typename piv_iterator,typename output_iterator,typename dim_type>
bool matrix_lu_solve(input_iterator1 A,piv_iterator piv,input_iterator2 b,output_iterator x,const dim_type& dim)
{
    typedef typename std::iterator_traits<input_iterator1>::value_type value_type;
    const size_t col_count = dim.col_count();
    const size_t matrix_size = dim.size();
    for (size_t i = 0; i < col_count;++i)
        *(x+i) = b[piv[i]];
    // Solve L*Y = B(piv)
    {
        for (size_t j = 0, k = 0;j < matrix_size;j += col_count,++k)
        {
			value_type x_k = *(x+k); 
            for (size_t i = j + col_count + k, m = k+1;i < matrix_size;i += col_count,++m)
                *(x+m) -= x_k*(A[i]);  // A[i][k]
        }
    }
    // Solve U*X = Y;
    {
        size_t j = matrix_size - 1;
        size_t diagonal_shift = col_count + 1;
        for (int k = dim.row_count()-1;k >= 0;--k,j -= diagonal_shift)
        {
            value_type Arowk_value = A[j];
            if (Arowk_value + value_type(1) == value_type(1))
                return false;
            value_type x_k(*(x+k) /= Arowk_value);
            input_iterator1 Arowi = A + k;
            for (size_t i = 0;i < k;++i,Arowi += col_count)
                *(x+i) -= x_k*(*Arowi);
        }
    }
    return true;
}

/**
    solve AX=B

*/
template<typename input_iterator1,typename input_iterator2,typename piv_iterator,typename output_iterator,typename dim_type>
bool matrix_lu_solve(input_iterator1 A,piv_iterator piv,input_iterator2 B,output_iterator X,const dim_type& dim,const dim_type& Bdim)
{
    typedef typename std::iterator_traits<input_iterator1>::value_type value_type;
    std::vector<value_type> b(Bdim.row_count());
    std::vector<value_type> x(Bdim.row_count());
    bool result = true;
    for (size_t col = 0;col < Bdim.col_count();++col)
    {
        for (size_t row = 0,index = col;row < Bdim.row_count();++row,index += Bdim.col_count())
            b[row] = B[index];
        if (matrix_lu_solve(A,piv,&*b.begin(),&*x.begin(),dim))
            for (size_t row = 0,index = col;row < Bdim.row_count();++row,index += Bdim.col_count())
                X[index] = x[row];
        else
            result = false;
    }
    return result;
}



/**
example:
\code
double sym[]={12, 8,
              9, 4};
matrix_inverse(sym,dim<2,2>());
\endcode
*/
template<typename input_iterator>
bool matrix_inverse(input_iterator iter,dim<2,2>)
{
    typedef typename std::iterator_traits<input_iterator>::value_type value_type;
    value_type det = matrix_determinant(iter,dim<2,2>());
    if (det+value_type(1) == value_type(1))
        return false;
    std::swap(*(iter),*(iter+3));
    *(iter  ) /= det;
    *(iter+1) /= -det;
    *(iter+2) /= -det;
    *(iter+3) /= det;
    return true;
}

template<typename input_iterator>
bool matrix_inverse(input_iterator iter,dim<3,3>)
{
    typedef typename std::iterator_traits<input_iterator>::value_type value_type;
    value_type det = matrix_determinant(iter,dim<3,3>());
    if (det+value_type(1) == value_type(1))
        return false;
    value_type temp[9];
    temp[0] = iter[4]*iter[8]-iter[5]*iter[7];
    temp[1] = iter[2]*iter[7]-iter[1]*iter[8];
    temp[2] = iter[1]*iter[5]-iter[2]*iter[4];
    temp[3] = iter[5]*iter[6]-iter[3]*iter[8];
    temp[4] = iter[0]*iter[8]-iter[2]*iter[6];
    temp[5] = iter[2]*iter[3]-iter[0]*iter[5];
    temp[6] = iter[3]*iter[7]-iter[4]*iter[6];
    temp[7] = iter[1]*iter[6]-iter[0]*iter[7];
    temp[8] = iter[0]*iter[4]-iter[1]*iter[3];
    iter[0] = temp[0]/det;
    iter[1] = temp[1]/det;
    iter[2] = temp[2]/det;
    iter[3] = temp[3]/det;
    iter[4] = temp[4]/det;
    iter[5] = temp[5]/det;
    iter[6] = temp[6]/det;
    iter[7] = temp[7]/det;
    iter[8] = temp[8]/det;
    return true;
}

template<typename input_iterator,typename dim_type>
bool matrix_inverse(input_iterator A_,dim_type dim)
{
	typedef typename std::iterator_traits<input_iterator>::value_type value_type;
    const size_t dimension = dim.row_count();
    const size_t matrix_size = dim.size();
	std::vector<value_type> buf(A_,A_+matrix_size);
	std::vector<unsigned int> piv(dimension);
    value_type* A = &(buf[0]);
    if(!matrix_lu_decomposition(A,piv.begin(),dim))
		return false;

    bool result = true;
    std::vector<value_type> x_buf(dimension);
    value_type* x = &(x_buf[0]);    
	for (size_t col = 0;col < dimension && result;++col)
    {
		
		// Solve L*Y = B(piv)
        {
			bool find = false;
            for (size_t j = 0, k = 0;j < matrix_size;j += dimension,++k)
            {
				if(find)
				{
					value_type x_k = *(x+k); 
            		for (size_t i = j + dimension + k, m = k+1;i < matrix_size;i += dimension,++m)
						*(x+m) -= x_k*(A[i]);  // A[i][k]
				}
				else
				{
					if(piv[k] != col)
						continue;
					*(x+k) = 1.0;
					for (size_t i = j + dimension + k, m = k+1;i < matrix_size;i += dimension,++m)
						*(x+m) -= A[i];  // A[i][k]
					find = true;
				}
            }
        }
        // Solve U*X = Y;
        {
            size_t j = matrix_size - 1;
            size_t diagonal_shift = dimension + 1;
            for (int k = dimension-1;k >= 0;--k,j -= diagonal_shift)
            {
                value_type Arowk_value = A[j];
                if (Arowk_value + value_type(1) == value_type(1))
				{
					result = false;
					break;
				}
                value_type x_k(*(x+k) /= Arowk_value);
                value_type* Arowi = A + k;
                for (size_t i = 0;i < k;++i,Arowi += dimension)
                    *(x+i) -= x_k*(*Arowi);
            }
        }

        for (size_t row = 0,index = col;row < dimension;++row,index += dimension)
			A_[index] = x[row];   
		std::fill(x,x+dimension,0.0);
    }
    return result;
}

template<typename input_iterator,typename output_iterator,typename dim_type>
bool matrix_inverse(input_iterator A_,output_iterator A,dim_type dim)
{
    std::copy(A_,A_+dim.size(),A);
    return matrix_inverse(A,dim);
}



/** Apply housholder reduction to make column (c) to be 0 below row (r)
		  |  a'				  |				a'
		A=|	   ,a1,a2,a3..	an|		, a0 = | |
		  |  x  			  |				x

	P=1-u*uT/uu_2

	u = x-|x|e0

	Px = |x|e0,

		    |I	0|    | a'| |  a00  |
	make P'=|	 |, P'|   |=|	    |
		    |0  P|    | x | | |x|e0 |
					|  a'				 |
	therefore P'A = |		,P'a1,P'a2...|
	                | |x|e0				 |

	result:
		vector u was stored in the place of x
		right side columns were changed due to P'
	*/
/*
template<typename input_iterator,typename dim_type>
typename std::iterator_traits<input_iterator>::value_type
household_col_reduction(input_iterator row,size_t row_index,size_t col_index,const dim_type& dim)
{
typedef typename std::iterator_traits<input_iterator>::value_type value_type;
// calculate all the parameter of Housholder reduction
input_iterator x = row + col_index;
value_type x2 = matrix_col_vector_length2(row_iterator(row.index(),x));
if (x2+value_type(1) == value_type(1)) // |x|=0
    return 0;
value_type x_norm = std::sqrt(x2);
//store u in the original place of x
// if x0 > 0, u = x-( |x|)e0, => Px = ( |x|)e0; u^2=2(|x|^2-( |x|)x0)
//    x0 < 0, u = x-(-|x|)e0, => Px = (-|x|)e0; u^2=2(|x|^2-(-|x|)x0)
//
if (x[0] < 0)
    x_norm = -x_norm;
value_type uu_2 = x2 - x_norm * x[0];
x[0] -= x_norm;

// calculate P'a1,P'a2..., upper part is not change for I in P'
// only lower part of vector a is changed
// Pa' = a'-u*uT*a'/uu_2

for (size_t i = col_index + 1; i < dim.col_count();++i)
{
    // perform a' <- a' - (u*a/uu_2)u
    row_iterator a(row.index(),row.iterator()+i);
    la::matrix_col_add(a,x,
                       -matrix_col_product<double*,dim_type>(a,x,dim_type())
                       /uu_2,dim_type());
}
return x_norm;
}*/

template <typename input_iterator,typename output_iterator>
void matrix_eigen_decomposition_sym(input_iterator A,
                                    output_iterator V,
                                    output_iterator d,dim<2,2>)
{
    double b = A[1];
    if (b + 1.0 == 1.0)
    {
        d[0] = A[0];
        d[1] = A[3];
        V[0] = 1.0;
        V[1] = 0.0;
        V[2] = 0.0;
        V[3] = 1.0;
        return;
    }
    double a = A[0];
    double b2 = b*b;
    double c = A[3];
    double a_c = a-c;
    double t = std::sqrt(a_c*a_c+4.0*b2);
    d[0] = (a+c+t)*0.5;
    d[1] = d[0]-t;

    double a_l1 = a-d[0];
    double a_l2 = a-d[1];
    double l1 = sqrt(b2+a_l1*a_l1);
    double l2 = sqrt(b2+a_l2*a_l2);
    V[0] = -b/l1;
    V[2] = a_l1/l1;
    V[1] = -b/l2;
    V[3] = a_l2/l2;
}


template<typename input_iterator,typename dim_type>
void matrix_col_swap(input_iterator i1,input_iterator i2,const dim_type& dim)
{
    size_t col_count = dim.col_count();
    // started from 1 for leap iterator problem
    for (size_t i = 1;i < dim.row_count();++i,i1 += col_count, i2 += col_count)
        std::swap(*i1,*i2);
    std::swap(*i1,*i2);
}

template <typename input_iterator,typename output_iterator,typename dym_type>
void matrix_eigenvalue(input_iterator A,output_iterator d,const dym_type& dimension)
{
    typedef typename std::iterator_traits<input_iterator>::value_type value_type;
    typedef value_type* iterator_type;
    const size_t size = dimension.size();
    const size_t dim = dimension.col_count();
    const size_t shift = dim + 1;
    std::vector<value_type> V_(size);
    std::vector<value_type> e_(dim);
    value_type* V = &*V_.begin();
    value_type* e = &*e_.begin();

    std::copy(A,A+size,V);
    std::fill(d,d+dim,value_type(0));

    //void tridiagonalize(void)
    {
        //  This is derived from the Algol procedures tred2 by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
        // Householder reduction to tridiagonal form.
        {
            iterator_type Vrowi = V+size-dim;//n-1 row
            for (size_t i = dim-1;i > 1;--i,Vrowi -= dim)
            {
                value_type h(0),g,f;
                // Generate Householder vector.u
                // x is the lower i-1 row vector of row i
                // h = |x|^2
                for (size_t k = 0; k < i; k++)
                    h += Vrowi[k] * Vrowi[k];
                if (h+value_type(1.0) == value_type(1.0))
                    //if(h < la::eps<value_type>::value)
                {
                    e[i] = Vrowi[i-1];
                    continue;
                }

                f = Vrowi[i-1];
                g = std::sqrt(h);			// g = |x|
                if (f >= value_type(0))
                    g = -g;	// choose sign of g accoring to V[i][i-1]
                e[i] = g;
                h -= f * g;	// h = 1/2|u|^2=1/2(|x-|x|e|^2)=1/2(|x|^2-2*|x|*x0+|x|^2)
                //   = |x|^2-|x|*x0;
                Vrowi[i-1] -= g; // Vrowi x becomes u, u = x-|x|e
                f = value_type(0);
                {
                    iterator_type Vrowj = V;// from the first row
                    for (size_t j = 0; j < i; ++j, Vrowj += dim)
                    {
                        size_t j_1 = j+1;
                        g = value_type(0);
                        iterator_type rowj_1 = V;
                        for (size_t k = 0;k < j_1;++k,rowj_1+=dim)
                            g += Vrowj[k]*Vrowi[k];
                        if (j_1 < i)
                        {
                            iterator_type Vrowk = rowj_1+j; //row j +1 , col j
                            for (size_t k = j_1;k < i;++k,Vrowk += dim)
                                g += (*Vrowk)*Vrowi[k];
                        }
                        e[j] = g/h;
                        f += e[j] * Vrowi[j];
                    }
                }
                d[i] = h;
                {
                    value_type hh = f / (h + h);
                    iterator_type Vrowj = V;// from the first row
                    for (size_t j = 0; j < i; ++j, Vrowj += dim)
                    {
                        f = Vrowi[j];
                        g = e[j]-hh * f;
                        e[j] = g;
                        for (size_t k = 0;k < j+1;++k)
                            Vrowj[k] -= (f * e[k] + g * Vrowi[k]);
                    }
                }
            }
        }

        e[1] = V[dim];
        iterator_type Vdia = V+size-1;
        output_iterator d_iter = d + dim -1;
        d[0] = V[0];
        while (Vdia != V)
        {
            *d_iter = *Vdia;
            --d_iter;
            Vdia -= shift;
        }
    }

    // Symmetric tridiagonal QL algorithm.

    //	void diagonalize(void)
    {
        using namespace std; // for hypot

        //  This is derived from the Algol procedures tql2, by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
        std::copy(e+1,e+dim,e);
        e[dim-1] = value_type(0);

        value_type p,r,b,f(0);
        for (int l = 0,iter = 0;l < dim && iter < 30;++iter)
        {
            size_t m = l;
            // Find small subdiagonal element
            for (;m < dim-1;++m)
            {
                /*
                tst1 = ((d[m+1] > 0) ?
                std::max(tst1,la::abs(d[m]) + d[m+1]):
                std::max(tst1,la::abs(d[m]) - d[m+1]);
                if(tst1+e[m] == tst1)
                break;*/
                if (d[m]+e[m] == d[m])
                    break;

            }
            // If m == l, d[l] is an eigenvalue, go for next
            if ((int)m == l)
            {
                ++l;
                iter = 0;
                continue;
            }

            // Compute implicit shift
            p = (d[l+1]-d[l])/(e[l]*value_type(2));
            r = std::sqrt(1+p*p);
            p = d[m]-d[l]+e[l]/(p+((p<0)? -r:r));

            value_type s(1),c(1),g(0);
            int i = m-1;
            for (; i >= l; i--)
            {
                f = s*e[i];
                b = c*e[i];
                e[i+1] = r = hypot(f,p);
                if (r+f == f && r+p == p)
                {
                    d[i-1] -= g;
                    e[m] = value_type(0);
                    break;
                }
                s = f/r;
                c = p/r;
                p = d[i+1]-g;
                r = (d[i]-p)*s+c*b*value_type(2);
                g = s*r;
                d[i+1] = p+g;
                p = c*r-b;
            } // i loop
            if (r != value_type(0) || i < l)
            {
                e[l] = p;
                e[m] = value_type(0);
                d[l] -= g;
            }
        } // l loop


    }
    std::sort(d,d+dim,std::greater<typename std::iterator_traits<output_iterator>::value_type>());
}

/**
A must be a symmetric matrix
Output V:eigenvectors, *stored in colomn major
Output d:eigenvalues
*/
template <typename input_iterator,typename output_iterator1,typename output_iterator2,typename dym_type>
void matrix_eigen_decomposition_sym(input_iterator A,
                                    output_iterator1 V,
                                    output_iterator2 d,const dym_type& dimension)
{
    typedef typename std::iterator_traits<input_iterator>::value_type value_type;
    const size_t size = dimension.size();
    const size_t dim = dimension.col_count();
    std::vector<value_type> e_(dim+1);
    value_type* e = &*e_.begin();

    std::copy(A,A+size,V);
    std::fill(d,d+dim,value_type(0));

    //void tridiagonalize(void)
    {
        //  This is derived from the Algol procedures tred2 by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
        // Householder reduction to tridiagonal form.
        {
            output_iterator1 Vrowi = V+size-dim;//n-1 row
            for (size_t i = dim-1;i > 1;--i,Vrowi -= dim)
            {
                value_type h(0),g,f;
                // Generate Householder vector.u
                // x is the lower i-1 row vector of row i
                // h = |x|^2
                for (size_t k = 0; k < i; k++)
                    h += Vrowi[k] * Vrowi[k];
                if (h+value_type(1.0) == value_type(1.0))
                    //if(h < la::eps<value_type>::value)
                {
                    e[i] = Vrowi[i-1];
                    continue;
                }

                f = Vrowi[i-1];
                g = std::sqrt(h);			// g = |x|
                if (f >= value_type(0))
                    g = -g;	// choose sign of g accoring to V[i][i-1]
                e[i] = g;
                h -= f * g;	// h = 1/2|u|^2=1/2(|x-|x|e|^2)=1/2(|x|^2-2*|x|*x0+|x|^2)
                //   = |x|^2-|x|*x0;
                Vrowi[i-1] -= g; // Vrowi x becomes u, u = x-|x|e
                f = value_type(0);
                {
                    output_iterator1 Vrowj = V;// from the first row
                    for (size_t j = 0; j < i; ++j, Vrowj += dim)
                    {
                        size_t j_1 = j+1;
                        Vrowj[i] = Vrowi[j]/h;
                        g = value_type(0);
                        output_iterator1 rowj_1 = V;
                        for (size_t k = 0;k < j_1;++k,rowj_1+=dim)
                            g += Vrowj[k]*Vrowi[k];
                        if (j_1 < i)
                        {
                            output_iterator1 Vrowk = rowj_1+j; //row j +1 , col j
                            for (size_t k = j_1;k < i;++k,Vrowk += dim)
                                g += (*Vrowk)*Vrowi[k];
                        }
                        e[j] = g/h;
                        f += e[j] * Vrowi[j];
                    }
                }
                d[i] = h;
                {
                    value_type hh = f / (h + h);
                    output_iterator1 Vrowj = V;// from the first row
                    for (size_t j = 0; j < i; ++j, Vrowj += dim)
                    {
                        f = Vrowi[j];
                        g = e[j]-hh * f;
                        e[j] = g;
                        for (size_t k = 0;k < j+1;++k)
                            Vrowj[k] -= (f * e[k] + g * Vrowi[k]);
                    }
                }
            }
        }

        e[0] = value_type(0);
        d[0] = V[0];
        e[1] = V[dim];
        d[1] = value_type(0);
        V[0] = value_type(1);
        // Accumulate transformations.
        // Also change V from column major to row major
        // Elements in V(j<i,k<i) is row major,
        // Elements in V(j>=i,k>=i) is column major,
        {
            output_iterator1 Vrowi = V;// from the second row
            for (size_t i = 1; i < dim; ++i)
            {
                Vrowi += dim;
                if (d[i] != value_type(0))
                {
                    for (output_iterator1 Vrowj = V;Vrowj != Vrowi;Vrowj += dim)
                    {
                        value_type g = vector_op_dot(Vrowi,Vrowi+i,Vrowj);
                        {
                            output_iterator1 Vrowk = V;
                            for (int k = 0;k < i;++k,Vrowk += dim)
                                Vrowj[k] -= g * Vrowk[i];
                        }
                    }
                }
                d[i] = Vrowi[i];
                Vrowi[i] = value_type(1);
                {
                    output_iterator1 Vrowk = V;
                    for (int k = 0;k < i;++k,Vrowk += dim)
                        Vrowk[i] = Vrowi[k] = value_type(0);

                }
            }
        }


    }

    // Symmetric tridiagonal QL algorithm.

    //	void diagonalize(void)
    {
        using namespace std; // for hypot

        //  This is derived from the Algol procedures tql2, by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.
        ++e;
        value_type p,r,b,f(0);
        for (int l = 0,iter = 0;l < dim && iter < 30;++iter)
        {
            size_t m = l;
            // Find small subdiagonal element
            for (;m < dim-1;++m)
            {
                /*
                tst1 = ((d[m+1] > 0) ?
                std::max(tst1,la::abs(d[m]) + d[m+1]):
                std::max(tst1,la::abs(d[m]) - d[m+1]);
                if(tst1+e[m] == tst1)
                break;*/
                if (d[m]+e[m] == d[m])
                    break;

            }
            // If m == l, d[l] is an eigenvalue, go for next
            if ((int)m == l)
            {
                ++l;
                iter = 0;
                continue;
            }

            // Compute implicit shift
            p = (d[l+1]-d[l])/(e[l]*value_type(2));
            r = std::sqrt(1+p*p);
            p = d[m]-d[l]+e[l]/(p+((p<0)? -r:r));

            value_type s(1),c(1),g(0);
            int i = m-1;
            output_iterator1 Vrowi = V+i*dim;
            output_iterator1 Vrowi_end = Vrowi+dim;
            do
            {
                f = s*e[i];
                b = c*e[i];
                e[i+1] = r = hypot(f,p);
                if (r+f == f && r+p == p)
                {
                    d[i-1] -= g;
                    e[m] = value_type(0);
                    break;
                }
                s = f/r;
                c = p/r;
                p = d[i+1]-g;
                r = (d[i]-p)*s+c*b*value_type(2);
                g = s*r;
                d[i+1] = p+g;
                p = c*r-b;
                // Accumulate transformation.
                vector_op_rot(Vrowi,Vrowi_end,Vrowi_end,c,-s);
                if (--i < l)
                    break;
                Vrowi_end = Vrowi;
                Vrowi -= dim;
            } // i loop
            while (1);

            if (r != value_type(0) || i < l)
            {
                e[l] = p;
                e[m] = value_type(0);
                d[l] -= g;
            }
        } // l loop


        // Sort eigenvalues and corresponding vectors.
        output_iterator1 Vrowi = V;
        for (size_t i = 0; i < dim-1; ++i,Vrowi += dim)
        {
            size_t k = std::max_element(d+i,d+dim)-d;
            if (k != i)
            {
                std::swap(d[k],d[i]);
                vector_op_swap(Vrowi,Vrowi+dim,V+k*dim);
            }
        }
    }
}

/*
Input:
    A: n-by-m matrix
	n <= m
Operation:
    calculate A=U'*S*V
Output:
    A: n singular vector having m dimensions (right side matrix)
	U: n singular vector having n dimensions (left side matrix)
	s: n singular values
*/

template <typename input_iterator,typename output_iterator1,typename output_iterator2,typename dym_type>
void matrix_svd(input_iterator A,output_iterator1 U,output_iterator2 s,dym_type dimension)
{
    typedef typename std::iterator_traits<input_iterator>::value_type value_type;
    int n = dimension.row_count();
    int m = dimension.col_count();
    int nu = n;
    std::vector<value_type> e_(n);
    value_type* e = &*e_.begin();
    value_type* e_end = e+n;
    std::vector<value_type> w_(m);
    value_type* w = &*w_.begin();
    value_type* w_end = w+m;
    value_type max_value = pow(2.0,-52.0);

    // Reduce A to bidiagonal form, storing the diagonal elements
    // in s and the super-diagonal elements in e.

    int nct = std::min(m-1,n);
    int nrt = std::max(0,nu-2);
    {
        int k = 0;
        input_iterator Arowk = A;
        input_iterator Arowk_end = Arowk+m;
        do
        {
            int k_1 = k+1;
            {
                input_iterator Arowk_k = Arowk+k;//A[k][k];

                value_type s_k = vector_op_norm2(Arowk_k,Arowk_end);
                if (s_k == 0.0)
                {
                    s[k] = 0.0;
                    int j = k_1;
                    if (j < n)
                    {
                        input_iterator Arowj_k = Arowk_end+k; // A[k][j];
                        do
                        {
                            e[j] = *Arowj_k;
                            if (++j == n)
                                break;
                            Arowj_k += m;
                        }
                        while (1);
                    }
                }
                else
                {
                    if (Arowk[k] < 0.0)
                        s_k = -s_k;
                    vector_op_scale(Arowk_k,Arowk_end,1.0/s_k);
                    *Arowk_k += 1.0;
                    s[k] = -s_k;
                    int j = k_1;
                    if (j < n)
                    {
                        input_iterator Arowj_k = Arowk_end+k; // A[k][j];
                        do
                        {
                            vector_op_aypx(Arowk_k,Arowk_end,vector_op_dot(Arowk_k,Arowk_end,Arowj_k)/-*Arowk_k,Arowj_k);
                            e[j] = *Arowj_k;
                            if (++j == n)
                                break;
                            Arowj_k += m;
                        }
                        while (1);
                    }
                }
            }
            //now U[k:m][k] stores in A[k:m][k];

            if (k < nrt)
            {
                value_type* e_k_1 = e + k_1;
                value_type e_k_value = vector_op_norm2(e_k_1,e_end);
                if (e_k_value != 0.0)
                {
                    if (*e_k_1 < 0.0)
                        e_k_value = -e_k_value;

                    vector_op_scale(e_k_1,e_end,1.0/e_k_value);
                    *e_k_1 += 1.0;
                    e_k_value = -e_k_value;
                    // Apply the transformation.
                    if (k+1 < m)
                    {
                        value_type* w_k_1 = w+k_1;
                        std::fill(w_k_1,w_end,0.0);

                        //if (j < n)   no need this check
                        {
                            int j = k_1;
                            input_iterator Arowj_k_1 = Arowk_end+k_1;
                            do
                            {
                                vector_op_axpy(w_k_1,w_end,e[j],Arowj_k_1);
                                if (++j == n)
                                    break;
                                Arowj_k_1 += m;
                            }
                            while (1);
                        }

                        //if (j < n)   no need this check
                        {
                            int j = k_1;
                            input_iterator Arowj_k_1 = Arowk_end+k_1;
                            value_type e_k_1_value = *e_k_1;
                            do
                            {
                                vector_op_aypx(w_k_1,w_end,-e[j]/e_k_1_value,Arowj_k_1);
                                if (++j == n)
                                    break;
                                Arowj_k_1 += m;
                            }
                            while (1);
                        }
                    }
                }

                e[k] = e_k_value;
                // store U[k+1:n][k] to A[k][k+1:n]
                std::copy(e+k_1,e_end,U+k*n+k_1);
            }
            max_value=std::max(max_value,std::abs(e[k])+std::abs(s[k]));
            if (++k == nct)
                break;
            Arowk = Arowk_end;
            Arowk_end += m;
        }
        while (1);
    }

    if (m == n)
        s[nct] = *(A+nct*m+nct);
    e[nrt] = A[(n-1)*m + nrt];

    {
        int k = nu-1;
        input_iterator Urowk = U+k*n;
        input_iterator Urowk_end = Urowk+n;
        while (1)
        {
            int k_1 = k+1;
            input_iterator Urowk_k_1 = Urowk+k_1;
            if (k < nrt && (e[k] != 0.0))
            {
                int j = k_1;
                input_iterator Urowj_k_1 = Urowk_end+k_1;
                do
                {
                    vector_op_aypx(Urowk_k_1,Urowk_end,vector_op_dot(Urowk_k_1,Urowk_end,Urowj_k_1)/-(*Urowk_k_1),Urowj_k_1);
                    if (++j == nu)
                        break;
                    Urowj_k_1 += n;
                }
                while (1);
            }
            std::fill(Urowk,Urowk_end,0.0);
            Urowk[k] = 1.0;
            if (--k < 0)
                break;
            Urowk_end = Urowk;
            Urowk -= n;
        }
    }

    // If required, generate U.
    {
        int k = nu-1;
        input_iterator Arowk = A+k*m;
        input_iterator Arowk_end = Arowk+m;
        if (k >= 0)
            while (1)
            {
                input_iterator Arowk_k = Arowk + k;
                if (s[k] != 0.0 && k != nct)
                {
                    int j = k + 1;
                    if (j < nu)
                    {
                        output_iterator1 Arowj_k = Arowk_end+k;
                        while (1)
                        {
                            vector_op_aypx(Arowk_k,Arowk_end,vector_op_dot(Arowk_k,Arowk_end,Arowj_k)/-(*Arowk_k),Arowj_k);
                            if (++j == nu)
                                break;
                            Arowj_k += m;
                        }
                    }
                    vector_op_negate(Arowk_k,Arowk_end);
                    *Arowk_k += 1.0;
                    std::fill(Arowk,Arowk_k,0.0);
                }
                else
                {
                    std::fill(Arowk,Arowk_end,0.0);
                    *Arowk_k = 1.0;
                }
                if (--k < 0)
                    break;
                Arowk_end = Arowk;
                Arowk -= m;
            }
    }
    for (int pp = nu-1;pp > 0;)
    {
        int pp_1 = pp-1;
        if (max_value + e[pp_1] == max_value)
        {
            e[pp_1] = 0.0;
            --pp;
            continue;
        }
        int k;
        for (k = pp-2; k >= 0; k--)
            if (max_value + e[k] == max_value)
            {
                e[k] = 0.0;
                break;
            }
        int ks;
        for (ks = pp; ks > k; ks--)
            if (max_value + s[ks] == max_value)// singularity
            {
                s[ks] = 0.0;
                if (ks == pp)
                {
                    ++k;
                    // Deflate negligible s(p).
                    value_type f(e[pp_1]);
                    e[pp_1] = 0.0;

                    output_iterator1 Urowpp = U+pp*n;
                    output_iterator1 Urowj = Urowpp-n;
                    value_type t,cs,sn;
                    int j = pp_1;
                    while (1)
                    {
                        t = hypot(s[j],f);
                        cs = s[j]/t;
                        sn = f/t;
                        s[j] = t;
                        vector_op_rot(Urowj,Urowj+n,Urowpp,cs,sn);
                        if (j-- != k)
                        {
                            f = -sn*e[j];
                            e[j] *= cs;
                        }
                        if (j < k)
                            break;
                        Urowj -= n;
                    }
                }
                else
                {
                    // Split at negligible s(k).
                    value_type f(e[ks]);
                    e[ks] = 0.0;
                    input_iterator Arowks = A+ks*m;
                    input_iterator Arowj = Arowks;
                    value_type t,cs,sn;
                    int j = ks;
                    while (1)
                    {
                        t = hypot(s[j],f);
                        cs = s[j]/t;
                        sn = f/t;
                        s[j] = t;
                        f = -sn*e[j];
                        e[j] = cs*e[j];
                        vector_op_rot(Arowj,Arowj+m,Arowks,cs,sn);
                        if (++j > pp)
                            break;
                        Arowj += m;
                    }
                }

                break;
            }
        if (ks == k)
        {
            // Perform one qr step.
            ++k;
            // Calculate the shift.
            value_type sp = s[pp];
            value_type sk = s[k];
            value_type ek = e[k];
            value_type b = e[pp_1];
            b *= b;
            b += (s[pp_1] + sp)*(s[pp_1] - sp);
            b /= 2.0;
            value_type c = sp*e[pp_1];
            value_type shift = 0.0;
            c *= c;
            if ((b != 0.0) || (c != 0.0))
            {
                shift = sqrt(b*b + c);
                if (b < 0.0)
                    shift = -shift;
                shift = c/(b + shift);
            }
            value_type f = (sk + sp)*(sk - sp) + shift;
            value_type g = sk*ek;
            value_type t,cs,sn;
            int j = k;
            input_iterator Arowj_1 = A+k*m+m;
            output_iterator1 Urowj_1 = U+k*n+n;
            while (1)
            {
                int j_1 = j+1;
                t = hypot(f,g);
                cs = f/t;
                sn = g/t;

                if (j != k)
                    e[j-1] = t;
                f = cs*s[j] + sn*e[j];
                e[j] = e[j]*cs-sn*s[j];
                g = sn*s[j_1];
                s[j_1] *= cs;

                vector_op_rot(Urowj_1-n,Urowj_1,Urowj_1,cs,sn);

                t = hypot(f,g);
                cs = f/t;
                sn = g/t;

                s[j] = t;
                f = cs*e[j] + sn*s[j_1];
                s[j_1] = s[j_1] * cs - sn*e[j];
                g = sn*e[j_1];
                e[j_1] *= cs;

                vector_op_rot(Arowj_1-m,Arowj_1,Arowj_1,cs,sn);
                if (++j == pp)
                    break;
                Arowj_1 += m;
                Urowj_1 += n;
            }
            e[pp_1] = f;

        }
    } // while

    // make singular value positive
    {
        output_iterator1 Urowk = U;
        output_iterator2 s_iter = s;
        output_iterator2 s_end = s+nu;
        while (1)
        {
            if (*s_iter < 0.0)
            {
                *s_iter = -*s_iter;
                vector_op_negate(Urowk,Urowk+n);
            }
            if (++s_iter == s_end)
                break;
            Urowk += n;
        }
    }
    // sort the values

    {
        input_iterator Arow = A;
        output_iterator1 Urow = U;
        output_iterator2 s_end = s+nu;
        for (output_iterator2 s_to = s+nu-1;s != s_to;++s,Arow += m,Urow += n)
        {
            output_iterator2 s_max= std::max_element(s,s_end);
            if (s_max == s)
                continue;
            std::swap(*s,*s_max);
            size_t dif = s_max-s;
            vector_op_swap(Urow,Urow+n,Urow+dif*n);
            vector_op_swap(Arow,Arow+m,Arow+dif*m);
        }
    }
    if (n != dimension.row_count())
    {
        matrix_inplace_transpose(A,dimension);
    }
}

template <typename input_iterator,typename output_iterator2,typename dym_type>
void matrix_svd(input_iterator A,output_iterator2 s,dym_type dimension)
{
    typedef typename std::iterator_traits<input_iterator>::value_type value_type;
    int n = dimension.row_count();
    int m = dimension.col_count();
    int nu = n;
    std::vector<value_type> e_(n);
    value_type* e = &*e_.begin();
    value_type* e_end = e+n;
    std::vector<value_type> w_(m);
    value_type* w = &*w_.begin();
    value_type* w_end = w+m;
    value_type max_value = pow(2.0,-52.0);

    // Reduce A to bidiagonal form, storing the diagonal elements
    // in s and the super-diagonal elements in e.

    int nct = std::min(m-1,n);
    int nrt = std::max(0,nu-2);
    {
        int k = 0;
        input_iterator Arowk = A;
        input_iterator Arowk_end = Arowk+m;
        do
        {
            int k_1 = k+1;
            {
                input_iterator Arowk_k = Arowk+k;//A[k][k];

                value_type s_k = vector_op_norm2(Arowk_k,Arowk_end);
                if (s_k == 0.0)
                {
                    s[k] = 0.0;
                    int j = k_1;
                    if (j < n)
                    {
                        input_iterator Arowj_k = Arowk_end+k; // A[k][j];
                        do
                        {
                            e[j] = *Arowj_k;
                            if (++j == n)
                                break;
                            Arowj_k += m;
                        }
                        while (1);
                    }
                }
                else
                {
                    if (Arowk[k] < 0.0)
                        s_k = -s_k;
                    vector_op_scale(Arowk_k,Arowk_end,1.0/s_k);
                    *Arowk_k += 1.0;
                    s[k] = -s_k;
                    int j = k_1;
                    if (j < n)
                    {
                        input_iterator Arowj_k = Arowk_end+k; // A[k][j];
                        do
                        {
                            vector_op_aypx(Arowk_k,Arowk_end,vector_op_dot(Arowk_k,Arowk_end,Arowj_k)/-*Arowk_k,Arowj_k);
                            e[j] = *Arowj_k;
                            if (++j == n)
                                break;
                            Arowj_k += m;
                        }
                        while (1);
                    }
                }
            }
            //now U[k:m][k] stores in A[k:m][k];

            if (k < nrt)
            {
                value_type* e_k_1 = e + k_1;
                value_type e_k_value = vector_op_norm2(e_k_1,e_end);
                if (e_k_value != 0.0)
                {
                    if (*e_k_1 < 0.0)
                        e_k_value = -e_k_value;

                    vector_op_scale(e_k_1,e_end,1.0/e_k_value);
                    *e_k_1 += 1.0;
                    e_k_value = -e_k_value;
                    // Apply the transformation.
                    if (k+1 < m)
                    {
                        value_type* w_k_1 = w+k_1;
                        std::fill(w_k_1,w_end,0.0);

                        //if (j < n)   no need this check
                        {
                            int j = k_1;
                            input_iterator Arowj_k_1 = Arowk_end+k_1;
                            do
                            {
                                vector_op_axpy(w_k_1,w_end,e[j],Arowj_k_1);
                                if (++j == n)
                                    break;
                                Arowj_k_1 += m;
                            }
                            while (1);
                        }

                        //if (j < n)   no need this check
                        {
                            int j = k_1;
                            input_iterator Arowj_k_1 = Arowk_end+k_1;
                            value_type e_k_1_value = *e_k_1;
                            do
                            {
                                vector_op_aypx(w_k_1,w_end,-e[j]/e_k_1_value,Arowj_k_1);
                                if (++j == n)
                                    break;
                                Arowj_k_1 += m;
                            }
                            while (1);
                        }
                    }
                }

                e[k] = e_k_value;
            }
            max_value=std::max(max_value,std::abs(e[k])+std::abs(s[k]));
            if (++k == nct)
                break;
            Arowk = Arowk_end;
            Arowk_end += m;
        }
        while (1);
    }

    if (m == n)
        s[nct] = *(A+nct*m+nct);
    e[nrt] = A[(n-1)*m + nrt];
    for (int pp = nu-1;pp > 0;)
    {
        int pp_1 = pp-1;
        if (max_value + e[pp_1] == max_value)
        {
            e[pp_1] = 0.0;
            --pp;
            continue;
        }
        int k;
        for (k = pp-2; k >= 0; k--)
            if (max_value + e[k] == max_value)
            {
                e[k] = 0.0;
                break;
            }
        int ks;
        for (ks = pp; ks > k; ks--)
            if (max_value + s[ks] == max_value)// singularity
            {
                if (ks == pp)
                {
                    ++k;
                    // Deflate negligible s(p).
                    value_type f(e[pp_1]);
                    e[pp_1] = 0.0;
                    value_type t,cs,sn;
                    for (int j = pp_1;j >= k;)
                    {
                        t = hypot(s[j],f);
                        cs = s[j]/t;
                        sn = f/t;
                        s[j] = t;
                        if (j-- != k)
                        {
                            f = -sn*e[j];
                            e[j] *= cs;
                        }
                    }
                }
                else
                {
                    // Split at negligible s(k).
                    value_type f(e[ks]);
                    e[ks] = 0.0;
                    value_type t,cs,sn;
                    for (int j = ks;j <= pp;++j)
                    {
                        t = hypot(s[j],f);
                        cs = s[j]/t;
                        sn = f/t;
                        s[j] = t;
                        f = -sn*e[j];
                        e[j] = cs*e[j];
                    }
                }

                break;
            }
        if (ks == k)
        {
            // Perform one qr step.
            ++k;
            // Calculate the shift.
            value_type sp = s[pp];
            value_type sk = s[k];
            value_type ek = e[k];
            value_type b = e[pp_1];
            b *= b;
            b += (s[pp_1] + sp)*(s[pp_1] - sp);
            b /= 2.0;
            value_type c = sp*e[pp_1];
            value_type shift = 0.0;
            c *= c;
            if ((b != 0.0) || (c != 0.0))
            {
                shift = sqrt(b*b + c);
                if (b < 0.0)
                    shift = -shift;
                shift = c/(b + shift);
            }
            value_type f = (sk + sp)*(sk - sp) + shift;
            value_type g = sk*ek;
            value_type t,cs,sn;
            for (int j = k;j < pp;++j)
            {
                int j_1 = j+1;
                t = hypot(f,g);
                cs = f/t;
                sn = g/t;

                if (j != k)
                    e[j-1] = t;
                f = cs*s[j] + sn*e[j];
                e[j] = e[j]*cs-sn*s[j];
                g = sn*s[j_1];
                s[j_1] *= cs;
                t = hypot(f,g);
                cs = f/t;
                sn = g/t;

                s[j] = t;
                f = cs*e[j] + sn*s[j_1];
                s[j_1] = s[j_1] * cs - sn*e[j];
                g = sn*e[j_1];
                e[j_1] *= cs;
            }
            e[pp_1] = f;

        }
    } // while
    for (size_t i = 0;i < nu;++i)
        if (s[i] < 0.0)
            s[i] = -s[i];
    std::sort(s,s+nu,std::greater<value_type>());
}






}//namespace math

#endif//MATRIX_OP_HPP
