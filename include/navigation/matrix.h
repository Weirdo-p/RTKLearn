/**************************************************************
 * @author XUZHUO WHU
 * 由于模板类声明和定义必须在一个文件中
 * 所以该文件包含Matrix的声明和定义
 * 2020 09 07
**************************************************************/

#include <iostream>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <omp.h>

using namespace std;

#ifndef _MATRIX_H_
#define _MATRIX_H_

const int Dynamic = -1;


/*******************************************************
 * try to complete matrix algorithm with style of Eigen
 * template parameters:
 * param Type  data type you want to use
 * param _0  rows of your matrix
 * param _1  cols of your matrix
*******************************************************/
template <class Type, int _0, int _1>
class Matrix
{
public:
    Matrix();

    Matrix(Matrix<Type, _0, _1> &matrix);

    Matrix(Matrix<Type, _0, _1> &&matrix);

    /* this function leads to deadly error */
    /* still working on it */
    /* 从监视上看似乎没有析构函数  程序也会释放空间 */

    /****************************************
     *  2020 09 11
     * 局部变量return以后与返回后的值共用一段地址
     * 故不需要释放局部变量
    ****************************************/

    /**************************************
     * 2020 09 22
     * 内存泄露的缘由在于局部变量未释放！！！
     * 使用C++智能指针可解决
     * 目前采用手动释放的方法
    **************************************/
    ~Matrix(){
        if (data)
            delete [] this->data;
        data = nullptr;
        rows = cols = 0;
    }

/* initialize method */
public:
    Matrix<Type, _0, _1> Identity();

    Matrix<Type, _0, _1> Zero();

/* basic functions */
public:
    int row() const;

    int col() const;

    Matrix<Type, _1, _0> transpose();

    Matrix<Type, _0, _1> inverse(int &flag);

    void SetElement(const int row, const int col, const Type num) const;

    Type* _2array();

/* override some operators */
public:
    /*********************
     * 返回引用可以便捷的赋值
     * 2020 10 10
    *********************/
    Type & operator() (const int row, const int col) const;

    Matrix(int row, int col);

    Matrix<Type, _0, _1> operator+(const Matrix<Type, _0, _1> &matrix);

    template<int srow, int scol>
    Matrix<Type, Dynamic, Dynamic> block(int nrow, int ncol);

    template<class Type_, class Param_Type, int _0_, int _1_>
    friend Matrix<Type_, _0_, _1_> operator+(Param_Type num, const Matrix<Type_, _0_, _1_> &matrix);

    template<class Type_, class Param_Type, int _0_, int _1_>
    friend Matrix<Type_, _0_, _1_> operator+(const Matrix<Type_, _0_, _1_> &matrix, Param_Type num);

    Matrix<Type, _0, _1> operator-(const Matrix<Type, _0, _1> &matrix);

    template<class Type_, class Param_Type, int _0_, int _1_>
    friend Matrix<Type_, _0_, _1_> operator-(Param_Type num, const Matrix<Type_, _0_, _1_> &matrix);

    template<class Type_, class Param_Type, int _0_, int _1_>
    friend Matrix<Type_, _0_, _1_> operator-(const Matrix<Type_, _0_, _1_> &matrix, Param_Type num);

    template<int _2>
    Matrix<Type, _0, _2> operator*(const Matrix<Type, _1, _2> &matrix);

    Matrix<Type, _0, _1> operator*(const int a);

    Matrix<Type, _0, _1> operator*(const double a);
    /* print the matrix */
    template<class Type1, int _0_, int _1_>
    friend ostream & operator<<(ostream &out, Matrix<Type1, _0_, _1_> &matrix);

    template<class Type1, int _0_, int _1_>
    friend istream & operator>>(istream &in, Matrix<Type1, _0_, _1_> &matrix);

    const Matrix<Type, _0, _1> & operator= (const Matrix<Type, _0, _1> matrix);



    void deleteMatrix();

    void resize(int row, int col);

    void SetSize(int rows, int cols);

    double norm();

    /* try to override = but failed */
    // const Type operator=(const int num);
public:
    template<class Type1, int _0_, int _1_>
    void block(int start_row, int start_col, Matrix<Type1, _0_, _1_> matrix);

private:
    /* rows of the matrix */
    int rows = -2;

    /* cols of the matrix */
    int cols = -2;

    /* data of the matrix */
    Type* data = nullptr;
};

/* constructor */
/* when this function called */
/* matrix will be initialized into null matrix(all elements are set 0) */
/* 2020 09 09 user should initialize the matrix manually */
template<class Type, int _0, int _1>
Matrix<Type, _0, _1>::Matrix() {
    if(_0 == Dynamic || _1 == Dynamic) {
        this->rows = _0;
        this->cols = _1;
        data = nullptr;
    }
    else {
        this->rows = _0;
        this->cols = _1;
        this->data = new Type[rows * cols];
    }
}
template<class Type, int _0, int _1>
Matrix<Type, _0, _1>::Matrix(Matrix<Type, _0, _1>& matrix) {
    if (matrix.data) {
        if (data)
            delete[] data;
        data = new Type [matrix.rows * matrix.cols];
        memcpy(data, matrix.data, sizeof(Type) * matrix.rows * matrix.cols);
        rows = matrix.rows; cols = matrix.cols;
    } else {
        data = nullptr;
        rows = matrix.rows; cols = matrix.cols;
    }
}
template<class Type, int _0, int _1>
Matrix<Type, _0, _1>::Matrix(Matrix<Type, _0, _1> &&matrix) {
    // move
    if (matrix.data) {
        data = matrix.data;
        rows = matrix.rows; cols = matrix.cols;
        matrix.data = nullptr;
    }
}

template<class Type_, int _0, int _1>
Matrix<Type_, _0, _1>::Matrix(int row, int col) {
    if(row < 0 || col < 0) {
        cerr << "dim error" << endl;
        exit(-1);
    }
    this->rows = row;
    this->cols = col;
    this->data = new Type_[row * col];
}

template<class Type, int _0, int _1>
void Matrix<Type, _0, _1>::resize(int row, int col)
{
    if(this->data != nullptr)
        delete[] this->data;
    this->rows = row;
    this->cols = col;
    data = new Type[row * col];
}


/* return row of this matrix */
template<class Type, int _0, int _1>
int Matrix<Type, _0, _1>::row() const
{
    return this->rows;
}

/* return col of this matrix */
template<class Type, int _0, int _1>
int Matrix<Type, _0, _1>::col() const
{
    return this->cols;
}

/* override () to get element in matrix */
template<class Type, int _0, int _1>
Type & Matrix<Type, _0, _1>::operator()(const int row, const int col) const
{
    if(row >= this->rows || col >= this->cols)
    {
        cerr << "assertion occurred! fail to find the element, please check" << endl;
    }

    // to get the position of the element
    int position = row * this->cols + col;

    // vector
    if(this->cols == 1)
        position = row;
    if(this->rows == 1)
        position = col;

    return (this->data[position]);
}

/* override + to plus two matrix */
template<class Type, int _0, int _1>
Matrix<Type, _0, _1> Matrix<Type, _0, _1>::operator+(const Matrix<Type, _0, _1> &matrix)
{
    if(matrix.row() == 0 || matrix.col() == 0)
    {
        cerr << "can not plus two empty matrix!" << endl;
        exit(-1);
    }
    if(matrix.row() != this->row() || 
       matrix.col() != this->col())
    {
        cerr << "can not plus matrix with different rows/cols" << endl;
        exit(-1);
    }

    Matrix<Type, _0, _1> result;
    if(_0 == Dynamic || _1 == Dynamic) 
        result.resize(matrix.row(), matrix.col());    

    for(int i = 0; i < this->row(); ++i)
        for(int j = 0; j < this->col(); ++j)
        {
            Type sum = this->operator()(i, j) + matrix(i, j);
            result.SetElement(i, j, sum);
        }
    return result;
}

// Param_Type + Matrix
template<class Type_, class Param_Type, int _0_, int _1_>
Matrix<Type_, _0_, _1_> operator+(Param_Type num, const Matrix<Type_, _0_, _1_> &matrix)
{
    Matrix<Type_, _0_, _1_> result;
    for(int i = 0; i < matrix.row(); ++i)
        for(int j = 0; j < matrix.col(); ++j)
            result.SetElement(i, j, matrix(i, j) + Type_(num));
    return result;
}

// Matrix + Type
template<class Type_, class Param_Type, int _0_, int _1_>
Matrix<Type_, _0_, _1_> operator+(const Matrix<Type_, _0_, _1_> &matrix, Param_Type num)
{
    Matrix<Type_, _0_, _1_> result;
    for(int i = 0; i < matrix.row(); ++i)
        for(int j = 0; j < matrix.col(); ++j)
            result.SetElement(i, j, matrix(i, j) + Type_(num));
    return result;
}

// Matrix - Type
template<class Type_, class Param_Type, int _0_, int _1_>
Matrix<Type_, _0_, _1_> operator-(const Matrix<Type_, _0_, _1_> &matrix, Param_Type num)
{
    Matrix<Type_, _0_, _1_> result;
    for(int i = 0; i < matrix.row(); ++i)
        for(int j = 0; j < matrix.col(); ++j)
            result.SetElement(i, j, matrix(i, j) - Type_(num));
    return result;
}

// Param_Type - Matrix
template<class Type_, class Param_Type, int _0_, int _1_>
Matrix<Type_, _0_, _1_> operator-(Param_Type num, const Matrix<Type_, _0_, _1_> &matrix)
{
    Matrix<Type_, _0_, _1_> result;
    for(int i = 0; i < matrix.row(); ++i)
        for(int j = 0; j < matrix.col(); ++j)
            result.SetElement(i, j, Type_(num) - matrix(i, j));
    return result;
}

/* override - to minus two matrix */
template<class Type, int _0, int _1>
Matrix<Type, _0, _1> Matrix<Type, _0, _1>::operator-(const Matrix<Type, _0, _1> &matrix)
{
    if(matrix.row() == 0 || matrix.col() == 0)
    {
        cerr << "can not plus two empty matrix!" << endl;
        exit(-1);
    }
    if(matrix.row() != this->row() || 
       matrix.col() != this->col())
    {
        cerr << "can not plus matrix with different rows/cols" << endl;
        exit(-1);
    }

    Matrix<Type, _0, _1> result;
    if(_0 == Dynamic || _1 == Dynamic) {
        result.resize(matrix.row(), matrix.col());
    }

    for(int i = 0; i < this->row(); ++i) {
        #pragma omp parallel for
        for(int j = 0; j < this->col(); ++j)
        {
            Type sum = this->operator()(i, j) - matrix(i, j);
            result.SetElement(i, j, sum);
        }
    }
    return result;
}

/* 赋值 */
/* 可以使用括号赋值 */
/* 此方法已弃用，但仍然兼容此方法 */
/* 2020 10 10 */
template<class Type, int _0, int _1>
void Matrix<Type, _0, _1>::SetElement(const int row, const int col, const Type num) const
{
    if(row < 0 || row >= this->row() ||
       col < 0 || col >= this->col())
    {
        cerr << "Assertion occurred! Please check!" << endl;
        exit(-1);
    }
    int position = row * this->cols + col;
    this->data[position] = num;
}

/* return the transpose of your matrix */
template<class Type, int _0, int _1>
Matrix<Type, _1, _0> Matrix<Type, _0, _1>::transpose()
{
    Matrix<Type, _1, _0> result;
    if(_1 == Dynamic || _0 == Dynamic)
        result.resize(this->cols, this->rows);
    #pragma omp parallel for
    for(int i = 0; i < result.row(); ++i)
    {
        #pragma omp parallel for
        for(int j = 0; j < result.col(); ++j)
        {
            Type num = this->operator()(j, i);
            result.SetElement(i, j, num);
        }
    }
    return result;
}

/* override * to multiply two matrix */
/* be attention to their rows and cols */
template<class Type, int _0, int _1>
template<int _2>
Matrix<Type, _0, _2> Matrix<Type, _0, _1>::operator*(const Matrix<Type, _1, _2> &matrix)
{
    if(matrix.row() != this->col())
    {
        cerr << "you can not mutiply two matrix which the col of the first does not match the row of the second" << endl;
        exit(-1);
    }
    Matrix<Type, _0, _2> result_;
    if(result_.row() == Dynamic)
        result_.resize(this->row(), matrix.col());

    #pragma omp parallel for
    for(int i = 0; i < this->row(); ++i)
    {
        #pragma omp parallel for
        for (int k = 0; k < matrix.col(); ++k)
        {
            Type sum = 0;
            for (int j = 0; j < this->col(); ++j)
            {
                sum = sum + this->operator()(i, j) * matrix(j, k);
            }
            result_.SetElement(i, k, sum);
        }
    }
    return result_;
}
template<class Type, int _0, int _1>
void Matrix<Type, _0, _1>::SetSize(int rows, int cols)
{
    this->rows = rows;
    this->cols = cols;
}


template<class Type, int _0, int _1>
const Matrix<Type, _0, _1> & Matrix<Type, _0, _1>::operator= (const Matrix<Type, _0, _1> matrix)
{
    if(this == &matrix)
    {
        return *this;
    }
    if(this->data != nullptr)
    {
        delete [] data;
        data = nullptr;
    }

    this->cols = matrix.col();
    this->rows = matrix.row();
    this->data = new Type[matrix.col() * matrix.row()];
    #pragma omp parallel for
    for(int i = 0; i < matrix.row(); ++i) {
        #pragma omp parallel for
       for(int j = 0; j < matrix.col(); ++j)
       {
            double tmp = matrix(i, j);
            this->operator()(i, j) = tmp;
       }
    }
    
    return (*this);
}

/********************************
 * 当模板类中有友元函数时
 * 友元函数的模板参数必须与类模板不同
 * 才会被认为为友元函数
 * 否则编译不通过
********************************/
template<class Type1, int _0_, int _1_>
ostream & operator<<(ostream &out, Matrix<Type1, _0_, _1_> &matrix)
{
    for(int i = 0; i < matrix.row(); ++i) {
        // format the out
        for(int j = 0; j < matrix.col(); ++j) {
            out << setw(15) << matrix.operator()(i, j) << " ";
        }
        out << endl;
    }
    return out;
}

template<class Type1, int _0_, int _1_>
istream & operator>>(istream &in, Matrix<Type1, _0_, _1_> &matrix)
{
    for(int i = 0; i < matrix.row(); ++i)
        for(int j = 0; j < matrix.col(); ++j)
        {
            Type1 temp;
            in >> temp;
            matrix.SetElement(i, j, temp);
        }
    if(!in)
    {
        cerr << "error occurred when input" << endl;
    }

    return in;
}

/* to initialize matrxi as identity */
template<class Type, int _0, int _1>
Matrix<Type, _0, _1> Matrix<Type, _0, _1>::Identity()
{
    // 这样的实现方法会导致过多占用内存
    // 故不使用静态成员的方法实现
    // 2020 09 10
    // Matrix<Type, _0, _1> Identity;

    // 2020 09 11 
    // 恢复使用静态方法  问题已解决

    if(this->row() != this->col())
    {
        cerr << "identity matrix can only be initialized by square matrix" << endl;
        return this->Zero();
    }
    #pragma omp parallel for
    for (int i = 0; i < this->row(); ++i) {
        #pragma omp parallel for
        for(int j = 0; j < this->col(); ++j)
            if(i == j)
                this->SetElement(i, j, 1);
            else
                this->SetElement(i, j, 0);
    }
    return (*this);
}

/* to initialize matrix as zero matrix */
template<class Type, int _0, int _1>
Matrix<Type, _0, _1> Matrix<Type, _0, _1>::Zero()
{
    // 这样的实现方法会导致过多占用内存
    // 故不使用静态成员的方法实现
    // 2020 09 10
    // Matrix<Type, _0, _1> Zero;

    // 2020 09 11 
    // 恢复使用静态方法  问题已解决
    // Matrix<Type, _0, _1> Zero_;
    #pragma omp parallel for
    for (int i = 0; i < this->row(); ++i) {
        #pragma omp parallel for
        for(int j = 0; j < this->col(); ++j) {
            this->SetElement(i, j, 0);
        }
    }
    return (*this);
}

/**************************************************************************************
 * caculate the inverse of this matrix
 * reference https://blog.csdn.net/qq_34122194/article/details/78055057
 * 
 * param [out] Identity_ the inverse matrix(should be initialized as identity matrix)
 * return bool true stands for successfully inverse 
 *************************************************************************************/
template<class Type, int _0, int _1>
Matrix<Type, _0, _1> Matrix<Type, _0, _1>::inverse(int &flag)
{
    // 此单位阵为待求逆矩阵的增广矩阵
    Matrix<Type, _0, _1> Identity_;
    Matrix<Type, _0, _1> copy;
    copy = (*this);

    if(_0 == Dynamic || _1 == Dynamic)
        Identity_.resize(copy.row(), copy.col());
    Identity_.Identity();

    if(this->col() != this->row())
    {
        cerr << "only square matrix can call this function(inverse)" << endl;
        flag = 0;
        copy.deleteMatrix();
        return Identity_;
    }
    // cout << copy << endl << endl;

    // 将主元置于对角线上
    for (int i = 0; i < copy.row() - 1; ++i)
    {
        // 目前主元位置
        // i行i列
        int posi = i;
        // 初始化主元绝对值大小
        Type value = abs(copy.operator()(i, i));

        // 对第 i+1行开始寻找真正主元
        #pragma omp parallel for
        for(int j = i + 1; j < copy.row(); ++j)
        {
            Type temp = abs(copy.operator()(j, i));
            if(temp > value)
            {
                posi = j;
                value = temp;
            }
        }
        // 若主元为0  则矩阵不可逆
        if (abs(value) <= 1e-10)
        {
            ofstream debug ("./inverse.txt", ios::app);
            cerr << "this matrix is not invertible" << endl << endl;
            debug << *this << endl;
            debug.close();
            flag = 0;
            copy.deleteMatrix();
            return Identity_;
        }
        // 如果主元不在第i行 则交换两行
        // 同时也要将增广矩阵（单位矩阵）交换相应的两行
        if(posi != i)
        {
            #pragma omp parallel for
            for(int j = 0; j < copy.col(); ++j)
            {
                // 先交换矩阵
                Type temp = copy.operator()(i, j);
                copy.SetElement(i, j, copy.operator()(posi, j));
                copy.SetElement(posi, j, temp);

                // 交换单位阵
                temp = Identity_(i, j);
                Identity_.SetElement(i, j, Identity_(posi, j));
                Identity_.SetElement(posi, j, temp);
            }
        }

        // 将下三角矩阵变为0 (应该可以减少一个循环)
        // 经过实践 并不能少循环
        #pragma omp parallel for
        for(int j = i + 1; j < copy.row(); ++j)
        {
            Type factor = copy.operator()(j, i) / copy.operator()(i, i);
            #pragma omp parallel for
            for(int k = 0; k < copy.col(); ++k)
            {
                Type temp = copy.operator()(j, k) - copy.operator()(i, k) * factor;
                copy.SetElement(j, k, temp);

                temp = Identity_(j, k) - Identity_(i, k) * factor;
                Identity_.SetElement(j, k, temp);
            }
        }
    }

    // 判断对角元素是否为0
    // 若为0则矩阵不可逆
    // 若不为0 则将对角置为1 
    // 且将上三角置为0
    for(int i = copy.row() - 1; i > -1; --i)
    {
        Type temp = copy.operator()(i, i);
        if(abs(temp) <= 1e-15)
        {
            cerr << "this matrix is not invertible" << endl;
            cout << *this << endl;
            flag = 0;
            copy.deleteMatrix();
            return Identity_;
        }
        // 将对角置1
        #pragma omp parallel for
        for(int j = 0; j < copy.col(); ++j)
        {
            Type deno = copy.operator()(i, j) / temp;
            copy.SetElement(i, j, deno);

            deno = Identity_(i, j) / temp;
            Identity_.SetElement(i, j, deno);
        }
        // 同样的方法上三角置0
        #pragma omp parallel for
        for(int j = 0; j < i; ++j)
        {
            Type deno = copy.operator()(j, i);
            #pragma omp parallel for
            for(int k = 0; k < copy.col(); ++k)
            {
                Type temp = copy.operator()(j, k) - copy.operator()(i, k) * deno;
                copy.SetElement(j, k, temp);
                temp = Identity_(j, k) - Identity_(i, k) * deno;
                Identity_.SetElement(j, k, temp);
            }
        }
    }
    flag = 1;
    copy.deleteMatrix();
    return Identity_;
}

template<class T, int _0, int _1>
void Matrix<T, _0, _1>::deleteMatrix()
{
    if(data == nullptr)
        return ;
    delete[] this->data;
    data = nullptr;
}

template<class T, int _0, int _1>
double Matrix<T, _0, _1>::norm()
{
    if(this->rows <= 0 || this->cols <= 0)
        return -1;
    double sum = 0;
    for(int i = 0; i < rows; ++i)
        for(int j = 0; j < cols; ++j)
            sum += this->operator()(i, j) * this->operator()(i, j);
    return sqrt(sum);
}

template <class Type, int _0, int _1>
template <class Type1, int _0_, int _1_>
void Matrix<Type, _0, _1>::block(int s_row, int s_col, Matrix<Type1, _0_, _1_> matrix) {
    if (matrix.col() <=0 || matrix.row() <= 0) return;
    if (s_row + matrix.rows > rows || s_col + matrix.cols > cols) return;
    #pragma omp parallel for
    for(int i = 0; i < matrix.rows; ++i) {
        #pragma omp parallel for
        for(int j = 0; j < matrix.cols; ++j)
            this->operator()(i + s_row, j + s_col) = matrix(i, j);
    }
}

template <class Type, int _0, int _1>
Matrix<Type, _0, _1> Matrix<Type, _0, _1>::operator*(const int a) {
    Matrix<Type, _0, _1> temp = *this;
    for (int i = 0; i < rows; ++i) {
        #pragma omp parallel for
        for(int j = 0; j < cols; ++j)
            temp(i, j) = this->operator()(i, j) * a;
    }
    return temp;
}

template <class Type, int _0, int _1>
Matrix<Type, _0, _1> Matrix<Type, _0, _1>::operator*(const double a) {
    Matrix<Type, _0, _1> temp = *this;
    for (int i = 0; i < rows; ++i) {
        #pragma omp parallel for
        for(int j = 0; j < cols; ++j)
            temp(i, j) = this->operator()(i, j) * a;
    }
    return temp;
}

template <class Type, int _0, int _1>
Type* Matrix<Type, _0, _1>::_2array() {
    if (!data) return nullptr;
    return data;
}

template <class Type, int _0, int _1>
template<int srow, int scol>
Matrix<Type, Dynamic, Dynamic> Matrix<Type, _0, _1>::block(int nrow, int ncol) {
    if (srow + nrow > rows || scol + ncol > cols) {
        cout << "error happened in Matrix::block" << endl;
        exit(-1);
    }
    Matrix<Type, Dynamic, Dynamic> result(nrow, ncol);
    for (int i_row = 0; i_row < nrow; ++ i_row) 
        for (int i_col = 0; i_col < ncol; ++ i_col) 
            result(i_row, i_col) = this->operator()(i_row + srow, i_col + scol);
    return result;
}

#endif

