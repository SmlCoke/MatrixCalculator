#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>

// Matrix 模板类声明

template<typename T>
class Matrix {
private:
    std::vector<std::vector<T>> data;
    int rows, cols;
public:
    Matrix(int r, int c);          // 构造函数：创建 r 行 c 列的矩阵，元素默认初始化为 T()
    int row_count() const;         // 获取行列数
    int col_count() const;
    T& at(int i, int j);           // 访问元素
    const T& at(int i, int j) const;   // 访问元素
    Matrix<T> operator+(const Matrix<T>& other) const;  // 重载加法
    Matrix<T> operator-(const Matrix<T>& other) const;  // 重载减法
    Matrix<T> operator-() const;                       // 重载负号
    void display() const;                              // 打印矩阵
    Matrix<T> operator*(const Matrix<T>& other) const;   // 重载乘法
    void swap_line(int i, int j);                     // 对换变换
    Matrix<T> to_hermite() const;                     // 化为Hermite Normal Form
    Matrix<T> to_hermite_mini() const;                // 化为Hermite Normal Form（简易版：不进行主元选取）
    T get_det() const;                                // 求行列式
    Matrix<T> to_inv() const;                         // 求逆 
};

// ================== 实现部分 ==================

template<typename T>
Matrix<T>::Matrix(int r, int c) : rows(r), cols(c), data(r, std::vector<T>(c, T())) {}
// vector构造函数： vector(len, type)
// T()：表示默认构造函数


template<typename T>
int Matrix<T>::row_count() const { return rows; }

template<typename T>
int Matrix<T>::col_count() const { return cols; }

template<typename T>
T& Matrix<T>::at(int i, int j) {
    if (i < 0 || i >= rows || j < 0 || j >= cols)
        throw std::out_of_range("Matrix indices out of range");
    return data[i][j];
}
 // 返回的是引用，可以通过此函数修改矩阵元素

template<typename T>
const T& Matrix<T>::at(int i, int j) const {
    if (i < 0 || i >= rows || j < 0 || j >= cols)
        throw std::out_of_range("Matrix indices out of range");
    return data[i][j];
}
// 返回的是常量引用，不能通过此函数修改矩阵元素
// 函数末尾const修饰，表示这是一个常成员函数，函数体内部不允许修改数据成员变量

template<typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& other) const {
    if (rows != other.rows || cols != other.cols)
        throw std::invalid_argument("Matrix dimensions must match for addition");
    Matrix<T> result(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result.at(i, j) = this->at(i, j) + other.at(i, j);
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& other) const {
    if (rows != other.rows || cols != other.cols)
        throw std::invalid_argument("Matrix dimensions must match for addition");
    Matrix<T> result(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result.at(i, j) = this->at(i, j) - other.at(i, j);
    return result;
}

template<typename T>
Matrix<T> Matrix<T>::operator-() const {
    Matrix<T> result(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result.at(i, j) = -this->at(i, j);
    return result;
}

template<typename T>
void Matrix<T>::display() const {
    for (const auto& row : data) {
        std::cout << "| ";
        for (const auto& elem : row) {
            std::cout << elem << " ";
        }
        std::cout << "|" << std::endl;
    }
    std::cout << std::endl;
}

template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const {
    if (this->cols != other.rows) {
        throw std::invalid_argument("The left multiplier and right multiplier rows and columns do not match");
    } else {
        Matrix<T> results(this->rows, other.cols);
        for (int i = 0; i < this->rows; i++) {
            for (int j = 0 ; j < other.cols; j++) {
                results.at(i, j) = T();
                for (int k = 0 ; k < this-> cols; k++) {
                    results.at(i, j) = results.at(i, j) + this->data[i][k]*other.data[k][j];
                }
            }
        }
        return results;
    }
}

template<typename T>
void Matrix<T>::swap_line(int i, int j) {
    std::vector<T> temp(cols, T());
    temp = data[i];
    data[i] = data[j];
    data[j] = temp;
}

template<typename T>
Matrix<T> Matrix<T>::to_hermite() const {
    Matrix<T> res = *this;   // 独立出一个新矩阵
    std::vector<int> pivot_index(cols,0);
    for (int i = 0 ; i < cols; i++)
    {
        pivot_index[i] = i;          // 理想主元列标，当不存在线性相关的向量时，该数组不会更新
    }
    for (int j = 0 ; j < cols; j ++) 
    {
        int i = pivot_index[j];  // 标记最佳主元
        int flag = pivot_index[j];    // 标记首个非零元，高斯元从此行开始
        if (i >= rows-1)  break;       // 越界，证明cols > rows
        while(i < rows && res.data[i][j] == T()) {   // 寻找首个非零元
            i++;
            flag++;
        }
        if (i == rows)    // 该列没有非零元
        {
            for (int k = j; k < cols; k ++)  pivot_index[k]--;  // 如果没有非零元，那么下一次高斯消元时应该从当前理应操作的行开始
        }
        else              // 该列有非零元
        {
            for (int k = i + 1; k < rows; k++)   // 为了避免数值误差/分数膨胀，选取的主元应该是该列绝对值最大的非零元
            {
                if (abs(res.at(k,j)) > abs(res.at(i,j)))
                {
                    i = k;
                }
            }
            res.swap_line(i, pivot_index[j]);   // 本来应该是j，但是由于存在全零列，所以应该是pivot_index[j]
            for (int k = flag + 1; k < rows; k++)
            {
                T factor = res.at(k,j)/res.at(pivot_index[j], j);   // 倍加因子
                for (int l = j; l < cols; l++)          // 高斯消元
                {
                    res.at(k,l) = res.at(k,l) - res.at(pivot_index[j], l)*factor;
                }
            }
        }
    }
    return res;
}

template<typename T>
Matrix<T> Matrix<T>::to_hermite_mini() const {
    Matrix<T> res = *this;   // 独立出一个新矩阵
    std::vector<int> pivot_index(cols,0);
    for (int i = 0 ; i < cols; i++)
    {
        pivot_index[i] = i;
    }
    for (int j = 0 ; j < cols; j ++) 
    {
        int i = pivot_index[j];  // 标记首个非零元
        if (i >= rows-1)  break;       // 越界，证明cols > rows
        while(i < rows && res.data[i][j] == T()) {
            i++;
        }
        if (i == rows)    // 该列没有非零元
        {
            for (int k = j; k < cols; k ++)  pivot_index[k]--;
        }
        else              // 该列有非零元
        {
            res.swap_line(i, pivot_index[j]);   // 本来应该是j，但是由于存在全零列，所以应该是pivot_index[j]
            for (int k = i + 1; k < rows; k++)
            {
                T factor = res.at(k,j)/res.at(pivot_index[j], j);    // 倍加因子
                for (int l = j; l < cols; l++)   // 高斯消元
                { 
                    res.at(k,l) = res.at(k,l) - res.at(pivot_index[j], l)*factor;
                }
            }
        }
    }
    return res;
}

template<typename T>
T Matrix<T>::get_det() const {
    if (rows != cols)
    {
        throw std::invalid_argument("Only square matrices can be used to find the determinant");
    }
    else
    {
        Matrix<T> res = *this;   // 独立出一个新矩阵
        int swap_times = 0;      // 纪录对换变换次数
        T det(1);
        std::vector<int> pivot_index(cols,0);
        for (int i = 0 ; i < cols; i++)
        {
            pivot_index[i] = i;          // 理想主元列标，当不存在线性相关的向量时，该数组不会更新
        }
        for (int j = 0 ; j < cols; j ++) 
        {
            int i = pivot_index[j];  // 标记最佳主元
            int flag = pivot_index[j];    // 标记首个非零元，高斯元从此行开始
            if (i >= rows-1)  break;       // 越界，证明cols > rows
            while(i < rows && res.data[i][j] == T()) {   // 寻找首个非零元
                i++;
                flag++;
            }
            if (i == rows)    // 该列没有非零元，证明向量组并非线性无关
            {
                return T(0);
            }
            else              // 该列有非零元
            {
                for (int k = i + 1; k < rows; k++)   // 为了避免数值误差/分数膨胀，选取的主元应该是该列绝对值最大的非零元
                {
                    if (abs(res.at(k,j)) > abs(res.at(i,j)))
                    {
                        i = k;
                    }
                }
                res.swap_line(i, pivot_index[j]);   // 本来应该是j，但是由于存在全零列，所以应该是pivot_index[j]
                swap_times += (i!=pivot_index[j]); // 两行不同时，符号取反;
                for (int k = flag + 1; k < rows; k++)
                {
                    T factor = res.at(k,j)/res.at(pivot_index[j], j);   // 倍加因子
                    for (int l = j; l < cols; l++)          // 高斯消元
                    {
                        res.at(k,l) = res.at(k,l) - res.at(pivot_index[j], l)*factor;
                    }
                }
            }
        }
        for (int i = 0 ; i < rows; i++)
        {
            det = det*res.at(i,i);
        }
        if (swap_times%2 == 1)
        {
            det = -det;
        }
        
        return det;
    } 
}

template<typename T>
Matrix<T> Matrix<T>::to_inv() const {
    // TODO: 实现逆矩阵
    throw std::logic_error("to_inv() not implemented");
}


