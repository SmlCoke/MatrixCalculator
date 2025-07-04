// Matrix.hpp
#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <vector>
#include <stdexcept>
#include <iostream>

template<typename T>
class Matrix {
private:
    std::vector<std::vector<T>> data;   // 基于vector实现的矩阵
    int rows, cols;                     // 行数与列数

public:
    // 构造函数：创建 r 行 c 列的矩阵，元素默认初始化为 T()
    Matrix(int r, int c) : rows(r), cols(c), data(r, std::vector<T>(c, T())) {}
    // vector构造函数： vector(len, type)
    // T()：表示默认构造函数

    // 获取行列数
    int row_count() const { return rows; }
    int col_count() const { return cols; }

    // 访问元素
    T& at(int i, int j)  {
        if (i < 0 || i >= rows || j < 0 || j >= cols)
            throw std::out_of_range("Matrix indices out of range");
        return data[i][j];
    }
    // 返回的是引用，可以通过此函数修改矩阵元素

    const T& at(int i, int j) const {
        if (i < 0 || i >= rows || j < 0 || j >= cols)
            throw std::out_of_range("Matrix indices out of range");
        return data[i][j];
    }
    // 返回的是常量引用，不能通过此函数修改矩阵元素
    // 函数末尾const修饰，表示这是一个常成员函数，函数体内部不允许修改数据成员变量

    


    // 加法操作
    Matrix<T> operator+(const Matrix<T>& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for addition");

        Matrix<T> result(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result.at(i, j) = this->at(i, j) + other.at(i, j);
        return result;
    }

    // 减法操作
    Matrix<T> operator-(const Matrix<T>& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for addition");
        Matrix<T> result(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result.at(i, j) = this->at(i, j) - other.at(i, j);
        return result;
    }

    // 取反操作
    Matrix<T> operator-()  const {
        Matrix<T> result(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                result.at(i, j) = -this->at(i, j); 
        return result;
    }

    // 输出显示（方便测试）
    void display() const {
        for (const auto& row : data) {
            std::cout << "|" << " ";
            for (const auto& elem : row) {
                std::cout << elem << " ";
            }
            std::cout << "|";
            std::cout << std::endl;
        }
        std::cout << std::endl; 
    }

    // 乘法操作
    Matrix<T> operator*(const Matrix<T>& other) const {
        if (this->cols != other.rows) {
            throw std::invalid_argument("The left multiplier and right multiplier rows and columns do not match");
        }
        else {
            Matrix<T> results(this->rows, other.cols);
            for (int i = 0; i < this->rows; i++) {
                for (int j = 0 ; j < other.cols; j++) {
                    results.at(i, j) = T();   // 不能写为results.at(i, j) = 0，否则会遇到拷贝构造函数
                    for (int k = 0 ; k < this-> cols; k++) {
                        results.at(i, j) = results.at(i, j) + this->data[i][k]*other.data[k][j];
                    }
                }
            }
            return results;
        }
        
    }

    void swap_line(int i, int j)   // 对换变换
    {
        if (i<0 || i >= rows || j < 0 || j >= rows)   throw std::out_of_range("Matrix indices out of range");
        std::vector<T> temp(cols, T());
        temp = data[i];
        data[i] = data[j];
        data[j] = temp;
    }

    void Multip(int index, const T& r, int j)  { // 对第index行，第j位元素及其之后的元素进行倍乘变换
        if (index < 0 || index  >= rows || j < 0 || j >= rows)   throw std::out_of_range("Matrix indices out of range");
        for (int k = j ; k < cols; k++ ) data[index][k] = data[index][k] * r;
    }


    void Gaussian_Elimination()   { // 高斯消元（只包含前向步骤）
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
            while(i < rows && data[i][j] == T()) {   // 寻找首个非零元
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
                    if (abs(data[k][j]) > abs(data[i][j]))
                    {
                        i = k;
                    }
                }
                *this.swap_line(i, pivot_index[j]);   // 本来应该是j，但是由于存在全零列，所以应该是pivot_index[j]
                for (int k = flag + 1; k < rows; k++)
                {
                    T factor = data[k][j]/ data[pivot_index[j]][j];   // 倍加因子
                    for (int l = j; l < cols; l++)          // 高斯消元
                    {
                        data[k][l] = data[k][l] - data[pivot_index[j]][l]*factor;
                    }
                }
            }
        }

    }
    Matrix<T> to_hermite() const {            // 化为厄尔米特标准型
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

    Matrix<T> to_hermite_mini() const {            // 化为厄尔米特标准型(简易版，即不选取主元)
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

    T get_det() const {
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
            if (swap_times%2 == 1)          // 兑换变换次数为奇数则行列式取反
            {
                det = -det;
            }
            
            return det;
        } 
    }

    Matrix<T> to_inv() const {
        if (this->rows != this->cols) {
            throw std::invalid_argument("Only square matrices have inverse_matrix");
        }
        
        // 也可以：Matrix<T> res = identity(rows); // 省略类名，类作用域内也可用
        
        Matrix<T> inv = Matrix<T>::identity(this->rows);  // 创建单位矩阵
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
            if (i == rows)    // 该列没有非零元，证明该矩阵行列式=0，无逆矩阵
            {
                throw std::invalid_argument("The determinant is 0, and there is no inverse matrix");
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
                res.swap_line(i, j);   // 本来应该是j，但是由于存在全零列，所以应该是pivot_index[j]
                inv.swap_line(i, j);
                for (int k = flag + 1; k < rows; k++)
                {
                    T factor = res.at(k,j)/res.at(pivot_index[j], j);   // 倍加因子
                    for (int l = j; l < cols; l++)          // 高斯消元
                    {
                        res.at(k,l) = res.at(k,l) - res.at(pivot_index[j], l)*factor;
                        
                    }
                    for (int l = 0; l < cols; l++) {
                        inv.at(k,l) = inv.at(k,l) - inv.at(pivot_index[j], l)*factor;
                    }
                }
            }  
        }  // 到此已经化为了Hermite Normal Form
        for (int i = 0; i < rows; i++)   {  // 将主元化为1
            T pivot = res.at(i, i).to_reciprocal();  // 缩放因子
            res.Multip(i, pivot,i);
            inv.Multip(i, pivot,0);
        }        
        // res.display();
        // inv.display();  
        for (int i = rows - 1; i >=0 ; i--)   {// 高斯消元的后向步骤
            for (int j = i - 1; j >= 0; j --) {
                T factor = res.at(j, i);            // 缩放因子
                res.at(j,i) = res.at(j,i) - factor;
                // res.display();
                for (int k = 0; k < rows; k++) {
                    inv.at(j, k) = inv.at(j, k) - factor*inv.at(i, k); 
                    // inv.display();
                }
            }
        }
        return inv;  
    }

    // 静态成员函数：生成 n×n 单位矩阵
    static Matrix<T> identity(int n) {
        Matrix<T> result(n, n);
        for (int i = 0; i < n; ++i) {
            result.at(i, i) = T(1);
        }
        return result;
    }

    // 静态成员函数：生成 r×c 全零矩阵
    static Matrix<T> zeros(int r, int c) {
        return Matrix<T>(r, c); // 默认构造已全零
    }
};

#endif

