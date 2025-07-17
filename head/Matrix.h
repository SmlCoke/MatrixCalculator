// Matrix.h
#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>
#include <iostream>

template<typename T>
class Matrix {
private:
    std::vector<std::vector<T>> data;   // 基于vector实现的矩阵
    int rows, cols;                     // 行数与列数

    
    // ---------- 缓存部分 ----------
    // 行列式
    mutable bool det_ready = false;
    mutable T cached_det;

    // 秩
    mutable bool rank_ready = false;
    mutable int cached_rank;


    // Hermite Normal Form
    mutable bool Hermite_ready = false;
    mutable std::vector<std::vector<T>> cached_Hermite;

    // inv
    mutable bool inv_ready = false;
    mutable std::vector<std::vector<T>> cached_inv;

    // RREF
    mutable bool RREF_ready = false;
    mutable std::vector<std::vector<T>> cached_RREF; 
    
    // 失效机制
    void invalidate_cache() const {
        det_ready = false;
        rank_ready = false;
        Hermite_ready = false;
        inv_ready = false;
        RREF_ready = false;
        // 其他缓存标志也在这里置为 false
    }
    

    // 对换变换
    void swap_line(int i, int j);   
    
    // 对第index行，第j位元素及其之后的元素进行倍乘变换
    void Multip(int index, const T& r, int j); 

    // 高斯消元（只包含前向步骤）
    // 带返回值，返回的是对换变换次数
    int Gaussian_Elimination();     // 私有函数，因为是直接操控数据成员 

public:
    // 构造函数：创建 r 行 c 列的矩阵，元素默认初始化为 T()
    Matrix(int r, int c) : rows(r), cols(c), data(r, std::vector<T>(c, T())) {
        cached_Hermite = std::vector(r, std::vector<T>(c, T()));
        cached_inv = std::vector(r, std::vector<T>(c, T()));
        invalidate_cache();
    }
    // vector构造函数： vector(len, type)
    // T()：表示默认构造函数

    // 获取行列数
    int row_count() const { return rows; }
    int col_count() const { return cols; }

    // 访问元素
    T& at(int i, int j)  {
        if (i < 0 || i >= rows || j < 0 || j >= cols)
            throw std::out_of_range("Matrix indices out of range");
        // invalidate_cache(); // 缓存失效   
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
    Matrix<T> operator+(const Matrix<T>& other) const;

    // 减法操作
    Matrix<T> operator-(const Matrix<T>& other) const;

    // 取反操作
    Matrix<T> operator-()  const;

    // 输出显示（方便测试）
    void display() const;

    // 乘法操作
    Matrix<T> operator*(const Matrix<T>& other) const;

    
    // 化为厄尔米特标准型
    std::vector<std::vector<T>> to_hermite() const;            
    
    // 化为厄尔米特标准型(简易版，即不选取主元)
    std::vector<std::vector<T>> to_hermite_mini() const;           


    // 求逆
    std::vector<std::vector<T>> to_inv() const;

    // 求简化阶梯型RREF
    std::vector<std::vector<T>> to_RREF() const;

    // 静态成员函数：生成 n×n 单位矩阵
    static Matrix<T> identity(int n);

    // 静态成员函数：生成 r×c 全零矩阵
    static Matrix<T> zeros(int r, int c);
    
    // 求行列式
    T get_det() const;

    // 求秩
    int get_rank() const;

    // 求Hermite Normal Form
    std::vector<std::vector<T>> get_Hermite() const;
    
    // 求逆矩阵
    std::vector<std::vector<T>> get_inv() const;

    // 求转置矩阵
    Matrix<T> get_Transpose() const;

    // 求行简化阶梯型
    std::vector<std::vector<T>> get_RREF() const;

    // 求解线性方程组Ax=b（方阵情况）
    std::vector<T> pha_solve(const std::vector<T> & Din) const;

    // 求零空间（即线性齐次方程组的解空间）
    std::vector<std::vector<T>> get_null() const;

    // 求解空间
    std::vector<std::vector<T>> solve(const std::vector<T>& Din) const;

    
};

// 基于矩阵内容构建矩阵对象
template<typename T>
Matrix<T> build_from_data(const std::vector<std::vector<T>>& m_data ) {
    if (m_data.size() == 0 || m_data[0].size() == 0) {
        throw std::invalid_argument("Cannot build matrix form empty");
    }
    Matrix<T> res(m_data.size(), m_data[0].size());
    for (int i = 0; i < res.row_count(); i++) {
        for (int j = 0; j < res.col_count(); j++) {
            res.at(i, j) = m_data[i][j];
        }
    }
    return res;
} 


// 加法操作
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


// 减法操作
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

// 取反操作
template<typename T>
Matrix<T> Matrix<T>::operator-()  const {
    Matrix<T> result(rows, cols);
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            result.at(i, j) = -this->at(i, j); 
    return result;
}

// 输出显示（方便测试）
template<typename T>
void Matrix<T>::display() const {
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
template<typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& other) const {
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

// 对换变换
template<typename T>
void Matrix<T>::swap_line(int i, int j)   
{
    if (i<0 || i >= rows || j < 0 || j >= rows)   throw std::out_of_range("Matrix indices out of range");
    std::vector<T> temp(cols, T());
    temp = data[i];
    data[i] = data[j];
    data[j] = temp;
    // invalidate_cache(); // 缓存失效
}

// 对第index行，第j位元素及其之后的元素进行倍乘变换
template<typename T>
void Matrix<T>::Multip(int index, const T& r, int j)  { 
    if (index < 0 || index  >= rows || j < 0 || j >= rows)   throw std::out_of_range("Matrix indices out of range");
    for (int k = j ; k < cols; k++ ) data[index][k] = data[index][k] * r;
    // invalidate_cache(); // 缓存失效
}

// 高斯消元（只包含前向步骤）
template<typename T>
int Matrix<T>::Gaussian_Elimination()   {         
    std::vector<int> pivot_index(cols,0);      // 理想主元列标
    int swap_times = 0; // 纪录对换变换次数
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
            this->swap_line(i, pivot_index[j]);   // 本来应该是j，但是由于存在全零列，所以应该是pivot_index[j]
            swap_times += (i != pivot_index[j]);
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
    
    return swap_times;
}

// 化为厄尔米特标准型
// 完成该步骤之后，已经可以求出行列式以及秩，因此本函数的功能不仅限于高斯消元，
// 同时还会设置缓存
template<typename T>
std::vector<std::vector<T>> Matrix<T>::to_hermite() const {            
    Matrix<T> res = *this;   // 独立出一个新矩阵
    res.Gaussian_Elimination();
    return res.data;
}

// 化为厄尔米特标准型(简易版，即不选取主元)
template<typename T>
std::vector<std::vector<T>> Matrix<T>::to_hermite_mini() const {            
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
    return res.data;
}




// 求逆
template<typename T>
std::vector<std::vector<T>> Matrix<T>::to_inv() const {
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
        T pivot = T(1)/res.at(i, i);  // 缩放因子
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
    return inv.data;  
}

// 静态成员函数：生成 n×n 单位矩阵
template<typename T>
Matrix<T> Matrix<T>::identity(int n) {
    Matrix<T> result(n, n);
    for (int i = 0; i < n; ++i) {
        result.at(i, i) = T(1);
    }
    return result;
}

// 静态成员函数：生成 r×c 全零矩阵
template<typename T> 
Matrix<T> Matrix<T>::zeros(int r, int c) {
    return Matrix<T>(r, c); // 默认构造已全零
}

// 重载输出数据流
template<typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T>& m) {
    for (int i = 0; i < m.row_count(); ++i) {
        os << "| ";
        for (int j = 0; j < m.col_count(); ++j) {
            os << m.at(i, j) << " ";
        }
        os << "|" << std::endl;
    }
    os << std::endl;
    return os;
}

// 重载std::vector<std::vector<T>>的输出数据流
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<std::vector<T>> & m_data) {
    int rows = m_data.size();
    int cols = m_data[0].size();
    for (int i = 0; i < rows; ++i) {
        os << "| ";
        for (int j = 0; j < cols; ++j) {
            os << m_data[i][j] << " ";
        }
        os << "|" << std::endl;
    }
    os << std::endl;
    return os;
}

// 重载std::vector<T>的输出数据流
template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T> & v_data) {
    int sizes = v_data.size();
    for (int i = 0; i < sizes; ++i) {
        os << v_data[i] << " ";
    }
    return os;
}

// 求Hermite Normal Form
template<typename T>
std::vector<std::vector<T>> Matrix<T>::get_Hermite() const { 
    if (Hermite_ready) return cached_Hermite;  // 惰性计算
    Hermite_ready = true;
    cached_Hermite = to_hermite();

    // 缓存行列式
    T det(1);
    if (rows == cols) {
        for (int i = 0; i < rows; i++ ) {
            det = det*cached_Hermite[i][i];
        }
    }
    det_ready = true;
    cached_det = det;

    // 缓存秩
    int rank = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) { 
            if (!(cached_Hermite[i][j] == T())) {
                rank++;
                break;
            }
        }
    }
    cached_rank = rank;
    rank_ready = true;

    return cached_Hermite;
}

// 求逆矩阵
template<typename T>
std::vector<std::vector<T>> Matrix<T>::get_inv() const {
    if (inv_ready) return cached_inv;   // 惰性计算
    inv_ready = true;
    cached_inv = to_inv();
    return cached_inv;

} 

// 求行列式
template<typename T>
T Matrix<T>::get_det() const {
    if (rows != cols)
    {
        throw std::invalid_argument("Only square matrices can be used to find the determinant");
    }
    else
    {
        if (det_ready) return cached_det;   // 惰性计算
        Matrix<T> res = *this;   // 独立出一个新矩阵
        int swap_times = 0;      // 纪录对换变换次数
        T det(1);
        swap_times = res.Gaussian_Elimination();
        for (int i = 0 ; i < rows; i++)
        {
            det = det*res.at(i,i);
        }
        if (swap_times%2 == 1)          // 兑换变换次数为奇数则行列式取反
        {
            det = -det;
        }
        cached_det = det;
        det_ready = true;

        // 缓存Hermite Norm Form
        Hermite_ready = true;
        cached_Hermite = res.data;

        // 缓存秩
        int rank = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) { 
                if (!(cached_Hermite[i][j] == T())) {
                    rank++;
                    break;
                }
            }
        }
        cached_rank = rank;
        rank_ready = true;

        return det;
    } 
}

// 求秩
template<typename T>
int Matrix<T>::get_rank() const {
    if (rank_ready) return cached_rank;
    std::vector<std::vector<T>> Hermite = get_Hermite(); // 调用厄尔米特标准型后，自动完成秩的缓存
    return cached_rank;
}

// 求转置矩阵
template<typename T>
Matrix<T> Matrix<T>::get_Transpose() const {
    Matrix<T> res(cols, rows);
    for (int i = 0; i < cols; i ++) {
        for (int j = 0 ; j < rows; j ++) {
            res.at(i, j) = data[j][i];
        }
    }
    return res;
}



// 求解线性方程组Ax=b(对于方阵)
template<typename T>
std::vector<T> Matrix<T>::pha_solve(const std::vector<T> & Din) const {
    if (this->rows != this->cols) throw std::invalid_argument("This func can only process phalanx case.");
    if (Din.size() != this->cols) throw std::invalid_argument("N of cols of Matrix != length of input.");
    std::vector<std::vector<T>>inv = this->get_inv();
    std::vector<T> solution(rows, T(0));
    for (int i = 0 ; i < rows; i ++) {
        for (int j = 0; j < cols; j++) {
            solution[i] = solution[i] + inv[i][j] * Din[j];
        }
    }
    return solution;
    
}


// 计算简化阶梯型RREF
template<typename T>
std::vector<std::vector<T>> Matrix<T>::to_RREF() const {
    // int flag = cols;  // 标记首个非零行
    std::vector<std::vector<T>> res = get_Hermite();  // 从厄尔米特标准型开始化简
    for (int i = rows - 1; i >=0; i--) {
        for (int j = i; j < cols; j++) {
            if (!(res[i][j]==T())) {// 找到非零元
                for (int l = j + 1; l < cols; l++) {
                    res[i][l] = res[i][l]/res[i][j];   // 将主元化为1
                }
                res[i][j] = T(1);
                for (int k = i - 1; k >=0; k--) {  // 高斯消元的后向步骤
                    T factor = res[k][j];  // 倍加因子
                    for (int l = j; l < cols; l++) {
                        res[k][l] = res[k][l] - factor * res[i][l];
                    }
                }
                break;
            }
            // if (j == cols - 1) { // 找到末尾还没找到非零元
            //     flag--;
            // }
        }
    }
    // cached_rank = flag;   // 首个全零行索引即非零行个数，也即矩阵的秩
    // rank_ready = true;
    return res;
}

// 求简化阶梯型RREF
template<typename T>
std::vector<std::vector<T>> Matrix<T>::get_RREF() const {
    if (RREF_ready) return cached_RREF;   // 惰性计算
    cached_RREF = to_RREF();
    RREF_ready = true;
    return cached_RREF;
}


// 求零空间（即线性齐次方程组的解空间）
template<typename T>
std::vector<std::vector<T>> Matrix<T>::get_null() const {
    int rank = get_rank();
    if (rank == cols) return std::vector(cols,std::vector(1,T()));            
    int n_null = cols - rank;        // 基础解系维度 = n - 列空间维度   
    std::vector<std::vector<T>> res = std::vector(n_null,std::vector(cols,T()));
    std::vector<std::vector<T>> rref = get_RREF();
    std::vector<int> pivot_flag(cols, 0);   // 标记该变量是否是主元
    std::vector<int> pivot_index(rank, 0);   // 标记该行的主元索引
    for (int i = 0; i < rank; i++) {
        int j = i;
        while(rref[i][j]==T()) j++;
        pivot_flag[j] = 1;   // 记录：该变量是主元
        pivot_index[i] = j;  // 记录：该行的主元位置
    }
    int res_index = 0; // 开始构造基础解系
    for (int j = 0 ; j < cols; j++) {
        if (pivot_flag[j] == 0) {
            res[res_index][j] = T(1);
            for (int i = 0 ; i < rank; i++) {
                res[res_index][pivot_index[i]] = -rref[i][j];
            }
            res_index++;
        }
    }
    return res;
}

// 求解空间
template<typename T>
std::vector<std::vector<T>> Matrix<T>::solve(const std::vector<T>& Din) const {
    if (Din.size() != rows) throw std::invalid_argument("b must have the same length as matrix.");
    Matrix<T> Aug_ma(rows, cols+1);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            Aug_ma.at(i,j) = at(i,j);
        }
        Aug_ma.at(i,cols) = Din[i];
    }
    // std::cout << "Aug_ma:\n" << Aug_ma << std::endl;
    std::vector<std::vector<T>> rref = Aug_ma.get_RREF();
    // std::cout << "Aug_rref:\n" << rref << std::endl;
    

    // 检测方程组是否无解
    bool flag = true;
    int rank = 0;
    for (int i = 0; i < rows; i++) {
        int j = 0;
        for (j = 0; j < cols; j++) {
            if (!(rref[i][j] == T())) {
                rank++;
                break;
            }
        }
        if (j == cols && !(rref[i][cols] == T())) {
            flag = false;
            break;
        }
     }

    if (!flag) return {};

    // 初始化解空间
    int n_null = cols - rank;
    std::vector<std::vector<T>> res(n_null+1, std::vector<T>(cols, T())); 
    std::vector<int> pivot_flag(cols, 0);   // 标记该变量是否是主元
    std::vector<int> pivot_index(rank, 0);   // 标记该行的主元索引
    for (int i = 0; i < rank; i++) {
        int j = i;
        while(rref[i][j]==T()) j++;
        pivot_flag[j] = 1;   // 记录：该变量是主元
        pivot_index[i] = j;  // 记录：该行的主元位置
    }

    // 特解
    for (int i = 0; i < rank; i++) {
        res[0][pivot_index[i]] = rref[i][cols];
    }

    // 齐次解
    int res_index = 1; // 开始构造基础解系
    for (int j = 0 ; j < cols; j++) {
        if (pivot_flag[j] == 0) {
            res[res_index][j] = T(1);
            for (int i = 0 ; i < rank; i++) {
                res[res_index][pivot_index[i]] = -rref[i][j];
            }
            res_index++;
        }
    }
    return res;
}

#endif