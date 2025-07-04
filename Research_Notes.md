# Research_Notes

## 一、C++ 基础

### static_cast<目标类型>(表达式)
static_cast是C++推荐的显式类型转换方式，比C风格强制类型转换更安全、可读性更好。
适用于基本数据类型之间的转换（如int转double）、指针类型之间的转换等。
编译时检查，能防止不安全的转换。

## 二、Class

### private, protected, public
1. public（公有）
含义：公有成员可以被**类的外部、类的成员函数、以及子类直接访问**。
用途：对外提供接口（如构造函数、常用操作函数）。
```C++
class A {
public:
    int x;
};
// 外部可以直接访问 a.x
```
2. protected（受保护）
含义：受保护成员不能被类的外部直接访问，但可以**被本类和其子类访问**。
用途：用于希望子类能访问但不希望外部访问的成员。
```C++
class A {
protected:
    int y;
};
class B : public A {
    void func() { y = 5; } // 子类可以访问
};
// 外部不能访问 a.y
```
3. private（私有）
含义：私有成员**只能被本类的成员函数**访问，不能被类的外部或子类访问。
用途：用于封装实现细节，保护数据安全。
```C++
class A {
private:
    int z;
};
// 外部和子类都不能访问 a.z
```

### const
Rational operator+(const Rational&) const;  
最后的const代表该成员函数**不可修改数据成员**


> me：
> T& at(int i, int j) {
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

> copilot：
> 1. 非const对象调用
> Matrix<int> m(2, 2);
m.at(0, 1) = 5; // 调用 T& at(int i, int j)
> 2. const对象调用
> const Matrix<int> m(2, 2);
int x = m.at(0, 1); // 调用 const T& at(int i, int j) const
> 3. 自动区分原理
如果对象是const，编译器只允许调用带const修饰的成员函数。
如果对象是非const，优先调用非const成员函数；如果没有非const版本，也可以调用const版本。




### 构造函数
构造函数参数默认值不要和无参构造混用，否则会二义性。
**推荐写法**：
> class Rational {
public:
    Rational(); // 只声明无参构造
    Rational(long long n, long long d); // 明确两个参数
    Rational(long long n); // 明确一个参数
    // ...其他成员...
};

**错误写法**：
> class Rational {
public:
    Rational(); // 只声明无参构造
    Rational(long long n = 0, long long d = 1); // 明确两个参数
    Rational(long long n = 0); // 明确一个参数
    // ...其他成员...
};



### 拷贝构造函数？
> me
> 构造基于有理数实现的矩阵类的这个函数：
// 乘法操作
Matrix<T> operator*(const Matrix<T>& other) {
if (this->cols != other.rows) {
throw std::invalid_argument("The left multiplier and right multiplier rows and columns do not match");
}
Matrix<T> results(this->rows, other.cols);
else {
for (int i = 0; i < this->rows; i++) {
results.at(i, j) = 0;
for (int j = 0 ; j < this->cols; j++) {
            }
        }
    }
}
是不是要求首先要实现Rational类的赋值构造函数？

> copilot:
> 不需要专门实现Rational类的赋值构造函数（即拷贝构造函数），因为C++会自动生成一个默认的拷贝构造函数和赋值操作符，**只要你的成员变量（如long long类型）本身支持拷贝和赋值即可**。
**你需要关注的点**：
Rational类的拷贝构造和赋值操作

1. 如果你没有自定义，编译器会自动生成，能正常用于矩阵乘法中的赋值和拷贝。
只有当你的类有特殊资源管理（如**指针、动态内存**）时，才需要自定义。
矩阵乘法实现

2. 你只需要保证Rational类**有无参构造（如Rational()），能被初始化为0**。
只要Rational支持加法、乘法、赋值，矩阵乘法就能正常实现。

### 重载abs函数
我定义了友元函数abs()之后，不会与std::abs冲突吗？
ChatGPT:
> 如果你定义了一个友元函数 abs(const Rational&)，不会与 std::abs 冲突，只要你不把它放进 std:: 命名空间。
在 C++ 中，当你调用 abs(x) 时，编译器使用：
🔍 Koenig 查找（ADL, Argument-Dependent Lookup）+ 普通作用域查找
>具体查找顺序是：
>1. 当前作用域中是否有 abs() 函数
>2. 函数参数类型所在的命名空间中是否有 abs()（即 ADL）
>3. 是否存在可用的 std::abs()（仅当没有更好匹配）
>4. 是否可以进行隐式类型转换

>这样定义以后：
abs(Rational(1, 2)) 会调用你定义的 abs
abs(3.14) 会调用 std::abs(double)
它们不会冲突，因为参数类型不同




## 三、cmake编译时遇到的一些破事
问题面板中的报错突然出现中文乱码，经过copilot后发现是选择的默认编译器时Visual Studio的MSVC

然后我将MINGW64移动到系统环境变量的最上端后，依旧无效：
> copilot：
> 为什么“PROBLEMS”面板依然乱码？
VS Code 的“PROBLEMS”面板显示的是自动编译任务（如 CMake Tools、VS Code 任务、F5 调试等）输出的内容。
如果这些任务用的是 MSVC（cl.exe），就会有中文乱码。
你在终端输入 g++ --version 只是测试了终端的 g++，和 VS Code 的自动编译任务用的编译器不是一回事。

然后**创建cmake项目时直接指定编译器**，最后成功

### 创建cmake项目的方法
1. 在项目根目录下新建 build 文件夹： 
> mkdir build
2. 进入 build 文件夹：
> cd build
3. 在 build 文件夹下执行 CMake，指定上一级目录为源码目录：
> cmake -G "MinGW Makefiles" ..
4. 编译项目：
> cmake --build .



## 四、泛型机制
❗模板类（如 template<typename T> class Matrix）不能将成员函数定义放到 .cpp 文件中使用普通方式编译链接。因为：
✅ 模板类必须在头文件中定义实现（或显式实例化）
这是 C++ 的语法机制决定的 —— 模板是在编译期展开的，只有看到完整定义编译器才能实例化出代码。