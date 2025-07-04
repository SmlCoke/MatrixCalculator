#ifndef RATIONAL_HPP
#define RATIONAL_HPP

#include <iostream>
#include <numeric>


// Rational operator+(const Rational&) const;  最后的const代表该成员函数不可修改数据成员
class Rational {
private:
    long long num, denom;
    void reduce();           // 约分

public:
    Rational(long long n = 0, long long d = 1);  // 构造函数：基于分子和分母  
    Rational operator+(const Rational&) const;  // 重载加法
    Rational operator-(const Rational&) const;  // 重载减法
    Rational operator-() const;                 // 重载负号 
    Rational operator*(const Rational&) const;  // 重载乘法
    Rational operator/(const Rational&) const;  // 重载除法
    Rational to_reciprocal() const;
    bool operator>(const Rational&) const;      // 重载大于
    bool operator<(const Rational&) const;      // 重载小于
    bool operator==(const Rational&) const;     // 重载等于

    double to_double() const;                   // 化为双精度浮点数
    friend std::ostream& operator<<(std::ostream&, const Rational&);  // 重载输出数据流
    friend Rational abs(const Rational&);       // 取绝对值
    static const Rational Zero;  // 声明静态常量成员
};

#endif