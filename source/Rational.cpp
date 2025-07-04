// src/Rational.cpp
#include "../head/Rational.h"

const Rational Rational::Zero(0);  // 初始化静态常量成员


// 构造函数：基于分子与分母
Rational::Rational(long long n, long long d) : num(n), denom(d) {
    if (d==0) {
        throw std::out_of_range("denom cannot be assigned as 0");
    }
    reduce();
}

// 约分
void Rational::reduce() {
    long long g = std::gcd(num, denom);
    if (denom < 0) g = -g;         // 保证分母是正数
    num /= g;
    denom /= g;
}


// 重载："+" 运算符
Rational Rational::operator+(const Rational& rhs) const {
    return Rational(this->num * rhs.denom + rhs.num * this->denom, denom * rhs.denom);
}

// 重载："-" 运算符（二元）
Rational Rational::operator-(const Rational& oper) const {
    return Rational(this->num * oper.denom - this->denom * oper.num, this->denom * oper.denom);
}

// 重载："-" 运算符（一元）
Rational Rational::operator-() const {
    return Rational(-this->num, this->denom);
}

// 重载: "*" 运算符
Rational Rational::operator*(const Rational& oper) const {
    return Rational(this->num * oper.num, this->denom * oper.denom);
}

// 重载： "/"运算符
Rational Rational::operator/(const Rational& oper) const {
    return Rational(this->num * oper.denom, this->denom * oper.num);
}
/*
// 重载： "=="运算符
bool Rational::operator==(const Rational& oper) const {
    if ((this->num*this->denom <0 && oper.num*oper.denom >0) || (this->num*this->denom >0 && oper.num*oper.denom <0)) {
        return false;
    }
    else {
        if ((std::abs(this->num) == std::abs(oper.num)) && (std::abs(this->denom) == std::abs(oper.denom)))
        {
            return true;
        }
        else if (this->num ==0 && oper.num==0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
     
}
*/
bool Rational::operator==(const Rational& oper) const {
    if (this->num * oper.num < 0 ) {
        return false;
    }
    else {
        if ((std::abs(this->num) == std::abs(oper.num)) && (std::abs(this->denom) == std::abs(oper.denom)))
        {
            return true;
        }
        else if (this->num ==0 && oper.num==0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
     
}
// 将分数转化为双精度浮点数
double Rational::to_double() const {
    return static_cast<double>(this->num)/this->denom;
}


// 重载输出数据流运算符
std::ostream& operator<<(std::ostream& os, const Rational& oper) {
    if (oper.num == 0) {    
        os << 0;    // 0 
        return os;
    }
    // 无需判断负号，分母一定是正的
    os << oper.num;           
    if (oper.denom==1) {
        return os;             // 整数
    }
    else {
        os << "/" << oper.denom;
        return os;
    }
}


// 求倒数
Rational Rational::to_reciprocal() const {
    if (this->num == 0) {
        throw std::out_of_range("0 has no reciprocal.");
    }
    return Rational(this->denom, this->num);  // 调用了构造函数，自动将可能存在的负号移动到分子
}

// 重载大于
bool Rational::operator>(const Rational& other) const {
    return this->num*other.denom > this->denom*other.num;
}

// 重载小于
bool Rational::operator<(const Rational& other) const {
    return this->num*other.denom < this->denom*other.num;
}

// 重载abs()
Rational abs(const Rational& other) {
    return (other.num < 0) ? Rational(-other.num, other.denom) : other;
}