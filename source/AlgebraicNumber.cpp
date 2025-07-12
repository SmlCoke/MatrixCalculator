#include "../head/AlgebraicNumber.h"



//-------------------------------------- 根式项 -----------------------------------------

// 缺省构造函数
RadicalTerm::RadicalTerm(): coefficient(0), base(1), exponent(1) {} // 默认 0×1^1


// 构造函数
RadicalTerm::RadicalTerm(const Rational& coeff, const Rational& base, const Rational& exp) {
    if (base<Rational(0) && exp.get_num() % 2 == 1 && exp.get_denom() % 2 == 0 ) {
        throw std::out_of_range("Negative number cannot be raised to an even power ");
    }
    this->coefficient = coeff;
    this->base = base;
    this->exponent = exp;
    simplify();
}


// 化简函数，将特殊案例的无理项化为标准情况：cofe*1^1
void RadicalTerm::simplify() {
    if (coefficient == Rational(0)) {
        base = Rational(1);
        exponent = Rational(1);
        return;
    }
    if (exponent == Rational(0)) {
        // a * b^0 = a * 1 = a
        base = Rational(1);
        exponent = Rational(1);
        return;
    }
    if (base == Rational(1)) {
        exponent = Rational(1); // 1^r = 1
    }
}


// 比较是否为同类项
bool RadicalTerm::is_like(const RadicalTerm& other) const {
    if (this->base == other.base && this->exponent == other.exponent) {
        return true;
    }
    return false;
}


// 计算近似值
double RadicalTerm::to_double() const {
    return this->coefficient.to_double() * std::pow(this->base.to_double(), this->exponent.to_double());
}

// 判断是否为有理数
bool RadicalTerm::is_constant() const {
    if (exponent == 0) return true;
    return false;
}

// 重载等于
bool RadicalTerm::operator==(const RadicalTerm& other) const {
    return coefficient == other.coefficient && base == other.base && exponent == other.exponent;
}

// 重载小于
bool RadicalTerm::operator<(const RadicalTerm& other) const {
    return to_double() < other.to_double();
}

// 重载大于
bool RadicalTerm::operator>(const RadicalTerm& other) const {
    return to_double() > other.to_double();
}

// 重载输出
std::ostream& operator<<(std::ostream& os, const RadicalTerm& rt) {
    if (rt.coefficient == 0 || rt.base == 0) {
        os << 0;
        return os;
    }
    if (rt.base == 1 || rt.exponent == 0) {
        os << rt.coefficient;
        return os;
    }
    os << rt.coefficient << "*(" << rt.base << ")^(" << rt.exponent << ")";
    return os; 
}

