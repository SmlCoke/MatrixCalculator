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


// 化简函数
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
