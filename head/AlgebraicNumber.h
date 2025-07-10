#ifndef ALGEBRAICNUMBER_H
#define ALGEBRAICNUMBER_H
#include "Rational.h"
#include <vector>
// “代数数”类
// 代数数的解析定义为任何整系数多项式的复数根，体现为有理数与根式项组的加和
// 本类仅仅实现实代数数


#pragma once
#include <iostream>
#include <cmath>

struct RadicalTerm {
    Rational coefficient;  // 系数 c
    Rational base;         // 被开方数/底数 b，通常为正
    Rational exponent;     // 指数/根次数的倒数 r，支持如1/2, 2/3, -1/4

    // 构造函数
    RadicalTerm();
    RadicalTerm(const Rational& coeff, const Rational& base, const Rational& exp);

    //  比较：判断是否为“同类项”（相同 base & exp）
    bool is_like(const RadicalTerm& other) const;

    //  简化：指数为0 ⇒ 变为常数；指数为1 ⇒ 合并到底数
    void simplify();

    //  double 近似值（用于排序或输出）
    double to_double() const;

    //  是否为常数项（即 exponent == 0）
    bool is_constant() const;

    //  重载运算符（用于排序、比较、等价）
    bool operator==(const RadicalTerm& other) const;
    bool operator<(const RadicalTerm& other) const;

    //  重载输出
    friend std::ostream& operator<<(std::ostream& os, const RadicalTerm& t);
};


class AlgebraicNumber {
private:
    std::vector<RadicalTerm> terms;  // 所有根式项线性组合

    void combine_like_terms();  // 合并相同 base 和 degree 的项
    void simplify_all();        // 所有项调用 simplify()

public:
    AlgebraicNumber();  // 默认 0
    AlgebraicNumber(const Rational& r);  // 纯有理数

    void add_term(const RadicalTerm& term);

    // 四则运算
    AlgebraicNumber operator+(const AlgebraicNumber& other) const;
    AlgebraicNumber operator-(const AlgebraicNumber& other) const;
    AlgebraicNumber operator*(const AlgebraicNumber& other) const;
    AlgebraicNumber operator/(const AlgebraicNumber& other) const;

    // 其他
    bool is_rational() const;
    double to_double() const;

    friend std::ostream& operator<<(std::ostream& os, const AlgebraicNumber& a);
};

#endif