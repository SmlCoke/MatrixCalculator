#include <iostream>
#include "head/Rational.h"
#include "head/Matrix.h"

int main() {
    Rational det;
    // 测试矩阵10：4x4整数方阵
    Matrix<Rational> M10(4, 4);
    M10.at(0, 0) = Rational(2); M10.at(0, 1) = Rational(0); M10.at(0, 2) = Rational(1); M10.at(0, 3) = Rational(3);
    M10.at(1, 0) = Rational(1); M10.at(1, 1) = Rational(2); M10.at(1, 2) = Rational(0); M10.at(1, 3) = Rational(1);
    M10.at(2, 0) = Rational(3); M10.at(2, 1) = Rational(1); M10.at(2, 2) = Rational(2); M10.at(2, 3) = Rational(2);
    M10.at(3, 0) = Rational(1); M10.at(3, 1) = Rational(2); M10.at(3, 2) = Rational(3); M10.at(3, 3) = Rational(1);
    std::cout << "matrix example 10 (4x4):\n" << M10 << std::endl;
    std::cout << "Hermite:\n" << M10.to_hermite() << std::endl;
    std::cout << "Hermite:\n" << M10.to_hermite() << std::endl;
    std::cout << "det = " << M10.get_det() << std::endl;
    std::cout << "det = " << M10.get_det() << std::endl;
    std::cout << "A * A':\n" << M10*M10.to_inv() << std::endl;
    
    return 0;
}