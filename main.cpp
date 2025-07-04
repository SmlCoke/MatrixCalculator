#include <iostream>
#include "head/Rational.h"
#include "head/Matrix.h"

int main() {
    Rational det;

    // 测试矩阵8：3x3整数方阵
    Matrix<Rational> M8(2,2);
    M8.at(0, 0) = Rational(1); M8.at(0, 1) = Rational(1); 
    M8.at(1, 0) = Rational(0); M8.at(1, 1) = Rational(1); 
    std::cout << "matrix example 8 (2x2):\n";
    M8.display();
    std::cout << "det = ";
    std::cout << M8.get_det() << std::endl;
    std::cout << "\n";
    (M8*M8.to_inv()).display();

    // 测试矩阵9：3x3分数方阵
    Matrix<Rational> M9(3, 3);
    M9.at(0, 0) = Rational(1, 2); M9.at(0, 1) = Rational(1, 3); M9.at(0, 2) = Rational(1, 4);
    M9.at(1, 0) = Rational(2, 3); M9.at(1, 1) = Rational(2, 5); M9.at(1, 2) = Rational(2, 7);
    M9.at(2, 0) = Rational(3, 4); M9.at(2, 1) = Rational(3, 8); M9.at(2, 2) = Rational(3, 11);
    std::cout << "matrix example 9 (3x3):\n";
    M9.display();
    std::cout << "det = ";
    std::cout << M9.get_det() << std::endl;
    std::cout << "\n";
    (M9.to_hermite()).display();
    (M9*M9.to_inv()).display();
    // 测试矩阵10：4x4整数方阵
    Matrix<Rational> M10(4, 4);
    M10.at(0, 0) = Rational(2); M10.at(0, 1) = Rational(0); M10.at(0, 2) = Rational(1); M10.at(0, 3) = Rational(3);
    M10.at(1, 0) = Rational(1); M10.at(1, 1) = Rational(2); M10.at(1, 2) = Rational(0); M10.at(1, 3) = Rational(1);
    M10.at(2, 0) = Rational(3); M10.at(2, 1) = Rational(1); M10.at(2, 2) = Rational(2); M10.at(2, 3) = Rational(2);
    M10.at(3, 0) = Rational(1); M10.at(3, 1) = Rational(2); M10.at(3, 2) = Rational(3); M10.at(3, 3) = Rational(1);
    std::cout << "matrix example 10 (4x4):\n";
    M10.display();
    std::cout << "det = ";
    std::cout << M10.get_det() << std::endl;
    std::cout << "\n";
    (M10*M10.to_inv()).display();
    Matrix<Rational> M11(2, 2);
    M11.at(0,0) = Rational(1);
    M11.at(0,1) = Rational(2);
    M11.at(1,0) = Rational(4);
    M11.at(1,1) = Rational(3);
    std::cout << "matrix example 11 (4x4):\n";
    M11.display();
    std::cout << "det = ";
    std::cout << M11.get_det() << std::endl;
    std::cout << "\n";
    (M11*M11.to_inv()).display();
    return 0;
}