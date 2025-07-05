#include <iostream>
#include "head/Rational.h"
#include "head/Matrix.h"

int main() {
    Rational det;
    
    std::vector<std::vector<Rational>> data(3, std::vector<Rational>(3, Rational(0)));
    data[0][0] = Rational(1); data[0][1] = Rational(1); data[0][2] = Rational(2);
    data[1][0] = Rational(1); data[1][1] = Rational(2); data[1][2] = Rational(2);
    data[2][0] = Rational(3); data[2][1] = Rational(1); data[2][2] = Rational(0);
    std::cout << "data:\n" <<  data << std::endl;
    // 测试矩阵10：4x4整数方阵
    Matrix<Rational> M10(4, 4);
    M10.at(0, 0) = Rational(2); M10.at(0, 1) = Rational(0); M10.at(0, 2) = Rational(1); M10.at(0, 3) = Rational(3);
    M10.at(1, 0) = Rational(1); M10.at(1, 1) = Rational(2); M10.at(1, 2) = Rational(0); M10.at(1, 3) = Rational(1);
    M10.at(2, 0) = Rational(3); M10.at(2, 1) = Rational(1); M10.at(2, 2) = Rational(2); M10.at(2, 3) = Rational(2);
    M10.at(3, 0) = Rational(1); M10.at(3, 1) = Rational(2); M10.at(3, 2) = Rational(3); M10.at(3, 3) = Rational(1);
    std::cout << "matrix example 10 (4x4):\n" << M10 << std::endl;
    std::cout << "Hermite:\n" << M10.get_Hermite() << std::endl;
    std::cout << "Hermite:\n" << M10.get_Hermite() << std::endl;
    std::cout << "det = " << M10.get_det() << std::endl;
    std::cout << "det = " << M10.get_det() << std::endl;
    std::cout << "inv:\n" << M10.get_inv() << std::endl;
    std::cout << "inv:\n" << M10.get_inv() << std::endl;
    std::cout << "rank:\n" << M10.get_rank() << std::endl;
    std::cout << "rank:\n" << M10.get_rank() << std::endl;

    
    return 0;
}