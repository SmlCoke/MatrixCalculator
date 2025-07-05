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
    Matrix<Rational> M9 = build_from_data(data);
    std::cout << "M9: \n" << M9 << std::endl;
    std::cout << "Hermite:\n" << M9.get_Hermite() << std::endl;
    std::cout << "Hermite:\n" << M9.get_Hermite() << std::endl;
    std::cout << "det = " << M9.get_det() << std::endl;
    std::cout << "det = " << M9.get_det() << std::endl;
    std::cout << "inv:\n" << M9.get_inv() << std::endl;
    std::cout << "inv:\n" << M9.get_inv() << std::endl;
    std::cout << "rank:\n" << M9.get_rank() << std::endl;
    std::cout << "rank:\n" << M9.get_rank() << std::endl;
    std::cout << "T:\n" << M9.get_Transpose() << std::endl;
    
    return 0;
}