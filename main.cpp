#include <iostream>
#include "head/Rational.h"
#include "head/Matrix.h"
#include "head/AlgebraicNumber.h"

int main() {
    /*
    Rational det;
    std::vector<std::vector<Rational>> data(3, std::vector<Rational>(3, Rational(0)));
    data[0][0] = Rational(1); data[0][1] = Rational(1); data[0][2] = Rational(2);
    data[1][0] = Rational(1); data[1][1] = Rational(2); data[1][2] = Rational(2);
    data[2][0] = Rational(3); data[2][1] = Rational(1); data[2][2] = Rational(0);
    std::cout << "data:\n" <<  data << std::endl;
    Matrix<Rational> M9 = build_from_data(data);
    //  std::vector<std::vector<double>> data(3, std::vector<double>(3, 0));
    //  data[0][0] = 1; data[0][1] = 1; data[0][2] = 2;
    //  data[1][0] = 1; data[1][1] = 2; data[1][2] = 2;
    //  data[2][0] = 3; data[2][1] = 1; data[2][2] = 0;
    //  std::cout << "data:\n" <<  data << std::endl;
    //  Matrix<double> M9 = build_from_data(data);
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
    std::vector<Rational> input(3,1);
    std::cout << "solution of [1, 1, 1]:\n" << M9.pha_solve(input) << std::endl;
    input[0]=Rational(4); input[1]=Rational(2); input[2]=Rational(0);
    std::cout << "solution of [4, 2, 0]:\n" << M9.pha_solve(input) << std::endl;
    */
    // std::vector<std::vector<Rational>> data(4, std::vector<Rational>(4, Rational(0)));
    // data[0][0] = Rational(1); data[0][1] = Rational(2); data[0][2] = Rational(3); data[0][3] = Rational(4); 
    // data[1][0] = Rational(0); data[1][1] = Rational(1); data[1][2] = Rational(2); data[1][3] = Rational(3);
    // data[2][0] = Rational(1); data[2][1] = Rational(1); data[2][2] = Rational(1); data[2][3] = Rational(1);
    // data[2][0] = Rational(0); data[2][1] = Rational(0); data[2][2] = Rational(1); data[2][3] = Rational(2);
    // Matrix<Rational> M = build_from_data(data);
    // std::cout << "RREF of M9:\n" << M.get_RREF() << std::endl;
    // std::cout << "Null_solu of M9: \n" << M.get_null() << std::endl;

    std::vector<std::vector<Rational>> data2(4, std::vector<Rational>(4, Rational(0)));
    data2[0][0] = Rational(0); data2[0][1] = Rational(1); data2[0][2] = Rational(3); data2[0][3] = Rational(0); 
    data2[1][0] = Rational(0); data2[1][1] = Rational(0); data2[1][2] = Rational(0); data2[1][3] = Rational(1);
    // data2[2][0] = Rational(-1); data2[2][1] = Rational(-1); data2[2][2] = Rational(-1); data2[2][3] = Rational(-1);
    // data2[2][0] = Rational(0); data2[2][1] = Rational(0); data2[2][2] = Rational(0); data2[2][3] = Rational(0);
    Matrix<Rational> M2 = build_from_data(data2);
    std::cout << "M2:\n" << M2 << std::endl;
    std::cout << "RREF of M2:\n" << M2.get_RREF() << std::endl;
    std::cout << "Null_solu of M2:\n" << M2.get_null() << std::endl;
    std::vector<Rational> b = std::vector<Rational>(4,Rational());
    b[0] = Rational(5); b[1] = Rational(13); b[2] = Rational(0); b[3] = Rational(0);
    std::cout << "solu of M2 and " << b << ":\n" << M2.solve(b) << std::endl;
    return 0;
}