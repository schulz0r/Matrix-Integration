#include <iostream>

#include <boost/math/quadrature/gauss.hpp>
#include <unsupported/Eigen/MatrixFunctions>
#include <Eigen/Geometry>

#include "matrices.h"

typedef Eigen::Matrix<double, 12, 12> GramianMatrix;

int main() {
    Matrices LTISystem;
    
    // the integral function
    std::function<GramianMatrix(double)> gramianIntegral = [&](const double & t) -> GramianMatrix {
        return (LTISystem.A * t).exp() * LTISystem.B * LTISystem.B.transpose() * (LTISystem.A.transpose() * t).exp();
    };

    // now do the integration
    boost::math::quadrature::gauss<double, 10> integrator;
    GramianMatrix W_t = integrator.integrate(gramianIntegral, 0.0, 1.0);

    std::cout << W_t << std::endl;
    
    return 0;
}
