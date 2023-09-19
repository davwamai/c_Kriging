#include "./headers/OrdinaryKriging.h"
#include <cmath>
#include <iostream>
#include <Eigen/Dense>  // For matrix operations
#include <nlopt.hpp>

OrdinaryKriging::OrdinaryKriging(const std::vector<std::vector<double>>& points, 
                                 const std::vector<double>& zvals, 
                                 const std::string& variogram)
    : points_(points), zvals_(zvals), variogram_(variogram) {
    zvals_org_ = zvals_;  // Copy of the original z values
    setupVariogram();
}

void OrdinaryKriging::setupVariogram() {
    if (variogram_ == "gaussian") {
        variogramFunction_ = [](double h, double a, double C) -> double {
            return C * (1 - std::exp(-1 * ((h / a) * (h / a))));
        };
    }
    //TODO: Add other variogram models

    else {
        std::cerr << "No valid variogram model selected" << std::endl;
        exit(1);
    }
}

void OrdinaryKriging::ManualParamSet(double C, double a, double nugget, double anisotropy_factor) {
    a_ = a;
    anisotropy_factor_ = anisotropy_factor;
    C_ = C;
    nugget_ = nugget;
}

void OrdinaryKriging::MatrixSetup() {
    // Compute the pairwise distance matrix
    size_t num_points = points_.size();
    std::vector<std::vector<double>> distances(num_points, std::vector<double>(num_points, 0.0));
    for (size_t i = 0; i < num_points; ++i) {
        for (size_t j = 0; j < num_points; ++j) {
            distances[i][j] = std::sqrt(std::pow(points_[i][0] - points_[j][0], 2) + 
                            std::pow(points_[i][1] - points_[j][1], 2) / (anisotropy_factor_ * anisotropy_factor_));
        }
    }
    // TODO: rest of the method
}

double OrdinaryKriging::SinglePoint(double Xo, double Yo, const std::vector<std::vector<double>>& training_points) {
    std::vector<std::vector<double>> points_to_use = (training_points.empty()) ? points_ : training_points;

    std::vector<double> distances_to_point0(points_to_use.size());
    for (size_t i = 0; i < points_to_use.size(); ++i) {
        distances_to_point0[i] = std::sqrt(std::pow(points_to_use[i][0] - Xo, 2) + 
                                           std::pow(points_to_use[i][1] - Yo, 2) / (anisotropy_factor_ * anisotropy_factor_));
    }

    Eigen::VectorXd vectorb(points_to_use.size() + 1);
    for (size_t i = 0; i < distances_to_point0.size(); ++i) {
        vectorb[i] = variogramFunction_(distances_to_point0[i], a_, C_);
    }
    vectorb[points_to_use.size()] = 1.0;

    Eigen::MatrixXd matrix_result = Eigen::MatrixXd::Map(&result_[0][0], result_.size(), result_[0].size());
    Eigen::VectorXd lamd = matrix_result.ldlt().solve(vectorb);
    lamd.conservativeResize(lamd.size() - 1);

    Eigen::VectorXd zvals_vector = Eigen::VectorXd::Map(&zvals_[0], zvals_.size());
    double zout = lamd.dot(zvals_vector);

    return zout;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<std::vector<double>>>
OrdinaryKriging::interpgrid(double grid_size) {
// Determine the range of x and y from the input points
double x_min = std::numeric_limits<double>::max();
double x_max = std::numeric_limits<double>::lowest();
double y_min = std::numeric_limits<double>::max();
double y_max = std::numeric_limits<double>::lowest();

for (const auto& point : points_) {
    x_min = std::min(x_min, point[0]);
    x_max = std::max(x_max, point[0]);
    y_min = std::min(y_min, point[1]);
    y_max = std::max(y_max, point[1]);
}

// Generate the X and Y grids based on the determined range and grid_size
std::vector<double> X, Y;
for (double x = x_min; x <= x_max; x += grid_size) {
    X.push_back(x);
}

    std::vector<std::vector<double>> Zout(X.size(), std::vector<double>(Y.size(), 0.0));
    for (size_t i = 0; i < X.size(); ++i) {
        for (size_t j = 0; j < Y.size(); ++j) {
            Zout[i][j] = SinglePoint(X[i], Y[j]);
        }
    }

    return {X, Y, Zout};
}

void OrdinaryKriging::AutoOptimize(const std::vector<std::pair<double, double>>& bounds) {
    nlopt::opt optimizer(nlopt::LD_LBFGS, 3);  // Using L-BFGS algorithm
    optimizer.set_min_objective(OrdinaryKriging::objfunctionWrapper, this);
    optimizer.set_xtol_rel(1e-4);

    std::vector<double> lbounds, ubounds;
    for (const auto& bound : bounds) {
        lbounds.push_back(bound.first);
        ubounds.push_back(bound.second);
    }
    optimizer.set_lower_bounds(lbounds);
    optimizer.set_upper_bounds(ubounds);

    std::vector<double> x = {a_, C_, nugget_};
    double minf;
    nlopt::result result = optimizer.optimize(x, minf);

    a_ = x[0];
    C_ = x[1];
    nugget_ = x[2];
}

double OrdinaryKriging::objfunction(const std::vector<double>& x, std::vector<double>& grad) {
    //TODO: Implement this objective function
}

double OrdinaryKriging::objfunctionWrapper(const std::vector<double>& x, std::vector<double>& grad, void* data) {
    return static_cast<OrdinaryKriging*>(data)->objfunction(x, grad);
}






