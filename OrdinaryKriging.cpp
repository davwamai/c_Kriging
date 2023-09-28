#include "./headers/OrdinaryKriging.h"
#include <cmath>
#include <numeric>
#include <iostream>
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
        std::cerr << "No valid Variogram model selected" << std::endl;
        exit(1);
    }
}

void OrdinaryKriging::ManualParamSet(double C, double a, double nugget, double anisotropy_factor) {
    a_ = a;
    anisotropy_factor_ = anisotropy_factor;
    C_ = C;
    nugget_ = nugget;
}

Eigen::MatrixXd OrdinaryKriging::MatrixSetup() {
    size_t n = points_.size();
    Eigen::MatrixXd A(n, n);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            if (i == j) {
                A(i, j) = variogramFunction_(0.0, a_, C_) + nugget_;  // Diagonal element (distance = 0), if I read the paper correctly...
            } else {
                double dx = points_[i][0] - points_[j][0];
                double dy = points_[i][1] - points_[j][1];
                double distance = std::sqrt(dx * dx + dy * dy);
                A(i, j) = variogramFunction_(distance, a_, C_);
            }
        }
    }

    return A;
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

//this shit is fuuuuuuucked lmao
void OrdinaryKriging::AutoOptimize(const std::vector<std::pair<double, double>>& bounds) {
    std::cout<< "Running Auto Optimize..." << std::endl;

    std::cout << "LD_LBFGS Optimizing..." << std::endl;
    nlopt::opt optimizer(nlopt::LD_LBFGS, 4);  // Using L-BFGS algorithm
    optimizer.set_min_objective(OrdinaryKriging::objfunctionWrapper, this);
    optimizer.set_xtol_rel(1e-4);

    std::vector<double> lbounds, ubounds;
    std::cout << "Unpacking Bounds..." << std::endl;
    for (const auto& bound : bounds) {
        lbounds.push_back(bound.first);
        ubounds.push_back(bound.second);
    }
    
    optimizer.set_lower_bounds(lbounds);
    optimizer.set_upper_bounds(ubounds);

    std::cout << "Unoptimized Parameters: " << std::endl;
    std::cout << "a: " << a_ << std::endl;
    std::cout << "C: " << C_ << std::endl;
    std::cout << "nugget: " << nugget_ << std::endl;
    std::cout << "anisotropy_factor: " << anisotropy_factor_ << std::endl;
    std::cout << "\n" << std::endl;

    std::vector<double> x = {a_, C_, nugget_, anisotropy_factor_};
    double minf;
    nlopt::result result = optimizer.optimize(x, minf); //failing here  

    a_ = x[0];
    C_ = x[1];
    nugget_ = x[2];
    anisotropy_factor_ = x[3];

    std::cout << "Optimized Parameters: " << std::endl;
    std::cout << "a: " << a_ << std::endl;
    std::cout << "C: " << C_ << std::endl;
    std::cout << "nugget: " << nugget_ << std::endl;
    std::cout << "anisotropy_factor: " << anisotropy_factor_ << std::endl;
    std::cout << "\n" << std::endl;

}

double OrdinaryKriging::objfunction(const std::vector<double>& x) {
        double sill = x[0];
        double range = x[1];
        double nugget = x[2];
        double anisotropy_factor = x[3];
        
        std::vector<double> predictions;
        std::vector<double> actuals;
        double sumZ = std::accumulate(zvals_.begin(), zvals_.end(), 0.0);
        double meanZ = sumZ / zvals_.size();

        // std::cout << "Performing LLO..." << std::endl;
        for (size_t i = 0; i < points_.size(); ++i) {
            // Create a new dataset excluding the i-th point
            std::vector<std::vector<double>> newPoints = points_;
            newPoints.erase(newPoints.begin() + i);
            std::vector<double> newZ = zvals_;
            newZ.erase(newZ.begin() + i);

            // Create a new OrdinaryKriging model with the new dataset
            OrdinaryKriging model(newPoints, newZ, "gaussian");
            model.ManualParamSet(x[1], x[0], x[2], x[3]);

            // Predict the value at the i-th point
            double prediction = model.Predict(points_[i]);
            predictions.push_back(prediction);
            actuals.push_back(zvals_[i]);
        }

        // Calculate R-squared value
        double numerator = 0.0;
        double denominator = 0.0;
        for (size_t i = 0; i < points_.size(); ++i) {
            numerator += std::pow(actuals[i] - predictions[i], 2);
            denominator += std::pow(actuals[i] - meanZ, 2);
        }
        double rSquared = 1.0 - (numerator / denominator);  

        // std::cout << "R squared: " << rSquared << std::endl;

        // Return 1 - R-squared to minimize it
        return 1.0 - rSquared;
}

//Nned a wrapper for nlopt, expects a double ret from OK class
//L-BFGS is gradient based but grad param can still be ignored, pass a gradient if using a derivative based method
//the 'real' wrapper is the Gradient Callback
double OrdinaryKriging::objfunctionWrapper(const std::vector<double>& x, std::vector<double>& grad, void* data) {
    OrdinaryKriging* kriging = reinterpret_cast<OrdinaryKriging*>(data);
    kriging->computeGradientCentralDifference(x, grad);
    return kriging->objfunction(x);
}

void OrdinaryKriging::computeGradientCentralDifference(const std::vector<double>& x, std::vector<double>& grad) {
    double delta = 1e-6; // small perturbation value
    double forwardObj, backwardObj;

    for (size_t i = 0; i < x.size(); i++) {
        std::vector<double> forwardX = x;
        std::vector<double> backwardX = x;

        forwardX[i] += delta;
        backwardX[i] -= delta;

        forwardObj = objfunction(forwardX);
        backwardObj = objfunction(backwardX);

        grad[i] = (forwardObj - backwardObj) / (2 * delta);
    }
}

double OrdinaryKriging::gradientCallback(const std::vector<double> &x, std::vector<double> &grad, void *data) {
    computeGradientCentralDifference(x, grad);
    return objfunctionWrapper(x, grad, data);
}

double OrdinaryKriging::Predict(const std::vector<double>& point) {
    // Ensures the point has the correct dimension
    if (point.size() != 2) {
        throw std::invalid_argument("Point dimension should be 2 for 2D kriging.");
    }

    // Calculate distances between the prediction point and all training points
    Eigen::VectorXd distances(points_.size());
    for (size_t i = 0; i < points_.size(); i++) {
        double dx = point[0] - points_[i][0];
        double dy = point[1] - points_[i][1];
        distances(i) = std::sqrt(dx * dx + dy * dy);
    }

    Eigen::MatrixXd A = MatrixSetup();
    Eigen::VectorXd b(points_.size());
    for (size_t i = 0; i < points_.size(); i++) {
        b(i) = variogramFunction_(distances(i), a_, C_);
    }

    // Solve for weights
    Eigen::VectorXd weights = A.colPivHouseholderQr().solve(b); //No clue how this works but Eigen said it would be better than solving in static

    // Compute prediction
    double prediction = 0.0;
    for (size_t i = 0; i < zvals_.size(); i++) {
        prediction += weights(i) * zvals_[i];
    }

    return prediction;
}









