#include "./headers/UniversalKriging.h"
#include <limits>

UniversalKriging::UniversalKriging(const std::vector<std::vector<double>>& points, 
                                   const std::vector<double>& zvals, 
                                   const std::string& variogram,
                                   const std::string& trend_function)
    : OrdinaryKriging(points, zvals, variogram), trend_function_(trend_function) {
    calculateTrendCoefficients();
}

void UniversalKriging::calculateTrendCoefficients() {
    //TODO: Implement this method
}

double UniversalKriging::SinglePointWithTrend(double Xo, double Yo, const std::vector<std::vector<double>>& training_points) {
    double zout = SinglePoint(Xo, Yo, training_points);
    if (trend_function_ == "first") {
        zout += trend_coeffs_[0] + trend_coeffs_[1] * Xo + trend_coeffs_[2] * Yo;
    }
    //TODO: Add other trend functions

    return zout;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<std::vector<double>>>
UniversalKriging::interpgridWithTrend(double grid_size) {
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
for (double y = y_min; y <= y_max; y += grid_size) {
    Y.push_back(y);
}

    std::vector<std::vector<double>> Zout(X.size(), std::vector<double>(Y.size(), 0.0));
    for (size_t i = 0; i < X.size(); ++i) {
        for (size_t j = 0; j < Y.size(); ++j) {
            Zout[i][j] = SinglePointWithTrend(X[i], Y[j]);
        }
    }

    return {X, Y, Zout};
}

void UniversalKriging::AutoOptimizeWithTrend(const std::vector<std::pair<double, double>>& bounds) {
    nlopt::opt optimizer(nlopt::LD_LBFGS, 3);  // Using L-BFGS algorithm
    optimizer.set_min_objective(UniversalKriging::objfunctionWithTrendWrapper, this);
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

double UniversalKriging::objfunctionWithTrendWrapper(const std::vector<double>& x, std::vector<double>& grad, void* data) {
    return static_cast<UniversalKriging*>(data)->objfunctionWithTrend(x, grad);
}

double UniversalKriging::objfunctionWithTrend(const std::vector<double>& x, std::vector<double>& grad) {
    //TODO: Implement this objective function with trend
}

