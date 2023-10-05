#ifndef ORDINARY_KRIGING_H
#define ORDINARY_KRIGING_H

#include <vector>
#include <string>
#include <functional>
#include <Eigen/Dense>  // For matrix operations

class OrdinaryKriging {
public:
    OrdinaryKriging(const std::vector<std::vector<double>>& points, 
                    const std::vector<double>& zvals, 
                    const std::string& variogram = "gaussian");

    void setupVariogram();
    void ManualParamSet(double C, double a, double nugget, double anisotropy_factor);
    Eigen::MatrixXd MatrixSetup();
    double SinglePoint(double Xo, double Yo, const std::vector<std::vector<double>>& training_points = {});
    std::tuple<std::vector<double>, std::vector<double>, std::vector<std::vector<double>>>
    interpgrid(double grid_size);
    void AutoOptimize(const std::vector<std::pair<double, double>>& bounds);
    double Predict(const std::vector<double>& point);

    std::vector<std::vector<double>> points_;
    std::vector<double> zvals_;
    std::vector<double> zvals_org_;
    std::string variogram_;
    double anisotropy_;
    std::function<double(double, double, double)> variogramFunction_;
    double a_;
    double anisotropy_factor_;
    double C_;
    double nugget_;

private:
    static double objfunctionWrapper(const std::vector<double>& x, std::vector<double>& grad, void* data);
    double objfunction(const std::vector<double>& x);
    void computeGradientCentralDifference(const std::vector<double>& x, std::vector<double>& grad);
    double gradientCallback(const std::vector<double> &x, std::vector<double> &grad, void *data);

    std::vector<std::vector<double>> result_;
};

#endif // ORDINARY_KRIGING_H





