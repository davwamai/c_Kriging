#pragma once
#include <vector>
#include <string>
#include <functional>
#include <Eigen/Dense>  // For matrix operations


class OrdinaryKriging {
public:
    OrdinaryKriging(const std::vector<std::vector<double>>& points, 
                    const std::vector<double>& zvals, 
                    const std::string& variogram = "gaussian");

public:
    std::vector<std::vector<double>> points_;
    std::vector<double> zvals_;
    std::vector<double> zvals_org_;
    std::string variogram_;
    double anisotropy_;
    std::function<double(double, double, double)> variogramFunction_;

    void setupVariogram();


public:
    void ManualParamSet(double C, double a, double nugget, double anisotropy_factor);

public:
    double a_;
    double anisotropy_factor_;
    double C_;
    double nugget_;

public:
    Eigen::MatrixXd MatrixSetup();

private:
    std::vector<std::vector<double>> result_;

public:
    double SinglePoint(double Xo, double Yo, const std::vector<std::vector<double>>& training_points = {});

public:
    std::tuple<std::vector<double>, std::vector<double>, std::vector<std::vector<double>>>
    interpgrid(double grid_size);

public:
    void AutoOptimize(const std::vector<std::pair<double, double>>& bounds);

private:
    static double objfunctionWrapper(const std::vector<double>& x, std::vector<double>& grad, void* data);
    double objfunction(const std::vector<double>& x);
    void computeGradientCentralDifference(const std::vector<double>& x, std::vector<double>& grad);
    double gradientCallback(const std::vector<double> &x, std::vector<double> &grad, void *data);

public: 
    double Predict(const std::vector<double>& point);

};




