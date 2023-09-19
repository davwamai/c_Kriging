#pragma once
#include "OrdinaryKriging.h"
#include <nlopt.hpp>


class UniversalKriging : public OrdinaryKriging {
public:
    UniversalKriging(const std::vector<std::vector<double>>& points, 
                     const std::vector<double>& zvals, 
                     const std::string& variogram = "gaussian",
                     const std::string& trend_function = "first");
    //TODO

private:
    std::string trend_function_;
    std::vector<double> trend_coeffs_;
    void calculateTrendCoefficients();

public:
    double SinglePointWithTrend(double Xo, double Yo, const std::vector<std::vector<double>>& training_points = {});
    //TODO

public:
    std::tuple<std::vector<double>, std::vector<double>, std::vector<std::vector<double>>>
    interpgridWithTrend(double grid_size);

public:
    void AutoOptimizeWithTrend(const std::vector<std::pair<double, double>>& bounds);

private:
    static double objfunctionWithTrendWrapper(const std::vector<double>& x, std::vector<double>& grad, void* data);
    double objfunctionWithTrend(const std::vector<double>& x, std::vector<double>& grad);
};


