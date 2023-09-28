#include <iostream>
#include "./headers/OrdinaryKriging.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
#include <cmath>
#include <random>

double variance(const std::vector<double>& zvals) {
    double mean = std::accumulate(zvals.begin(), zvals.end(), 0.0) / zvals.size();
    double sq_sum = std::inner_product(zvals.begin(), zvals.end(), zvals.begin(), 0.0);
    double var = sq_sum / zvals.size() - mean * mean;
    return var;
}

double maxPairwiseDistance(const std::vector<std::vector<double>>& points) {
    double maxDist = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            double dist = 0.0;
            for (size_t k = 0; k < points[i].size(); ++k) {
                dist += (points[i][k] - points[j][k]) * (points[i][k] - points[j][k]);
            }
            maxDist = std::max(maxDist, std::sqrt(dist));
        }
    }
    return maxDist;
}

int main() {
    
    int xi_col = 0;
    int yi_col = 1;
    int zi_col = 4;
    std::ifstream file("surface_roughness.csv");
    std::vector<std::vector<double>> loadeddata;
    std::string line;
    std::getline(file, line); // skip header
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream ss(line);
        std::string value;
        while (std::getline(ss, value, ',')) {
            row.push_back(std::stod(value));
        }
        loadeddata.push_back(row);
    }
    std::vector<std::vector<double>> points;
    std::vector<double> zi;
    for (const auto& row : loadeddata) {
        points.push_back({row[xi_col], row[yi_col]});
        zi.push_back(row[zi_col]);
    }

    // Instantiate model
    OrdinaryKriging krigingModel(points, zi, "gaussian");

    // Optimize models parameters
    std::vector<std::pair<double, double>> bounds = {{15, 30}, {1900, 2000}, {0.001, 10}, {0.1, 10}}; // Find some comfortable bounds

    std::vector<double> InitialParams = {
        variance(zi),
        maxPairwiseDistance(points) / 2,
        0.001,
        1.0
    };

    double a = InitialParams[0];
    double C = InitialParams[1];
    double nugget = InitialParams[2];
    double anisotropy_factor = InitialParams[3];

    krigingModel.ManualParamSet(C, a, nugget, anisotropy_factor); // Set initial parameters, LMAO

    krigingModel.AutoOptimize(bounds);

    // Predict targets for the points in the dataset
    std::vector<double> predictedZvals;
    for (const auto& point : points) {
        double predictedZ = krigingModel.Predict(point);
        predictedZvals.push_back(predictedZ);
    }

    //calculates r squared value for the model
    double sumZ = std::accumulate(zi.begin(), zi.end(), 0.0);
    double meanZ = sumZ / zi.size();
    double ss_tot = 0.0;
    double ss_res = 0.0;
    for (size_t i = 0; i < zi.size(); i++) {
        ss_tot += std::pow(zi[i] - meanZ, 2);
        ss_res += std::pow(zi[i] - predictedZvals[i], 2);
    }
    double r_squared = 1 - (ss_res / ss_tot);
    std::cout << "ss_tot: " << ss_tot << std::endl;
    std::cout << "ss_res: " << ss_res << std::endl;

    // Compute MSE
    double mse = 0.0;
    for (size_t i = 0; i < zi.size(); i++) {
        mse += std::pow(zi[i] - predictedZvals[i], 2);
    }
    mse /= zi.size();

    std::cout << "Overall MSE: " << mse << std::endl;
    std::cout << "Overall r^2: " << r_squared << std::endl;

    return 0;
}



