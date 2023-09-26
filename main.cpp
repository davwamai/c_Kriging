#include <iostream>
#include "./headers/OrdinaryKriging.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
#include <cmath>


int main() {
    // // Read data from CSV
    // const std::string filename = "surface_roughness.csv";
    // std::ifstream file(filename);
    //     std::string line;
    //     std::vector<std::vector<double>> points;
    //     std::vector<double> zvals;

    //     while (std::getline(file, line)) {
    //         std::stringstream ss(line);
    //         std::string cell;
    //         std::vector<double> point(2);  // x and y values
    //         double z;  // target value
    //         int colIndex = 0;

    //         while (std::getline(ss, cell, ',')) {
    //             if (colIndex == 0) {
    //                 point[0] = std::stod(cell);
    //             } else if (colIndex == 1) {
    //                 point[1] = std::stod(cell);
    //             } else if (colIndex == 4) {
    //                 z = std::stod(cell);
    //             }
    //             colIndex++;
    //         }
    //         points.push_back(point);
    //         zvals.push_back(z);
    //     }

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
    std::vector<std::pair<double, double>> bounds = {{0.1, 10}, {0.1, 10}, {0.1, 5}}; // Find some comfortable bounds
    krigingModel.ManualParamSet(1.0, 1.0, 1.0, 1.0); // Set initial parameters
    krigingModel.AutoOptimize(bounds);

    // Predict targets for the points in the dataset
    std::vector<double> predictedZvals;
    for (const auto& point : points) {
        double predictedZ = krigingModel.Predict(point);
        predictedZvals.push_back(predictedZ);
    }

    // Compute MSE
    double mse = 0.0;
    for (size_t i = 0; i < zi.size(); i++) {
        mse += std::pow(zi[i] - predictedZvals[i], 2);
        std::cout << "Actual: " << zi[i] << " Predicted: " << predictedZvals[i] << std::endl;
        std::cout << "MSE at point " << i << ": " << std::pow(zi[i] - predictedZvals[i], 2) << std::endl << std::endl; // Print MSE at each point
    }
    mse /= zi.size();

    std::cout << "Total MSE: " << mse << std::endl;

    return 0;
}



