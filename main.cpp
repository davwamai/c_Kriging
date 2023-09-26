#include <iostream>
#include "./headers/OrdinaryKriging.h"
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <tuple>
#include <cmath>


int main() {
    // Read data from CSV
    const std::string filename = "path_to_your_dataset.csv";
    std::ifstream file(filename);
        std::string line;
        std::vector<std::vector<double>> points;
        std::vector<double> zvals;

        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string cell;
            std::vector<double> point(2);  // x and y values
            double z;  // target value
            int colIndex = 0;

            while (std::getline(ss, cell, ',')) {
                if (colIndex == 0) {
                    point[0] = std::stod(cell);
                } else if (colIndex == 1) {
                    point[1] = std::stod(cell);
                } else if (colIndex == 4) {
                    z = std::stod(cell);
                }
                colIndex++;
            }
            points.push_back(point);
            zvals.push_back(z);
        }


    // Instantiate model
    OrdinaryKriging krigingModel(points, zvals, "gaussian");

    // Optimize models parameters
    std::vector<std::pair<double, double>> bounds = {{0.1, 10}, {0.1, 10}, {0, 5}}; // Find some comfortable bounds
    krigingModel.AutoOptimize(bounds);

    // Predict targets for the points in the dataset
    std::vector<double> predictedZvals;
    for (const auto& point : points) {
        double predictedZ = krigingModel.Predict(point);
        predictedZvals.push_back(predictedZ);
    }

    // Compute MSE
    double mse = 0.0;
    for (size_t i = 0; i < zvals.size(); i++) {
        mse += std::pow(zvals[i] - predictedZvals[i], 2);
    }
    mse /= zvals.size();

    std::cout << "MSE: " << mse << std::endl;

    return 0;
}


