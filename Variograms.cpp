#include "./headers/Variograms.h"
//Will more than likely abandon this and implement all the Varios in OK.cpp
namespace Variogram {

double Gaussian(double h, double a, double C) {
    return C * (1 - std::exp(-1 * ((h / a) * (h / a))));
}

double Spherical(double h, double a, double C) {
    //TODO
}

double Exponential(double h, double a, double C) {
    //TODO
}

//TODO: Add other variogram models

}  // namespace Variogram
