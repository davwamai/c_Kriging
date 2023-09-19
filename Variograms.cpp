#include "./headers/Variograms.h"

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
