#ifndef VARIOGRAMS_H
#define VARIOGRAMS_H
#include <cmath>

namespace Variogram {

double Gaussian(double h, double a, double C);
double Spherical(double h, double a, double C);
double Exponential(double h, double a, double C);

//TODO: Add other variogram models

}  // namespace Variogram
#endif // VARIOGRAMS_H
