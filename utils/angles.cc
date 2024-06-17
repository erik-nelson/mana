#include "utils/angles.h"

#include <cmath>

namespace mana {

double Normalized(double radians) {
  radians = fmod(radians, 2.0 * M_PI);
  return (radians < 0) ? (radians + 2.0 * M_PI) : radians;
}

double RadToDeg(double radians) { return radians * 180.0 / M_PI; }

double DegToRad(double degrees) { return degrees / 180.0 * M_PI; }

double AngleBetween(double beg, double end) {
  // Normalize angles to [0, 2*pi).
  beg = Normalized(beg);
  end = Normalized(end);

  // Compute the difference.
  double delta = end - beg;

  // Adjust for wrap-around.
  if (delta > M_PI) {
    delta -= 2.0 * M_PI;
  } else if (delta < -M_PI) {
    delta += 2.0 * M_PI;
  }
  return delta;
}

double InterpolateAngles(double beg, double end, double fraction) {
  return Normalized(beg + fraction * AngleBetween(beg, end));
}

void CartesianToSpherical(double x, double y, double z, double& rho,
                          double& theta, double& phi) {
  rho = std::sqrt(x * x + y * y + z * z);
  theta = std::atan2(y, x);
  phi = std::acos(z / rho);
}

void SphericalToCartesian(double rho, double theta, double phi, double& x,
                          double& y, double& z) {
  const double sin_phi = std::sin(phi);
  const double cos_phi = std::cos(phi);
  const double sin_theta = std::sin(theta);
  const double cos_theta = std::cos(theta);
  x = rho * sin_phi * cos_theta;
  y = rho * sin_phi * sin_theta;
  z = rho * cos_phi;
}

}  // namespace mana