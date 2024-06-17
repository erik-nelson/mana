#pragma once

namespace mana {

// Normalize the provided angle to lie in [0, 2*pi).
double Normalized(double radians);

// Convert radians to degrees.
double RadToDeg(double radians);

// Convert degrees to radians.
double DegToRad(double degrees);

// Compute the smallest angle between two other angles, accounting for
// wrap-around. Values in [0, pi).
double AngleBetween(double beg, double end);

// Function to interpolate between two angles, accounting for wrap-around.
// Provided fraction should be from [0, 1].
double InterpolateAngles(double beg, double end, double fraction);

// Convert a 3D point from cartesian to spherical coordinates.
void CartesianToSpherical(double x, double y, double z, double& rho,
                          double& theta, double& phi);

// Convert a 3D point from spherical to cartesian coordinates.
void SphericalToCartesian(double rho, double theta, double phi, double& x,
                          double& y, double& z);

}  // namespace mana