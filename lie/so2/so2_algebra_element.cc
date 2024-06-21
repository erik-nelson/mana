
#include "lie/so2/so2_algebra_element.h"

#include "lie/so2/so2_group_element.h"
#include "utils/angles.h"

namespace mana {

so2AlgebraElement::so2AlgebraElement(Scalar angle_radians)
    : so2AlgebraElement(Vector::Constant(angle_radians)) {}

/*static*/ so2AlgebraElement so2AlgebraElement::FromRadians(
    Scalar angle_radians) {
  return so2AlgebraElement(angle_radians);
}

/*static*/ so2AlgebraElement so2AlgebraElement::FromDegrees(
    Scalar angle_degrees) {
  return so2AlgebraElement(DegToRad(angle_degrees));
}

so2AlgebraElement::Scalar so2AlgebraElement::AngleRadians() const {
  return coordinates_(0);
}

so2AlgebraElement::Scalar so2AlgebraElement::AngleDegrees() const {
  return RadToDeg(AngleRadians());
}

so2AlgebraElement::so2AlgebraElement(Vector coordinates)
    : coordinates_(std::move(coordinates)) {}

so2AlgebraElement::so2AlgebraElement(const Matrix matrix) {
  coordinates_(0) = matrix(1, 0);
}

so2AlgebraElement::Matrix so2AlgebraElement::AsMatrixImpl() const {
  const Scalar angle_radians = AngleRadians();
  Matrix matrix;
  matrix << 0, -angle_radians, angle_radians, 0;
  return matrix;
}

}  // namespace mana