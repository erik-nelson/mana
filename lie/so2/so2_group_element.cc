#include "lie/so2/so2_group_element.h"

#include "so2_group_element.h"
#include "utils/angles.h"

namespace mana {

SO2GroupElement::SO2GroupElement() : SO2GroupElement(1, 0) {}

SO2GroupElement::SO2GroupElement(Scalar angle_radians)
    : SO2GroupElement(std::cos(angle_radians), std::sin(angle_radians)) {}

/*static*/ SO2GroupElement SO2GroupElement::FromRadians(Scalar angle_radians) {
  return SO2GroupElement(angle_radians);
}

/*static*/ SO2GroupElement SO2GroupElement::FromDegrees(Scalar angle_degrees) {
  return SO2GroupElement(DegToRad(angle_degrees));
}

SO2GroupElement::Scalar SO2GroupElement::AngleRadians() const {
  return std::acos(cos_theta_);
}

SO2GroupElement::Scalar SO2GroupElement::AngleDegrees() const {
  return RadToDeg(AngleRadians());
}

SO2GroupElement::EmbeddingPoint SO2GroupElement::AsMatrix() const {
  return Point();
}

/*static*/ SO2GroupElement SO2GroupElement::IdentityImpl() {
  return SO2GroupElement();
}

SO2GroupElement SO2GroupElement::InverseImpl() const {
  return SO2GroupElement(cos_theta_, -sin_theta_);
}

SO2GroupElement SO2GroupElement::ComposeImpl(const SO2GroupElement& rhs) const {
  // Complex multiplication: (a + bi) * (c + di) = (ac - bd) + (ad + bc)i.
  // (where a = cos(lhs.t), b = sin(lhs.t), c = cos(rhs.t), d = sin(rhs.t)).
  Scalar cos_theta = cos_theta_ * rhs.cos_theta_ - sin_theta_ * rhs.sin_theta_;
  Scalar sin_theta = cos_theta_ * rhs.sin_theta_ - sin_theta_ * rhs.cos_theta_;
  const Scalar magnitude = cos_theta * cos_theta + sin_theta * sin_theta;
  if (std::abs(magnitude - 1) >= Constants<Scalar>::kEpsilon) {
    cos_theta /= magnitude;
    sin_theta /= magnitude;
  }
  return SO2GroupElement(cos_theta, sin_theta);
}

/*static*/ SO2GroupElement SO2GroupElement::FromPointImpl(
    const EmbeddingPoint& point) {
  assert(IsValid(point));
  return SO2GroupElement(point(0, 0), point(1, 0));
}

/*static*/ SO2GroupElement::EmbeddingPoint SO2GroupElement::ProjectImpl(
    const EmbeddingPoint& point) {
  Eigen::JacobiSVD<EmbeddingPoint> svd(
      point, Eigen::ComputeFullU | Eigen::ComputeFullV);
  EmbeddingPoint U = svd.matrixU();
  EmbeddingPoint V = svd.matrixV();
  if (U.determinant() * V.determinant() < 0) {
    U.col(1) *= -1;
  }
  return (U * V.transpose()).eval();
}

/*static*/ bool SO2GroupElement::IsValidImpl(const EmbeddingPoint& point,
                                             Scalar tolerance) {
  // Matrix must of the form [a, -b; b a], where a^2+b^2 = 1. This is equivalent
  // to: r11 ~= r22 r12 ~= -r21 r.col(0) * r.col(1) ~= 0
  if ((point(0, 0) - point(1, 1)) >= tolerance) return false;
  if ((point(0, 1) + point(1, 0)) >= tolerance) return false;
  if (point.col(0).dot(point.col(1)) >= tolerance) return false;
  return true;
}

SO2GroupElement::EmbeddingPoint SO2GroupElement::PointImpl() const {
  // Points in the embedding space are 2x2 rotation matrices.
  Eigen::Matrix<Scalar, 2, 2> point;
  point << cos_theta_, -sin_theta_, sin_theta_, cos_theta_;
  return point;
}

SO2GroupElement::TangentVector SO2GroupElement::LogImpl() const {
  return TangentVector(AngleRadians());
}

/*static*/ SO2GroupElement SO2GroupElement::ExpImpl(
    const TangentVector& coordinate) {
  const Scalar angle = coordinate(0);
  return SO2GroupElement(angle);
}

SO2GroupElement::Jacobian SO2GroupElement::AdjointImpl() const {
  return Jacobian::Constant(1);
}

SO2GroupElement::SO2GroupElement(Scalar cos_theta, Scalar sin_theta)
    : cos_theta_(cos_theta), sin_theta_(sin_theta) {
  assert(std::abs(cos_theta_ * cos_theta_ + sin_theta_ * sin_theta - 1) <
         Constants<Scalar>::kEpsilon);
}

}  // namespace mana