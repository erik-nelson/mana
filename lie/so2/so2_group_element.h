#pragma once

#include "lie/base/lie_group_element.h"

namespace mana {

struct SO2GroupElement;

// Specialization of manifold traits for SO2.
template <>
struct ManifoldTraits<SO2GroupElement> {
  using Element = SO2GroupElement;
  using Scalar = double;
  using Chart = ManifoldChart<SO2GroupElement>;
  using Geodesic = ManifoldGeodesic<SO2GroupElement>;
  using TangentVector = Eigen::Vector<Scalar, 1>;      // Angles in R.
  using EmbeddingPoint = Eigen::Matrix<Scalar, 2, 2>;  // 2D rotation matrices.
  static constexpr int Dimension = 1;
  static constexpr int EmbeddingDimension = 4;
};

class SO2GroupElement : public LieGroupElement<SO2GroupElement> {
  using Base = LieGroupElement<SO2GroupElement>;

 public:
  // Trait typedefs.
  using Scalar = typename Base::Scalar;
  using Chart = typename Base::Chart;
  using Geodesic = typename Base::Geodesic;
  using TangentVector = typename Base::TangentVector;
  using EmbeddingPoint = typename Base::EmbeddingPoint;
  using GroupElement = typename Base::GroupElement;
  using AlgebraElement = typename Base::AlgebraElement;
  using Jacobian = typename Base::Jacobian;
  static constexpr int Dimension = Base::Dimension;
  static constexpr int EmbeddingDimension = Base::EmbeddingDimension;

  // Default construct to angle=0.
  SO2GroupElement();

  // Construct from angle.
  explicit SO2GroupElement(Scalar angle_radians);

  // Syntactic sugar: Static construct from radians.
  static SO2GroupElement FromRadians(Scalar angle_radians);

  // Syntactic sugar: Static construct from degrees.
  static SO2GroupElement FromDegrees(Scalar angle_degrees);

  // Retrieve underlying angle.
  Scalar AngleRadians() const;
  Scalar AngleDegrees() const;

  // Get the rotation matrix for this element (returns a 2x2 matrix).
  EmbeddingPoint AsMatrix() const;

  // Implement `GroupElement<>` interface.
  static SO2GroupElement Identity();
  SO2GroupElement Inverse() const;
  SO2GroupElement Compose(const SO2GroupElement& rhs) const;

  // Implement `ManifoldElement<>` interface.
  static SO2GroupElement FromPointImpl(const EmbeddingPoint& point);
  static EmbeddingPoint ProjectImpl(const EmbeddingPoint& point);
  static bool IsValidImpl(const EmbeddingPoint& point, Scalar tolerance);
  EmbeddingPoint PointImpl() const;

  // Implement `LieGroupElement<>` interface.
  TangentVector LogImpl() const;
  static SO2GroupElement ExpImpl(const TangentVector& coordinate);
  Jacobian AdjointImpl() const;

 private:
  // Private ctor assumes c^2 + s^2 = 1.
  SO2GroupElement(Scalar cos_theta, Scalar sin_theta);

  // Underlying storage is sin(theta), cos(theta), where theta is the angle of
  // rotation. This is a happy medium between storing the full rotation matrix
  // (large storage overhead), and storing the angle (costly to perform sin/cos
  // when casting back to matrix).
  Scalar cos_theta_, sin_theta_;
};

}  // namespace mana