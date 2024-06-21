#pragma once

#include "lie/base/lie_algebra_element.h"
#include "lie/so2/so2_group_element.h"
#include "so2_group_element.h"

namespace mana {

struct so2AlgebraElement;

// Specialization of algebra traits for so2.
template <>
struct AlgebraTraits<so2AlgebraElement> {
  using Element = so2AlgebraElement;
  using Scalar = double;
  using Vector = Eigen::Vector<Scalar, 1>;
  static constexpr int Dimension = 1;
};

// Specialization of Lie algebra traits for so2.
template <>
struct LieAlgebraTraits<so2AlgebraElement> {
  using AlgebraElement = so2AlgebraElement;
  using GroupElement = SO2GroupElement;
  using Matrix =
      Eigen::Matrix<typename AlgebraTraits<so2AlgebraElement>::Scalar, 2, 2>;
};

class so2AlgebraElement : public LieAlgebraElement<so2AlgebraElement> {
  using Base = LieAlgebraElement<so2AlgebraElement>;

 public:
  // Trait typedefs.
  using Scalar = typename Base::Scalar;
  using Vector = typename Base::Vector;
  using Matrix = typename Base::Matrix;
  using GroupElement = typename Base::GroupElement;
  using AlgebraElement = typename Base::AlgebraElement;
  static constexpr int Dimension = Base::Dimension;

  // Construct from angle.
  explicit so2AlgebraElement(Scalar angle_radians);

  // Syntactic sugar: Static construct from radians.
  static so2AlgebraElement FromRadians(Scalar angle_radians);

  // Syntactic sugar: Static construct from degrees.
  static so2AlgebraElement FromDegrees(Scalar angle_degrees);

  // Retrieve underlying angle.
  Scalar AngleRadians() const;
  Scalar AngleDegrees() const;

  // --------------------------------------------------------------------------
  // Implement `AlgebraElement<>` interface.
  explicit so2AlgebraElement(Vector coordinates);

  // Implement `LieAlgebraElement<>` interface.
  explicit so2AlgebraElement(const Matrix& matrix);
  Matrix AsMatrixImpl() const;

  /* The following methods are inherited from `LieAlgebraElement<>`:
   *
   *  Vector& Coordinates();
   *
   *  const Vector& Coordinates() const;
   *
   *  static so2AlgebraElement Hat(Vector coordinate);
   *
   *  const Vector& Vee() const;
   *
   *  so2AlgebraElement Compose(const so2AlgebraElement& rhs) const;
   *
   *  so2AlgebraElement operator+(const so2AlgebraElement& rhs) const;
   *
   *  so2AlgebraElement& operator+=(const so2AlgebraElement& rhs);
   *
   *  so2AlgebraElement operator-(const so2AlgebraElement& rhs) const;
   *
   *  so2AlgebraElement& operator-=(const so2AlgebraElement& rhs);
   *
   *  so2AlgebraElement operator*(const so2AlgebraElement& rhs) const;
   *
   *  so2AlgebraElement& operator*=(const so2AlgebraElement& rhs);
   *
   *  so2AlgebraElement operator-() const;
   *
   *  so2AlgebraElement operator*(Scalar rhs) const;
   *
   *  so2AlgebraElement& operator*=(Scalar rhs);
   *
   *  so2AlgebraElement operator/(Scalar rhs) const;
   *
   *  so2AlgebraElement& operator/=(Scalar rhs);
   *
   *  bool operator==(const so2AlgebraElement& rhs) const;
   *
   *  bool operator!=(const so2AlgebraElement& rhs) const;
   *
   *  so2AlgebraElement operator*(Scalar lhs, const so2AlgebraElement& rhs);
   *
   *  so2AlgebraElement Bracket(const so2AlgebraElement& rhs) const;
   *
   *  static so2AlgebraElement log(const SO2GroupElement& group);
   *
   *  SO2GroupElement exp() const;
   */

 private:
  Vector coordinates_;
};

}  // namespace mana