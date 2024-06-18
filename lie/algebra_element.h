#pragma once

#include <Eigen/Dense>
#include <utility>

#include "utils/crtp.h"

namespace mana {

// Traits template for an algebra.
template <typename Derived>
struct AlgebraTraits {};

// Base CRTP class for an element of an algebra. An algebra is a vector space
// with a bilinear product defined (i.e. a product on two vectors from the
// vector space that produces a third vector).
// See https://en.wikipedia.org/wiki/Algebra_over_a_field.
//
// Note that elements of vector spaces are not generically coordinate vectors.
// For instance, polynomial functions with coefficients from a field form a
// vector space, whose elements are polynomial functions (not vectors, in the
// Eigen::VectorXd sense). This type is meant to represent the former, but its
// underlying storage is the latter. That is, we only store the relevant
// coefficients of the coordinate vector.
//
// Derived group algebra element types should inherit from this like so:
//
//   class SO3AlgebraElement : public AlgebraElement<SO3AlgebraElement> {
//     ...
//   };
//
// Derived classes must implement these methods.
// - explicit Element(Vector coordinates); // (constructor)
// - Element ComposeImpl(const Element& rhs) const;
template <typename Derived>
class AlgebraElement : public Crtp<Derived> {
 public:
  // An element of the algebra (identical to `Derived`).
  using Element = typename AlgebraTraits<Derived>::Element;
  // The scalar type used to represent coordinates.
  using Scalar = typename AlgebraTraits<Derived>::Scalar;
  // A coordinate vector for this algebra.
  using Vector = typename AlgebraTraits<Derived>::Vector;
  // The dimension of the algebra's vector space.
  static constexpr int Dimension = AlgebraTraits<Derived>::Dimension;

  // Get the vector this algebra element corresponds to.
  Vector& Coordinates() { return coordinates_; }
  const Vector& Coordinates() const { return coordinates_; }

  // The hat operator builds an algebra element from a coordinate vector.
  static Element Hat(Vector coordinate) {
    return Element(std::move(coordinate));
  }

  // The vee operator builds a coordinate vector from this algebra element. This
  // is synonmous with `this->Coordinates()`.
  const Vector& Vee() const { return Coordinates(); }

  // The bilinear product of this algebra.
  Element Compose(const Element& rhs) const {
    return Crtp<Derived>::get().ComposeImpl(rhs);
  }

  // Vector addition: a + b.
  Element operator+(const Element& rhs) const {
    return Element(coordinates_ + rhs.coordinates_);
  }

  // Vector in-place addition: a += b.
  Element& operator+=(const Element& rhs) {
    coordinates_ += rhs.coordinates_;
    return Crtp<Derived>::get();
  }

  // Vector subtraction: a - b.
  Element operator-(const Element& rhs) const {
    return Element(coordinates_ - rhs.coordinates_);
  }

  // Vector in-place subtraction: a -= b.
  Element& operator-=(const Element& rhs) {
    coordinates_ -= rhs.coordinates_;
    return Crtp<Derived>::get();
  }

  // Vector product: a * b.
  Element operator*(const Element& rhs) const { return Compose(rhs); }

  // Vector in-place product: a *= b.
  Element& operator*=(const Element& rhs) {
    Element result = Compose(rhs);
    coordinates_ = std::move(result.coordinates_);
    return Crtp<Derived>::get();
  }

  // Vector negation: -a.
  Element operator-() const { return Element(-coordinates_); }

  // Vector-scalar multiplication: v * s.
  Element operator*(Scalar rhs) const { return Element(coordinates_ * rhs); }

  // Vector-scalar in-place multiplication: v *= s.
  Element& operator*=(Scalar rhs) {
    coordinates_ *= rhs;
    return Crtp<Derived>::get();
  }

  // Vector-scalar division: v / s.
  Element operator/(Scalar rhs) const { return operator*(1.0 / rhs); }

  // Vector-scalar in-place division: v /= s.
  Element& operator/=(Scalar rhs) { return operator*=(1.0 / rhs); }

  // Scalar-vector multiplication: s * v.
  friend Element operator*(Scalar lhs, const Element& rhs) {
    return Element(lhs * rhs.coordinates_);
  }

  // Equality comparison.
  bool operator==(const Element& rhs) const {
    return coordinates_ == rhs.coordinates_;
  }

  // Inequality comparison.
  bool operator!=(const Element& rhs) const { return !(*this == rhs); }

 protected:
  // Construct from coordinate vector.
  explicit AlgebraElement(Vector coordinates)
      : coordinates_(std::move(coordinates)) {}

  // The underlying coordinate vector for this algebra element.
  Vector coordinates_;
};

}  // namespace mana
