#pragma once

#include <Eigen/Dense>
#include <utility>

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
// Derived classes must implement these methods:
// - explicit Element(Vector coordinates); // (constructor)
// - Element ComposeImpl(const Element& rhs) const;
template <typename Derived>
class AlgebraElement {
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
  Vector& Coordinates();
  const Vector& Coordinates() const;

  // The hat operator builds an algebra element from a coordinate vector.
  static Element Hat(Vector coordinate);

  // The vee operator builds a coordinate vector from this algebra element. This
  // is synonmous with `this->Coordinates()`.
  const Vector& Vee() const;

  // The bilinear product of this algebra.
  Element Compose(const Element& rhs) const;

  // Vector addition: a + b.
  Element operator+(const Element& rhs) const;

  // Vector in-place addition: a += b.
  Element& operator+=(const Element& rhs);

  // Vector subtraction: a - b.
  Element operator-(const Element& rhs) const;

  // Vector in-place subtraction: a -= b.
  Element& operator-=(const Element& rhs);

  // Vector product: a * b.
  Element operator*(const Element& rhs) const;

  // Vector in-place product: a *= b.
  Element& operator*=(const Element& rhs);

  // Vector negation: -a.
  Element operator-() const;

  // Vector-scalar multiplication: v * s.
  Element operator*(Scalar rhs) const;

  // Vector-scalar in-place multiplication: v *= s.
  Element& operator*=(Scalar rhs);

  // Vector-scalar division: v / s.
  Element operator/(Scalar rhs) const;

  // Vector-scalar in-place division: v /= s.
  Element& operator/=(Scalar rhs);

  // Equality comparison.
  bool operator==(const Element& rhs) const;

  // Inequality comparison.
  bool operator!=(const Element& rhs) const;

  // Scalar-vector multiplication: s * v.
  friend Element operator*(Scalar lhs, const Element& rhs) { return rhs * lhs; }

 protected:
  // Construct from coordinate vector.
  explicit AlgebraElement(Vector coordinates);

  // The underlying coordinate vector for this algebra element.
  Vector coordinates_;

 private:
  // CRTP helpers.
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }
};

template <typename Derived>
typename AlgebraElement<Derived>::Vector&
AlgebraElement<Derived>::Coordinates() {
  return coordinates_;
}

template <typename Derived>
const typename AlgebraElement<Derived>::Vector&
AlgebraElement<Derived>::Coordinates() const {
  return coordinates_;
}

template <typename Derived>
/*static*/ typename AlgebraElement<Derived>::Element
AlgebraElement<Derived>::Hat(Vector coordinate) {
  return Element(std::move(coordinate));
}

template <typename Derived>
const typename AlgebraElement<Derived>::Vector& AlgebraElement<Derived>::Vee()
    const {
  return Coordinates();
}

template <typename Derived>
typename AlgebraElement<Derived>::Element AlgebraElement<Derived>::Compose(
    const Element& rhs) const {
  return derived().ComposeImpl(rhs);
}

template <typename Derived>
typename AlgebraElement<Derived>::Element AlgebraElement<Derived>::operator+(
    const Element& rhs) const {
  return Element(coordinates_ + rhs.coordinates_);
}

template <typename Derived>
typename AlgebraElement<Derived>::Element& AlgebraElement<Derived>::operator+=(
    const Element& rhs) {
  coordinates_ += rhs.coordinates_;
  return derived();
}

template <typename Derived>
typename AlgebraElement<Derived>::Element AlgebraElement<Derived>::operator-(
    const Element& rhs) const {
  return Element(coordinates_ - rhs.coordinates_);
}

template <typename Derived>
typename AlgebraElement<Derived>::Element& AlgebraElement<Derived>::operator-=(
    const Element& rhs) {
  coordinates_ -= rhs.coordinates_;
  return derived();
}

template <typename Derived>
typename AlgebraElement<Derived>::Element AlgebraElement<Derived>::operator*(
    const Element& rhs) const {
  return Compose(rhs);
}

template <typename Derived>
typename AlgebraElement<Derived>::Element& AlgebraElement<Derived>::operator*=(
    const Element& rhs) {
  Element result = Compose(rhs);
  coordinates_ = std::move(result.coordinates_);
  return derived();
}

template <typename Derived>
typename AlgebraElement<Derived>::Element AlgebraElement<Derived>::operator-()
    const {
  return Element(-coordinates_);
}

template <typename Derived>
typename AlgebraElement<Derived>::Element AlgebraElement<Derived>::operator*(
    Scalar rhs) const {
  return Element(coordinates_ * rhs);
}

template <typename Derived>
typename AlgebraElement<Derived>::Element& AlgebraElement<Derived>::operator*=(
    Scalar rhs) {
  coordinates_ *= rhs;
  return derived();
}

template <typename Derived>
typename AlgebraElement<Derived>::Element AlgebraElement<Derived>::operator/(
    Scalar rhs) const {
  return operator*(1.0 / rhs);
}

template <typename Derived>
typename AlgebraElement<Derived>::Element& AlgebraElement<Derived>::operator/=(
    Scalar rhs) {
  return operator*=(1.0 / rhs);
}

template <typename Derived>
bool AlgebraElement<Derived>::operator==(const Element& rhs) const {
  return coordinates_ == rhs.coordinates_;
}

template <typename Derived>
bool AlgebraElement<Derived>::operator!=(const Element& rhs) const {
  return !(*this == rhs);
}

template <typename Derived>
AlgebraElement<Derived>::AlgebraElement(Vector coordinates)
    : coordinates_(std::move(coordinates)) {}

}  // namespace mana
