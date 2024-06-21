#pragma once

#include <Eigen/Dense>

#include "lie/base/algebra_element.h"

namespace mana {

// Traits template for a Lie algebra.
template <typename Derived>
struct LieAlgebraTraits {};

// Base CRTP class for an element of a Lie algebra. Lie algebra elements are
// vectors in the tangent space at the identity of a corresponding Lie group,
// that also have a product operation defined on them (the Lie bracket).
//
// Derived classes must implement all methods required of a `AlgebraElement`, as
// well as these methods:
// - explicit LieAlgebraElement(const Matrix& matrix); // (constructor)
// - Matrix AsMatrixImpl() const;
template <typename Derived>
class LieAlgebraElement : public AlgebraElement<Derived> {
 public:
  // Traits inherited from the fact that we are an algebra.
  using Scalar = typename AlgebraTraits<Derived>::Scalar;
  using Vector = typename AlgebraTraits<Derived>::Vector;
  static constexpr int Dimension = AlgebraTraits<Derived>::Dimension;

  // Traits inherited from the fact that we are a Lie algebra.
  using AlgebraElement = typename LieAlgebraTraits<Derived>::AlgebraElement;
  using GroupElement = typename LieAlgebraTraits<Derived>::GroupElement;
  using Matrix = typename LieAlgebraTraits<Derived>::Matrix;

  // Return the matrix representation of this Lie algebra element.
  Matrix AsMatrix() const;

  // The Lie bracket of this element and `rhs`.
  AlgebraElement Bracket(const AlgebraElement& rhs) const;

  // Lower-case `log` map: construct a Lie algebra element from corresponding
  // Lie group element.
  static AlgebraElement log(const GroupElement& group);

  // Lower-case `exp` map: return the Lie group element for this Lie algebra
  // element.
  GroupElement exp() const;

  // Override base class's ComposeImpl method. (Matrix) Lie algebras have a well
  // defined composition available.
  Matrix ComposeImpl(const AlgebraElement& rhs) const;

 private:
  // CRTP helpers.
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }
};

template <typename Derived>
typename LieAlgebraElement<Derived>::Matrix
LieAlgebraElement<Derived>::AsMatrix() const {
  return derived().AsMatrixImpl();
}

template <typename Derived>
typename LieAlgebraElement<Derived>::AlgebraElement
LieAlgebraElement<Derived>::Bracket(const AlgebraElement& rhs) const {
  return Compose(rhs) - rhs.Compose(derived());
}

template <typename Derived>
/*static*/ typename LieAlgebraElement<Derived>::AlgebraElement
LieAlgebraElement<Derived>::log(const GroupElement& group) {
  return group.log();
}

template <typename Derived>
typename LieAlgebraElement<Derived>::GroupElement
LieAlgebraElement<Derived>::exp() const {
  return GroupElement::exp(derived());
}

template <typename Derived>
typename LieAlgebraElement<Derived>::Matrix
LieAlgebraElement<Derived>::ComposeImpl(const AlgebraElement& rhs) const {
  return AlgebraElement(AsMatrix() * rhs.AsMatrix());
}

}  // namespace mana