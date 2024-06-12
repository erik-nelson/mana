#pragma once

#include <iostream>
#include <type_traits>

#include "base/common.h"
#include "lie/constants.h"
#include "lie/lie_algebra_element.h"
#include "lie/traits.h"

namespace mana {

// Base class / abstract API for an element of a Lie group. Templated on the
// derived Lie group type.
template <typename T>
class LieGroupElementBase {
 public:
  MANA_INHERIT_TRAITS(Traits<T>);

  // Return the identity element for this Lie group.
  static GroupElement Identity();

  // Return the inverse of this Lie group element.
  GroupElement Inverse() const;

  // Invert this Lie group element in place.
  void Invert();

  // Return the corresponding Lie algebra coordinate vector (upper-case log).
  AlgebraElement Log() const;

  // Return the corresponding Lie algebra element (lower-case log).
  AlgebraElementStorage log() const;

  // Return the Lie group element for the provided Lie algebra coordinate
  // vector (upper-case exp).
  static GroupElement Exp(const AlgebraElement& element);

  // Return the Lie group element for the provided Lie algebra element
  // (lower-case exp).
  static GroupElement exp(const AlgebraElementStorage& element);

  // Compose this Lie group element with the provided `rhs` Lie group element.
  GroupElement Compose(const GroupElement& rhs) const;

  // Return the adjoint of this Lie group element.
  Jacobian Adjoint() const;

  // TODO(erik): Act? How do we define a manifold? Should this class also
  // inherit from Manifold base?

  // Access underlying data.
  GroupElementStorage& Storage();
  const GroupElementStorage& Storage() const;

  // --------------------------------------------------------------------------
  // Induced API (built from methods above).

  // The (upper-case versions of) right- and left- plus and minus operators (see
  // Eqs. 25-28). These functions operate on / return Lie algebra coordinate
  // vectors. The `Plus` and `Minus` methods default to right-plus and
  // right-minus.
  GroupElement RightPlus(const AlgebraElement& rhs) const;
  GroupElement LeftPlus(const AlgebraElement& rhs) const;
  GroupElement Plus(const AlgebraElement& rhs) const;
  AlgebraElement RightMinus(const GroupElement& rhs) const;
  AlgebraElement LeftMinus(const GroupElement& rhs) const;
  AlgebraElement Minus(const GroupElement& rhs) const;

  // Inner between: computes inv(lhs) * rhs.
  GroupElement BetweenInner(const GroupElement& rhs) const;

  // Outer between: computes lhs * inv(rhs).
  GroupElement BetweenOuter(const GroupElement& rhs) const;

  // Return whether this Lie group element is approximately equal to `rhs`, up
  // to the provided scalar tolerance. This check is performed in the Lie
  // algebra.
  bool EqualTo(const GroupElement& rhs,
               Scalar epsilon = Constants<Scalar>::kEps) const;

  // --------------------------------------------------------------------------
  // Operator overloads.

  // Addition: uses right-plus operations.
  GroupElement operator+(const AlgebraElement& rhs) const;

  // In-place addition: uses right-plus operations.
  GroupElement& operator+=(const AlgebraElement& rhs) const;

  // Subtraction: uses left-minus operation.
  AlgebraElement operator-(const GroupElement& rhs) const;

  // Composition.
  GroupElement operator*(const GroupElement& rhs) const;

  // In-place composition.
  GroupElement& operator*=(const GroupElement& rhs) const;

  // Approximate equality.
  bool operator==(const GroupElement& rhs) const;

  // Approximate inequality.
  bool operator!=(const GroupElement& rhs) const;

  // Multiplication with Lie algebra elements (not coordinate vectors): G * A.
  AlgebraElementStorage operator*(const AlgebraElementStorage& rhs) const;

  // Multiplication with Lie algebra elements (not coordinate vectors): A * G.
  inline friend AlgebraElementStorage operator*(
      const AlgebraElementStorage& lhs, const GroupElement& rhs) {
    return GroupElement(lhs).Compose(rhs).Storage();
  }

  // Ostreamable.
  inline friend std::ostream& operator<<(std::ostream& os,
                                         const GroupElement& rhs) {
    return (os << rhs.Storage());
  }

 protected:
  LieGroupElementBase() = default;

  // Cast to underlying derived type.
  T& derived() noexcept { return static_cast<T&>(*this); }
  const T& derived() const noexcept { return static_cast<const T&>(*this); }
};

// Impl -----------------------------------------------------------------------
template <typename T>
/*static*/ typename LieGroupElementBase<T>::GroupElement
LieGroupElementBase<T>::Identity() {
  return T::Identity();
}

template <typename T>
typename LieGroupElementBase<T>::GroupElement LieGroupElementBase<T>::Inverse()
    const {
  return derived().Inverse();
}

template <typename T>
void LieGroupElementBase<T>::Invert() {
  derived().Invert();
}

template <typename T>
typename LieGroupElementBase<T>::AlgebraElement LieGroupElementBase<T>::Log()
    const {
  return derived().Log();
}

template <typename T>
typename LieGroupElementBase<T>::AlgebraElementStorage
LieGroupElementBase<T>::log() const {
  return Log().Hat();
}

template <typename T>
/*static*/ typename LieGroupElementBase<T>::GroupElement
LieGroupElementBase<T>::Exp(const AlgebraElement& element) {
  return element.Exp();
}

template <typename T>
/*static*/ typename LieGroupElementBase<T>::GroupElement
LieGroupElementBase<T>::exp(const AlgebraElementStorage& element) {
  return AlgebraElement::Vee(element).Exp();
}

template <typename T>
typename LieGroupElementBase<T>::GroupElement LieGroupElementBase<T>::Compose(
    const GroupElement& rhs) const {
  return derived().Compose(rhs);
}

template <typename T>
typename LieGroupElementBase<T>::Jacobian LieGroupElementBase<T>::Adjoint()
    const {
  return derived().Adjoint();
}

template <typename T>
typename LieGroupElementBase<T>::GroupElementStorage&
LieGroupElementBase<T>::Storage() {
  return derived().Storage();
}

template <typename T>
const typename LieGroupElementBase<T>::GroupElementStorage&
LieGroupElementBase<T>::Storage() const {
  return derived().Storage();
}

template <typename T>
typename LieGroupElementBase<T>::GroupElement LieGroupElementBase<T>::RightPlus(
    const AlgebraElement& rhs) const {
  return Compose(rhs.Exp());
}

template <typename T>
typename LieGroupElementBase<T>::GroupElement LieGroupElementBase<T>::LeftPlus(
    const AlgebraElement& rhs) const {
  return rhs.Exp().Compose(derived());
}

template <typename T>
typename LieGroupElementBase<T>::GroupElement LieGroupElementBase<T>::Plus(
    const AlgebraElement& rhs) const {
  return RightPlus(rhs);
}

template <typename T>
typename LieGroupElementBase<T>::AlgebraElement
LieGroupElementBase<T>::RightMinus(const GroupElement& rhs) const {
  return Inverse().Compose(rhs).Log();
}

template <typename T>
typename LieGroupElementBase<T>::AlgebraElement
LieGroupElementBase<T>::LeftMinus(const GroupElement& rhs) const {
  return Compose(rhs.Inverse()).Log();
}

template <typename T>
typename LieGroupElementBase<T>::AlgebraElement LieGroupElementBase<T>::Minus(
    const GroupElement& rhs) const {
  return RightMinus(rhs);
}

template <typename T>
typename LieGroupElementBase<T>::GroupElement
LieGroupElementBase<T>::BetweenInner(const GroupElement& rhs) const {
  return Inverse().Compose(rhs);
}

template <typename T>
typename LieGroupElementBase<T>::GroupElement
LieGroupElementBase<T>::BetweenOuter(const GroupElement& rhs) const {
  return Compose(rhs.Inverse());
}

template <typename T>
bool LieGroupElementBase<T>::EqualTo(const GroupElement& rhs,
                                     Scalar epsilon) const {
  return Minus(rhs).EqualTo(AlgebraElement::Identity(), epsilon);
}

template <typename T>
typename LieGroupElementBase<T>::GroupElement LieGroupElementBase<T>::operator+(
    const AlgebraElement& rhs) const {
  return Plus(rhs);
}

template <typename T>
typename LieGroupElementBase<T>::GroupElement&
LieGroupElementBase<T>::operator+=(const AlgebraElement& rhs) const {
  derived() = Plus(rhs);
  return derived();
}

template <typename T>
typename LieGroupElementBase<T>::AlgebraElement
LieGroupElementBase<T>::operator-(const GroupElement& rhs) const {
  return Minus(rhs);
}

template <typename T>
typename LieGroupElementBase<T>::GroupElement LieGroupElementBase<T>::operator*(
    const GroupElement& rhs) const {
  return Compose(rhs);
}

template <typename T>
typename LieGroupElementBase<T>::GroupElement&
LieGroupElementBase<T>::operator*=(const GroupElement& rhs) const {
  derived() = Compose(rhs);
  return derived();
}

template <typename T>
bool LieGroupElementBase<T>::operator==(const GroupElement& rhs) const {
  return EqualTo(rhs);
}

template <typename T>
bool LieGroupElementBase<T>::operator!=(const GroupElement& rhs) const {
  return !EqualTo(rhs);
}

template <typename T>
typename LieGroupElementBase<T>::AlgebraElementStorage
LieGroupElementBase<T>::operator*(
    const typename LieGroupElementBase<T>::AlgebraElementStorage& rhs) const {
  return Compose(GroupElement(rhs)).Storage();
}

}  // namespace mana