#pragma once

#include <iostream>
#include <type_traits>

#include "base/common.h"
#include "lie/constants.h"
#include "lie/lie_group_element.h"
#include "lie/traits.h"

namespace mana {

// Base class / abstract API for an element of a Lie algebra. Internally this is
// storing a coordinate vector for the corresponding Lie algebra element (not an
// element of the vector space). Templated on the derived Lie algebra type.
template <typename T>
class LieAlgebraElementBase {
 public:
  MANA_INHERIT_TRAITS(Traits<T>);

  // Return the identity element for this Lie algebra.
  static AlgebraElement Identity();

  // The hat operator converts from a Lie algebra vector to an element.
  AlgebraElementStorage Hat() const;

  // The vee operator converts from a Lie algebra element to a vector.
  static AlgebraElement Vee(const AlgebraElementStorage& element);

  // Return the corresponding Lie group element (upper-case Exp).
  GroupElement Exp() const;

  // The (upper-case) right- and left- Jacobians of the Exp map (see Eqs. 41 and
  // 44).
  Jacobian RightJacobian() const;
  Jacobian LeftJacobian() const;

  // Access underlying data.
  AlgebraVectorStorage& Storage();
  const AlgebraVectorStorage& Storage() const;

  // --------------------------------------------------------------------------
  // Induced API (built from methods above).

  // Return the Lie algebra coordinate vector for the provided Lie group element
  // (upper-case Log).
  static AlgebraElement Log(const GroupElement& element);

  // Access to coordinates of underlying data (checked, use operator[] for
  // unchecked access).
  Scalar& At(size_t idx);
  Scalar At(size_t idx) const;

  // Evaluate the Lie bracket of this and the provided element.
  AlgebraElement Compose(const AlgebraElement& rhs) const;

  // The (upper-case versions of) the right- and left- plus and minus operators
  // (see Eqs. 25-28). The `plus` and `minus` methods default to right-plus and
  // right-minus.
  GroupElement RightPlus(const GroupElement& rhs) const;
  GroupElement LeftPlus(const GroupElement& rhs) const;
  GroupElement Plus(const GroupElement& rhs) const;

  static AlgebraElement LeftMinus(const GroupElement& lhs,
                                  const GroupElement& rhs);
  static AlgebraElement RightMinus(const GroupElement& lhs,
                                   const GroupElement& rhs);
  static AlgebraElement Minus(const GroupElement& lhs, const GroupElement& rhs);

  // Return whether this Lie algebra element is approximately equal to `rhs`, up
  // to the provided scalar tolerance.
  bool EqualTo(const AlgebraElement& rhs,
               Scalar epsilon = Constants<Scalar>::kEpsilon) const;

  // --------------------------------------------------------------------------
  // Operator overloads.

  // Negation.
  AlgebraElement operator-() const;

  // Coordinate vector addition.
  AlgebraElement operator+(const AlgebraElement& rhs) const;

  // In-place coordinate vector addition.
  AlgebraElement& operator+=(const AlgebraElement& rhs);

  // Coordinate vector subtraction.
  AlgebraElement operator-(const AlgebraElement& rhs) const;

  // In-place coordinate vector subtraction.
  AlgebraElement& operator-=(const AlgebraElement& rhs);

  // Coordinate vector scaling.
  AlgebraElement operator*(Scalar rhs) const;
  AlgebraElement operator/(Scalar rhs) const;

  // In-place coordinate vector scaling.
  AlgebraElement& operator*=(Scalar rhs);
  AlgebraElement& operator/=(Scalar rhs);

  // Lie bracket and in-place Lie bracket.
  AlgebraElement operator*(const AlgebraElement& rhs) const;
  AlgebraElement& operator*=(const AlgebraElement& rhs);

  // Coordinate access.
  Scalar& operator[](size_t idx);
  Scalar operator[](size_t idx) const;

  // Approximate equality.
  bool operator==(const AlgebraElement& rhs) const;

  // Approximate inequality.
  bool operator!=(const AlgebraElement& rhs) const;

  // Ostreamable.
  inline friend std::ostream& operator<<(std::ostream& os,
                                         const AlgebraElement& rhs) {
    return (os << rhs.Storage());
  }

 protected:
  LieAlgebraElementBase() = default;

  // Cast to underlying derived type.
  T& derived() noexcept { return static_cast<T&>(*this); }
  const T& derived() const noexcept { return static_cast<const T&>(*this); }
};

// Impl -----------------------------------------------------------------------

// Multiplication with adjoint matrices: A * x.
template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement operator*(
    const typename LieAlgebraElementBase<T>::Jacobian& lhs,
    const LieAlgebraElementBase<T>& rhs) {
  return typename LieAlgebraElementBase<T>::AlgebraElement(
      (lhs * rhs.Storage()).eval());
}

template <typename T>
/*static*/ typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::Identity() {
  return T::Identity();
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElementStorage
LieAlgebraElementBase<T>::Hat() const {
  return derived().Hat();
}

template <typename T>
/*static*/ typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::Vee(const AlgebraElementStorage& element) {
  return T::Vee(element);
}

template <typename T>
typename LieAlgebraElementBase<T>::GroupElement LieAlgebraElementBase<T>::Exp()
    const {
  return derived().Exp();
}

template <typename T>
typename LieAlgebraElementBase<T>::Jacobian
LieAlgebraElementBase<T>::RightJacobian() const {
  return derived().RightJacobian();
}

template <typename T>
typename LieAlgebraElementBase<T>::Jacobian
LieAlgebraElementBase<T>::LeftJacobian() const {
  return derived().LeftJacobian();
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraVectorStorage&
LieAlgebraElementBase<T>::Storage() {
  return derived().Storage();
}

template <typename T>
const typename LieAlgebraElementBase<T>::AlgebraVectorStorage&
LieAlgebraElementBase<T>::Storage() const {
  return derived().Storage();
}

template <typename T>
/*static*/ typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::Log(const GroupElement& element) {
  return element.Log();
}

template <typename T>
typename LieAlgebraElementBase<T>::Scalar& LieAlgebraElementBase<T>::At(
    size_t idx) {
  return Storage().coeffRef(idx);
}

template <typename T>
typename LieAlgebraElementBase<T>::Scalar LieAlgebraElementBase<T>::At(
    size_t idx) const {
  return Storage().coeff(idx);
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::Compose(const AlgebraElement& rhs) const {
  const AlgebraElementStorage lhs_element = Hat();
  const AlgebraElementStorage rhs_element = rhs.Hat();
  return Vee(lhs_element * rhs_element - rhs_element * lhs_element);
}

template <typename T>
typename LieAlgebraElementBase<T>::GroupElement
LieAlgebraElementBase<T>::RightPlus(const GroupElement& rhs) const {
  return rhs.RightPlus(derived());
}

template <typename T>
typename LieAlgebraElementBase<T>::GroupElement
LieAlgebraElementBase<T>::LeftPlus(const GroupElement& rhs) const {
  return rhs.LeftPlus(derived());
}

template <typename T>
typename LieAlgebraElementBase<T>::GroupElement LieAlgebraElementBase<T>::Plus(
    const GroupElement& rhs) const {
  return rhs.Plus(derived());
}

template <typename T>
/*static*/ typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::LeftMinus(const GroupElement& lhs,
                                    const GroupElement& rhs) {
  return lhs.LeftMinus(rhs);
}

template <typename T>
/*static*/ typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::RightMinus(const GroupElement& lhs,
                                     const GroupElement& rhs) {
  return lhs.RightMinus(rhs);
}

template <typename T>
/*static*/ typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::Minus(const GroupElement& lhs,
                                const GroupElement& rhs) {
  return lhs.Minus(rhs);
}

template <typename T>
bool LieAlgebraElementBase<T>::EqualTo(const AlgebraElement& rhs,
                                       Scalar epsilon) const {
  return (Storage() - rhs.Storage()).array().abs().maxCoeff() <= epsilon;
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::operator-() const {
  return AlgebraElement(-Storage());
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::operator+(const AlgebraElement& rhs) const {
  return AlgebraElement(Storage() + rhs.Storage());
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement&
LieAlgebraElementBase<T>::operator+=(const AlgebraElement& rhs) {
  derived() = this->operator+(rhs);
  return derived();
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::operator-(const AlgebraElement& rhs) const {
  return AlgebraElement(Storage() - rhs.Storage());
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement&
LieAlgebraElementBase<T>::operator-=(const AlgebraElement& rhs) {
  derived() = this->operator-(rhs);
  return derived();
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::operator*(Scalar rhs) const {
  return AlgebraElement(Storage() * rhs);
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::operator/(Scalar rhs) const {
  return AlgebraElement(Storage() * (1.0 / rhs));
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement&
LieAlgebraElementBase<T>::operator*=(Scalar rhs) {
  derived() = this->operator*(rhs);
  return derived();
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement&
LieAlgebraElementBase<T>::operator/=(Scalar rhs) {
  derived() = this->operator/(rhs);
  return derived();
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement
LieAlgebraElementBase<T>::operator*(const AlgebraElement& rhs) const {
  return Compose(rhs);
}

template <typename T>
typename LieAlgebraElementBase<T>::AlgebraElement&
LieAlgebraElementBase<T>::operator*=(const AlgebraElement& rhs) {
  derived() = Compose(rhs);
  return derived();
}

template <typename T>
typename LieAlgebraElementBase<T>::Scalar& LieAlgebraElementBase<T>::operator[](
    size_t idx) {
  return Storage()[idx];
}

template <typename T>
typename LieAlgebraElementBase<T>::Scalar LieAlgebraElementBase<T>::operator[](
    size_t idx) const {
  return Storage()[idx];
}

template <typename T>
bool LieAlgebraElementBase<T>::operator==(const AlgebraElement& rhs) const {
  return EqualTo(rhs);
}

template <typename T>
bool LieAlgebraElementBase<T>::operator!=(const AlgebraElement& rhs) const {
  return !EqualTo(rhs);
}

}  // namespace mana