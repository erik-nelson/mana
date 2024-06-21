#include <Eigen/Dense>
#include <cassert>
#include <type_traits>

#include "lie/base/algebra_element.h"
#include "lie/base/group_element.h"
#include "lie/base/lie_algebra_element.h"
#include "lie/base/manifold_element.h"

namespace mana {

// Base CRTP class for an element of a Lie group. Lie group elements are
// elements of both a group as well as a (differentiable, extrinsic) manifold.
//
// Derived classes must implement all methods required of a `GroupElement` and a
// `ManifoldElement` (see group_element.h and manifold_element.h), as well as
// these methods:
// - TangentVector LogImpl() const;
// - static GroupElement ExpImpl(const TangentVector& coordinate);
// - Jacobian AdjointImpl() const;
template <typename Derived>
class LieGroupElement : public GroupElement<Derived>,
                        public ManifoldElement<Derived> {
 public:
  // Traits inherited from the fact that we are a manifold.
  using Scalar = typename ManifoldTraits<Derived>::Scalar;
  using Chart = typename ManifoldTraits<Derived>::Chart;
  using Geodesic = typename ManifoldTraits<Derived>::Geodesic;
  using TangentVector = typename ManifoldTraits<Derived>::TangentVector;
  using EmbeddingPoint = typename ManifoldTraits<Derived>::EmbeddingPoint;
  static constexpr int Dimension = ManifoldTraits<Derived>::Dimension;
  static constexpr int EmbeddingDimension =
      ManifoldTraits<Derived>::EmbeddingDimension;

  // Traits specific to Lie groups.
  using GroupElement = Derived;
  using AlgebraElement = typename LieAlgebraElement<Derived>::AlgebraElement;
  using Jacobian = Eigen::Matrix<Scalar, Dimension, Dimension>;

  // Trait checks: Ensure this Lie group is compatible with the associated
  // algebra type.
  static_assert(std::is_same_v<AlgebraElement::Scalar, Scalar>);
  static_assert(std::is_same_v<AlgebraElement::Vector, TangentVector>);
  static_assert(std::is_same_v<AlgebraElement::Dimension, Dimension>);

  // Lower-case `log` map: return the Lie algebra element for this Lie group
  // element.
  AlgebraElement log() const;

  // Upper-case `Log` map: return the Lie algebra coordinate vector for this Lie
  // group element.
  TangentVector Log() const;

  // Lower-case `exp` map: construct a Lie group element from corresponding Lie
  // algebra element.
  static GroupElement exp(const AlgebraElement& algebra);

  // Upper-case `Exp` map: construct a Lie group element from corresponding Lie
  // algebra coordinate vector.
  static GroupElement Exp(const TangentVector& coordinate);

  // Return the adjoint of this Lie group element. The adjoint is a square
  // matrix that maps coordinate vectors from the tangent space of this group
  // element to the Lie algebra (the tangent space at the identity element).
  Jacobian Adjoint() const;

  // The (lower-case versions of) right- and left- plus and minus operators (see
  // Eqs. 25-28 in A micro Lie theory for state estimation in robotics). These
  // functions operate on / return Lie algebra elements (not tangent vectors).
  // The `plus` and `minus` methods default to right-plus and right-minus.
  GroupElement rplus(const AlgebraElement& rhs) const;
  GroupElement lplus(const AlgebraElement& lhs) const;
  GroupElement plus(const AlgebraElement& rhs) const;
  AlgebraElement rminus(const GroupElement& rhs) const;
  AlgebraElement lminus(const GroupElement& rhs) const;
  AlgebraElement minus(const GroupElement& rhs) const;

  // The (upper-case versions of) right- and left- plus and minus operators (see
  // Eqs. 25-28 in A micro Lie theory for state estimation in robotics). These
  // functions operate on / return tangent vectors in the Lie algebra. The
  // `Plus` and `Minus` methods default to right-plus and right-minus.
  GroupElement Rplus(const TangentVector& rhs) const;
  GroupElement Lplus(const TangentVector& lhs) const;
  GroupElement Plus(const TangentVector& rhs) const;
  TangentVector Rminus(const GroupElement& rhs) const;
  TangentVector Lminus(const GroupElement& rhs) const;
  TangentVector Minus(const GroupElement& rhs) const;

  // Override base class's DistanceToImpl and InterpolateImpl methods. Lie
  // groups have a well defined standard distance and interpolation available.
  Scalar DistanceToImpl(const GroupElement& rhs) const;
  GroupElement InterpolateImpl(const GroupElement& rhs, Scalar fraction) const;

 private:
  // CRTP helpers.
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }
};

template <typename Derived>
typename LieGroupElement<Derived>::AlgebraElement
LieGroupElement<Derived>::log() const {
  return AlgebraElement::Hat(Log());
}

template <typename Derived>
typename LieGroupElement<Derived>::TangentVector LieGroupElement<Derived>::Log()
    const {
  return derived().LogImpl();
}

template <typename Derived>
/*static*/ typename LieGroupElement<Derived>::GroupElement
LieGroupElement<Derived>::exp(const AlgebraElement& algebra) {
  return Exp(algebra.Vee());
}

template <typename Derived>
/*static*/ typename LieGroupElement<Derived>::GroupElement
LieGroupElement<Derived>::Exp(const TangentVector& coordinate) {
  return derived().ExpImpl(coordinate);
}

template <typename Derived>
typename LieGroupElement<Derived>::Jacobian LieGroupElement<Derived>::Adjoint()
    const {
  return derived().AdjointImpl();
}

template <typename Derived>
typename LieGroupElement<Derived>::GroupElement LieGroupElement<Derived>::rplus(
    const AlgebraElement& rhs) const {
  return Rplus(rhs.Vee());
}

template <typename Derived>
typename LieGroupElement<Derived>::GroupElement LieGroupElement<Derived>::lplus(
    const AlgebraElement& lhs) const {
  return Lplus(lhs.Vee());
}

template <typename Derived>
typename LieGroupElement<Derived>::GroupElement LieGroupElement<Derived>::plus(
    const AlgebraElement& rhs) const {
  return Plus(rhs.Vee());
}

template <typename Derived>
typename LieGroupElement<Derived>::AlgebraElement
LieGroupElement<Derived>::rminus(const GroupElement& rhs) const {
  return BetweenInner(rhs).log();
}

template <typename Derived>
typename LieGroupElement<Derived>::AlgebraElement
LieGroupElement<Derived>::lminus(const GroupElement& rhs) const {
  return rhs.BetweenOuter(*this).log();
}

template <typename Derived>
typename LieGroupElement<Derived>::AlgebraElement
LieGroupElement<Derived>::minus(const GroupElement& rhs) const {
  return rminus(rhs);
}

template <typename Derived>
typename LieGroupElement<Derived>::GroupElement LieGroupElement<Derived>::Rplus(
    const TangentVector& rhs) const {
  return Compose(Exp(rhs));
}

template <typename Derived>
typename LieGroupElement<Derived>::GroupElement LieGroupElement<Derived>::Lplus(
    const TangentVector& lhs) const {
  return Exp(lhs).Compose(*this);
}

template <typename Derived>
typename LieGroupElement<Derived>::GroupElement LieGroupElement<Derived>::Plus(
    const TangentVector& rhs) const {
  return Rplus(rhs);
}

template <typename Derived>
typename LieGroupElement<Derived>::TangentVector
LieGroupElement<Derived>::Rminus(const GroupElement& rhs) const {
  return rminus(rhs).Vee();
}

template <typename Derived>
typename LieGroupElement<Derived>::TangentVector
LieGroupElement<Derived>::Lminus(const GroupElement& rhs) const {
  return lminus(rhs).Vee();
}

template <typename Derived>
typename LieGroupElement<Derived>::TangentVector
LieGroupElement<Derived>::Minus(const GroupElement& rhs) const {
  return minus(rhs).Vee();
}

template <typename Derived>
typename LieGroupElement<Derived>::Scalar
LieGroupElement<Derived>::DistanceToImpl(const GroupElement& rhs) const {
  return Minus().norm();
}

template <typename Derived>
typename LieGroupElement<Derived>::GroupElement
LieGroupElement<Derived>::InterpolateImpl(const GroupElement& rhs,
                                          Scalar fraction) const {
  return Plus(fraction * rhs.Log());
}

}  // namespace mana