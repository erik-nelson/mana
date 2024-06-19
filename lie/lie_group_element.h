#include <Eigen/Dense>
#include <cassert>
#include <type_traits>

#include "lie/algebra_element.h"
#include "lie/group_element.h"
#include "lie/lie_algebra_element.h"
#include "lie/manifold_element.h"

namespace mana {

// Base CRTP class for an element of a Lie group. Lie group elements are
// elements of both a group as well as a (differentiable, extrinsic) manifold.
//
// Derived classes must implement all methods required of a `GroupElement` and a
// `ManifoldElement` (see group_element.h and manifold_element.h), as well as
// these methods:
// - AlgebraElement logImpl() const;
// - TangentVector LogImpl() const;
// - static GroupElement expImpl(const AlgebraElement& algebra);
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
  AlgebraElement log() const { return derived().logImpl(); }

  // Upper-case `Log` map: return the Lie algebra coordinate vector for this Lie
  // group element.
  TangentVector Log() const { return derived().LogImpl(); }

  // Lower-case `exp` map: construct a Lie group element from corresponding Lie
  // algebra element.
  static GroupElement exp(const AlgebraElement& algebra) {
    return derived().expImpl(algebra);
  }

  // Upper-case `Exp` map: construct a Lie group element from corresponding Lie
  // algebra coordinate vector.
  static GroupElement Exp(const TangentVector& coordinate) {
    return derived().ExpImpl(coordinate);
  }

  // Return the adjoint of this Lie group element. The adjoint is a square
  // matrix that maps coordinate vectors from the tangent space of this group
  // element to the Lie algebra (the tangent space at the identity element).
  Jacobian Adjoint() const { return derived().AdjointImpl(); }

  // The (lower-case versions of) right- and left- plus and minus operators (see
  // Eqs. 25-28 in A micro Lie theory for state estimation in robotics). These
  // functions operate on / return Lie algebra elements (not tangent vectors).
  // The `plus` and `minus` methods default to right-plus and right-minus.
  GroupElement rplus(const AlgebraElement& rhs) const {
    return Compose(exp(rhs));
  }
  GroupElement lplus(const AlgebraElement& lhs) const {
    return exp(lhs).Compose(*this);
  }
  GroupElement plus(const AlgebraElement& rhs) const { return rplus(rhs); }

  AlgebraElement rminus(const GroupElement& rhs) const {
    return BetweenInner(rhs).log();
  }
  AlgebraElement Lminus(const GroupElement& rhs) const {
    return rhs.BetweenOuter(*this).log();
  }
  AlgebraElement minus(const GroupElement& rhs) const { return rminus(rhs); }

  // The (upper-case versions of) right- and left- plus and minus operators (see
  // Eqs. 25-28 in A micro Lie theory for state estimation in robotics). These
  // functions operate on / return tangent vectors in the Lie algebra. The
  // `Plus` and `Minus` methods default to right-plus and right-minus.
  GroupElement Rplus(const TangentVector& rhs) const {
    return Compose(Exp(rhs));
  }
  GroupElement Lplus(const TangentVector& lhs) const {
    return Exp(lhs).Compose(*this);
  }
  GroupElement Plus(const TangentVector& rhs) const { return Rplus(rhs); }

  TangentVector Rminus(const GroupElement& rhs) const {
    return BetweenInner(rhs).Log();
  }
  TangentVector Lminus(const GroupElement& rhs) const {
    return rhs.BetweenOuter(*this).Log();
  }
  TangentVector Minus(const GroupElement& rhs) const { return Rminus(rhs); }

 private:
  // CRTP helpers.
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }
};

}  // namespace mana