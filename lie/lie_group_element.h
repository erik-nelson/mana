#include <Eigen/Dense>

#include "lie/algebra_element.h"
#include "lie/group_element.h"
#include "lie/lie_algebra_element.h"
#include "lie/manifold_element.h"

namespace mana {

// Base CRTP class for an element of a Lie group. Lie group elements are
// elements of both a group as well as a (differentiable, extrinsic) manifold.
template <typename Derived>
class LieGroupElement : public GroupElement<Derived>,
                        public ManifoldElement<Derived> {
 public:
  using Algebra = LieAlgebraElement<Derived>;

  Algebra log() const { return derived().logImpl(); }

  Derived exp(const Algebra& algebra) const {
    return derived().expImpl(algebra);
  }

  Algebra::Vector Log() const { return derived().LogImpl(); }

  Derived Exp(const Algebra::Vector& algebra) const {
    return derived().ExpImpl(algebra);
  }

  Jacobian Adjoint() const;
};

}  // namespace mana