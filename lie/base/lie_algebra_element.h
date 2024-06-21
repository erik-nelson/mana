
#include <Eigen/Dense>

#include "lie/base/algebra_element.h"

namespace mana {

// Forward declaration.
template <typename Derived>
class LieGroupElement;

// Base CRTP class for an element of a Lie algebra. Lie algebra elements are
// vectors in the tangent space at the identity of a corresponding Lie group,
// that also have a product operation defined on them (the Lie bracket).
template <typename Derived>
class LieAlgebraElement : public AlgebraElement<Derived> {
 public:
  // Traits inherited from the fact that we are an algebra.
  using Scalar = typename AlgebraTraits<Derived>::Scalar;
  using Vector = typename AlgebraTraits<Derived>::Vector;
  static constexpr int Dimension = AlgebraTraits<Derived>::Dimension;

  // Traits specific to Lie algebras.
  using AlgebraElement = Derived;
  using GroupElement = LieGroupElement<Derived>;
};

}  // namespace mana