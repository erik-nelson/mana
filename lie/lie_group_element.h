#include <Eigen/Dense>

#include "lie/algebra_element.h"
#include "lie/group_element.h"
#include "lie/manifold_element.h"

namespace mana {

template <typename Derived>
class LieGroupElement : public GroupElement<Derived>,
                        public ManifoldElement<Derived> {
 public:
#if 0
  LieAlgebraElement Log() const { return derived().LogImpl(); }

  Derived Exp(const LieAlgebraElement& algebra_element) const {
    return derived().ExpImpl(algebra_element);
  }

  Derived Geodesic(const Derived& other, double t) const {
    return derived().geodesic_impl(other, t);
  }

  double Distance(const Derived& other) const {
    return derived().DistanceImpl(other);
  }

  double InnerProduct(const Eigen::VectorXd& vec1,
                      const Eigen::VectorXd& vec2) const {
    return derived().InnerProductImpl(vec1, vec2);
  }

  std::vector<Eigen::VectorXd> TangentSpaceBasis() const {
    return derived().TangentSpaceBasisImpl();
  }
#endif
};

}  // namespace mana