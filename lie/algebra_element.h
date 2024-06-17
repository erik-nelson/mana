#pragma once

#include <Eigen/Dense>

#include "utils/crtp.h"

namespace mana {

template <typename Derived, typename CoordinateType>
class AlgebraElementBase : public Crtp<Derived> {
 public:
  Eigen::MatrixXd& Data() { return derived().DataImpl(); }
  const Eigen::MatrixXd& Data() const { return derived().DataImpl(); }

  Derived Hat(const CoordinateType& coords) const {
    return derived().HatImpl(coords);
  }

  CoordinateType Vee(const Derived& algebra_element) const {
    return derived().VeeImpl(algebra_element);
  }

 private:
  // Derived classes must implement these methods.
  Eigen::MatrixXd DataImpl() const;
  Derived HatImpl(const CoordinateType& coords) const;
  CoordinateType VeeImpl(const Derived& algebra_element) const;
};

}  // namespace mana
