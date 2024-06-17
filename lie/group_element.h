#pragma once

#include "utils/crtp.h"

namespace mana {

// Base CRTP class for an element of a group. Derived group element types should
// inherit from this like so:
//
//   class CyclicGroupElement : public GroupElement<CyclicGroupElement> {
//     ...
//   };
//
template <typename Derived>
class GroupElement : public Crtp<Derived> {
 public:
  // Identity element of this group.
  static Derived Identity() { return Derived::IdentityImpl(); }

  // Inverse of this group element.
  Derived Inverse() const { return Crtp<Derived>::get().InverseImpl(); }

  // Composition of group elements.
  Derived Compose(const Derived& rhs) const {
    return Crtp<Derived>::get().ComposeImpl(rhs);
  }

 private:
  // Derived classes must implement these methods.
  static Derived IdentityImpl();
  Derived InverseImpl() const;
  Derived ComposeImpl(const Derived& rhs) const;
};

}  // namespace mana