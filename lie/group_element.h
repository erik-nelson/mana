#pragma once

namespace mana {

// Base CRTP class for an element of a group.
//
// Derived classes must implement these methods.
// - static Derived IdentityImpl();
// - Derived InverseImpl() const;
// - Derived ComposeImpl(const Derived& rhs) const;
template <typename Derived>
class GroupElement {
 public:
  // Identity element of this group.
  static Derived Identity() { return Derived::IdentityImpl(); }

  // Inverse of this group element.
  Derived Inverse() const { return derived().InverseImpl(); }

  // Composition of group elements.
  Derived Compose(const Derived& rhs) const {
    return derived().ComposeImpl(rhs);
  }

  // Helper for X^{-1} * Y.
  Derived BetweenInner(const Derived& rhs) const {
    return Inverse().Compose(rhs);
  }

  // Helper for X * Y^{-1}.
  Derived BetweenOuter(const Derived& rhs) const {
    return Compose(rhs.Inverse());
  }

 private:
  // CRTP helpers.
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }
};

}  // namespace mana