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
  static Derived Identity();

  // Inverse of this group element.
  Derived Inverse() const;

  // Composition of group elements.
  Derived Compose(const Derived& rhs) const;

  // Helper for X^{-1} * Y.
  Derived BetweenInner(const Derived& rhs) const;

  // Helper for X * Y^{-1}.
  Derived BetweenOuter(const Derived& rhs) const;

 private:
  // CRTP helpers.
  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }
};

template <typename Derived>
/*static*/ Derived GroupElement<Derived>::Identity() {
  return Derived::IdentityImpl();
}

template <typename Derived>
Derived GroupElement<Derived>::Inverse() const {
  return derived().InverseImpl();
}

template <typename Derived>
Derived GroupElement<Derived>::Compose(const Derived& rhs) const {
  return derived().ComposeImpl(rhs);
}

template <typename Derived>
Derived GroupElement<Derived>::BetweenInner(const Derived& rhs) const {
  return Inverse().Compose(rhs);
}

template <typename Derived>
Derived GroupElement<Derived>::BetweenOuter(const Derived& rhs) const {
  return Compose(rhs.Inverse());
}

}  // namespace mana