#pragma once

namespace mana {

// Helper for implemented Curiously Recurring Template Pattern (CRTP). Defines
// methods that cast the CRTP base type to the derived type (via
// `Crtp<Derived>::get()` methods).
//
// Example usage:
//
//     // The CRTP base class inherits from Crtp<Derived>, giving it access to
//     // the `get()` method used to cast itself to the `Derived` type.
//     template <typename Derived>
//     class Base : public Crtp<Derived> {
//      public:
//       void Foo() {
//         // Cast `Base` to `Derived`.
//         Derived& derived = Crtp<Derived>::get();
//
//         // Call `FooImpl` on underlying derived type.
//         derived.FooImpl();
//       }
//     };
//
//     class Derived : public Base<Derived> {
//      public:
//       void FooImpl() {
//         ...
//       }
//     };
//
template <typename Derived>
class Crtp {
 protected:
  Derived& get() { return static_cast<Derived&>(*this); }
  const Derived& get() const { return static_cast<const Derived&>(*this); }
};

}  // namespace mana