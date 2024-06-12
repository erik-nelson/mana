#pragma once

namespace mana {

template <typename DerivedT>
struct Traits;

// Helper macro for inheriting traits for a Lie group element, algebra element,
// or algebra vector.
#define MANA_INHERIT_TRAITS(X)                                   \
  static constexpr size_t Dims = X::Dims;                        \
  static constexpr size_t Dofs = X::Dofs;                        \
  using Scalar = typename X::Scalar;                             \
  using Jacobian = typename X::Jacobian;                         \
  using GroupElement = typename X::GroupElement;                 \
  using AlgebraElement = typename X::AlgebraElement;             \
  using GroupElementStorage = typename X::GroupElementStorage;   \
  using AlgebraVectorStorage = typename X::AlgebraVectorStorage; \
  using AlgebraElementStorage = typename X::AlgebraElementStorage

}  // namespace mana
