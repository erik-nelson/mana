#pragma once

namespace mana {

template <typename ScalarT>
struct Constants {
  static constexpr ScalarT kEps = ScalarT{1e-8};
};

}  // namespace mana