#include <Eigen/Dense>

#include "lie/lie_algebra_element.h"
#include "lie/lie_group_element.h"

namespace mana {

// A basic Lie group element derived type representing the reals under
// addition. The corresponding Lie algebra is the vector space of reals
// under addition. Note that the reals actually form two different Lie groups -
// the group can be taken as the reals (minus 0) under multiplication, or the
// reals under addition. Both Lie groups have the same Lie algebra: the reals
// (as a 1D vector space).
template <typename ScalarT>
class R1GroupElement : public LieGroupElementBase<R1GroupElement<ScalarT>> {
  // Inherit traits from base class.
  using Base = LieGroupElementBase<R1GroupElement<ScalarT>>;

 public:
  MANA_INHERIT_TRAITS(Base);

  // Construct from scalar.
  explicit R1GroupElement(ScalarT value)
      : value_(GroupElementStorage::Constant(value)) {}

  // Construct from 1x1 matrix.
  explicit R1GroupElement(
      GroupElementStorage value = GroupElementStorage::Constant(0))
      : value_(std::move(value)) {}

  // Override base API --------------------------------------------------------
  static GroupElement Identity() { return GroupElement(0); }
  GroupElement Inverse() const { return GroupElement(-value_); }
  void Invert() { value_ *= ScalarT{-1}; }
  AlgebraElement Log() const { return AlgebraElement(value_); }
  GroupElement Compose(const GroupElement& rhs) const {
    return GroupElement(value_ + rhs.value_);
  }
  Jacobian Adjoint() const { return Jacobian::Constant(1); }
  GroupElementStorage& Storage() { return value_; }
  const GroupElementStorage& Storage() const { return value_; }

  // Helper to access single underlying element.
  const ScalarT& Value() const { return value_[0]; }

 private:
  GroupElementStorage value_;
};

// The corresponding Lie algebra element derived type representing the reals
// under addition.
template <typename ScalarT>
class R1AlgebraElement
    : public LieAlgebraElementBase<R1AlgebraElement<ScalarT>> {
  using Base = LieAlgebraElementBase<R1AlgebraElement<ScalarT>>;

 public:
  // Inherit traits from base class.
  MANA_INHERIT_TRAITS(Base);

  // Construct from scalar.
  explicit R1AlgebraElement(ScalarT value)
      : value_(AlgebraVectorStorage::Constant(value)) {}

  // Construct from 1x1 matrix.
  explicit R1AlgebraElement(
      AlgebraVectorStorage value = AlgebraVectorStorage::Constant(0))
      : value_(std::move(value)) {}

  // Override base API.
  static AlgebraElement Identity() { return AlgebraElement(0); }
  AlgebraElementStorage Hat() const { return AlgebraElementStorage(value_); }
  static AlgebraElement Vee(const AlgebraElementStorage& element) {
    return AlgebraElement(element);
  }
  GroupElement Exp() const { return GroupElement(value_); }
  Jacobian RightJacobian() const { return Jacobian(1); }
  Jacobian LeftJacobian() const { return Jacobian(1); }
  AlgebraVectorStorage& Storage() { return value_; }
  const AlgebraVectorStorage& Storage() const { return value_; }

  // Helper to access single underlying element.
  const ScalarT& Value() const { return value_[0]; }

 private:
  AlgebraVectorStorage value_;
};

// Define traits for our `R1GroupElement` type.
template <typename ScalarT>
struct Traits<R1GroupElement<ScalarT>> {
  static constexpr size_t Dims = 1;
  static constexpr size_t Dofs = 1;
  using Scalar = ScalarT;
  using Jacobian = Eigen::Matrix<ScalarT, 1, 1>;
  using GroupElement = R1GroupElement<ScalarT>;
  using AlgebraElement = R1AlgebraElement<ScalarT>;
  using GroupElementStorage = Eigen::Vector<ScalarT, 1>;
  using AlgebraVectorStorage = Eigen::Vector<ScalarT, 1>;
  using AlgebraElementStorage = Eigen::Vector<ScalarT, 1>;
};

// Copy traits to `RealsUnderAdditionAlgebraElement` type.
template <typename ScalarT>
struct Traits<R1AlgebraElement<ScalarT>>
    : public Traits<R1GroupElement<ScalarT>> {};

}  // namespace mana