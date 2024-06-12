#include <Eigen/Dense>
#include <type_traits>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
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
  explicit R1AlgebraElement(ScalarT value = 0)
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

TEST(LieGroupElement, Construction) {
  using T = Traits<R1GroupElement<double>>;

  // Default constructor yields identity.
  {
    R1GroupElement<double> value;
    EXPECT_EQ(value.Storage(), T::GroupElementStorage::Constant(0));
  }

  // Construct from scalar.
  {
    R1GroupElement<double> value(3);
    EXPECT_EQ(value.Storage(), T::GroupElementStorage::Constant(3));
  }

  // Construct from underlying storage type.
  {
    R1GroupElement<double> value(T::GroupElementStorage::Constant(5));
    EXPECT_EQ(value.Storage(), T::GroupElementStorage::Constant(5));
  }
}

TEST(LieGroupElement, Identity) {
  const R1GroupElement<double> identity = R1GroupElement<double>::Identity();
  const Eigen::Matrix<double, 1, 1> identity_value = identity.Storage();
  EXPECT_EQ(identity_value[0], 0);
}

TEST(LieGroupElement, Inverse) {
  R1GroupElement<double> value(3);

  // Inverting yields expected value (negation, via addition).
  const R1GroupElement<double> value_inv = value.Inverse();
  EXPECT_EQ(value_inv.Storage()[0], -3);

  // Inverting twice yields the original value.
  const R1GroupElement<double> value_inv2 = value_inv.Inverse();
  EXPECT_EQ(value_inv2.Storage()[0], 3);
  EXPECT_EQ(value, value_inv2);

  // Inverting in place.
  value.Invert();
  EXPECT_EQ(value.Storage()[0], -3);
  EXPECT_EQ(value, value_inv);

  // Inverting in place twice.
  value.Invert();
  EXPECT_EQ(value.Storage()[0], 3);
  EXPECT_EQ(value, value_inv2);
}

TEST(LieGroupElement, Log) {
  R1GroupElement<double> value(3);

  // (Upper-case) Log yields a Lie algebra coordinate vector.
  R1AlgebraElement<double> log1 = value.Log();
  EXPECT_EQ(log1.At(0), 3);

  // (Lower-case) log yields a Lie algebra element.
  R1AlgebraElement<double>::AlgebraElementStorage log2 = value.log();
  EXPECT_EQ(log2[0], 3);
}

TEST(LieGroupElement, Exp) {
  // (Upper-case) Exp creates a group element from a Lie algebra coordinate
  // vector.
  {
    const R1AlgebraElement<double> algebra(3);
    const R1GroupElement<double> group = R1GroupElement<double>::Exp(algebra);
    EXPECT_EQ(group.Storage()[0], 3);
  }
  // (Lower-case) exp creates a group element from a Lie algebra element.
  {
    const R1AlgebraElement<double>::AlgebraElementStorage algebra(3);
    const R1GroupElement<double> group = R1GroupElement<double>::exp(algebra);
    EXPECT_EQ(group.Storage()[0], 3);
  }
}

TEST(LieGroupElement, Compose) {
  // Check that we can compose two elements.
  const R1GroupElement<double> value1(3);
  const R1GroupElement<double> value2(5);
  const R1GroupElement<double> value3 = value1.Compose(value2);
  EXPECT_EQ(value3.Storage()[0], 8);

  // R1 is Abelian --> same as value2 * value1.
  const R1GroupElement<double> value4 = value2.Compose(value1);
  EXPECT_EQ(value4.Storage()[0], 8);
}

TEST(LieGroupElement, Adjoint) {
  using T = Traits<R1GroupElement<double>>;

  const R1GroupElement<double> group(5);
  const T::Jacobian adjoint = group.Adjoint();
  EXPECT_EQ(adjoint[0], 1);

  // Check adjoint relationship:
  // Adj_{X} * v = vee(X * hat(v) * inv(X))
  //
  // Note that implicitly these operator*'s should be performing the addition
  // operation - they are group composition, not matrix multiplication.
  const R1AlgebraElement<double> algebra(3);
  const R1AlgebraElement<double> lhs = adjoint * algebra;
  const R1AlgebraElement<double> rhs =
      R1AlgebraElement<double>::Vee(group * algebra.Hat() * group.Inverse());
  EXPECT_EQ(lhs, rhs);
}

TEST(LieGroupElement, Plus) {
  const R1GroupElement<double> group(5);
  const R1AlgebraElement<double> algebra(3);

  // Check that the plus operators are implemented (commutatively) on both
  // the group and algebra classes.
  EXPECT_EQ(algebra.RightPlus(group), group.RightPlus(algebra));
  EXPECT_EQ(algebra.LeftPlus(group), group.LeftPlus(algebra));
  EXPECT_EQ(algebra.Plus(group), group.Plus(algebra));

  // Check values.All the same since R1 is commutative under addition.
  EXPECT_EQ(group.RightPlus(algebra), R1GroupElement<double>(8));
  EXPECT_EQ(group.LeftPlus(algebra), R1GroupElement<double>(8));
  EXPECT_EQ(group.Plus(algebra), R1GroupElement<double>(8));
}

TEST(LieGroupElement, Minus) {
  const R1GroupElement<double> group1(5);
  const R1GroupElement<double> group2(3);

  // Check that the minus operators are implemented (commutatively) on both
  // the group and algebra classes.
  EXPECT_EQ(R1AlgebraElement<double>::LeftMinus(group1, group2),
            group1.LeftMinus(group2));
  EXPECT_EQ(R1AlgebraElement<double>::RightMinus(group1, group2),
            group1.RightMinus(group2));
  EXPECT_EQ(R1AlgebraElement<double>::Minus(group1, group2),
            group1.Minus(group2));

  // Check values.
  EXPECT_EQ(group1.RightMinus(group2), R1AlgebraElement<double>(2));
  EXPECT_EQ(group1.LeftMinus(group2), R1AlgebraElement<double>(2));
  EXPECT_EQ(group1.Minus(group2), R1AlgebraElement<double>(2));
}

TEST(LieGroupElement, Between) {
  R1GroupElement<double> group1(5);
  R1GroupElement<double> group2(3);

  EXPECT_EQ(group1.BetweenInner(group2), R1GroupElement<double>(-2));
  EXPECT_EQ(group2.BetweenInner(group1), R1GroupElement<double>(2));
  EXPECT_EQ(group1.BetweenOuter(group2), R1GroupElement<double>(2));
  EXPECT_EQ(group2.BetweenOuter(group1), R1GroupElement<double>(-2));
}

TEST(LieGroupElement, Equal) {
  ;
  ;
}

TEST(LieGroupElement, Operators) {
  ;
  ;
}

TEST(LieAlgebraElement, Construction) {
  ;
  ;
}

}  // namespace mana