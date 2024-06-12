#include <Eigen/Dense>

#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "lie/lie_algebra_element.h"
#include "lie/lie_group_element.h"
#include "lie/test_utils.h"

namespace mana {

// The following tests use the R1GroupElement and R1AlgebraElement types from
// `test_utils.h`. This is a very simple Lie group implementation representing
// the real numbers under addition.

TEST(LieGroupElement, Construction) {
  // Default constructor yields identity.
  {
    R1GroupElement<double> value;
    EXPECT_EQ(value.Value(), 0);
  }

  // Construct from scalar.
  {
    R1GroupElement<double> value(3);
    EXPECT_EQ(value.Value(), 3);
  }

  // Construct from underlying storage type.
  {
    R1GroupElement<double> value(
        Traits<R1GroupElement<double>>::GroupElementStorage::Constant(5));
    EXPECT_EQ(value.Value(), 5);
  }
}

TEST(LieGroupElement, Identity) {
  const R1GroupElement<double> identity = R1GroupElement<double>::Identity();
  EXPECT_EQ(identity.Value(), 0);
}

TEST(LieGroupElement, Inverse) {
  R1GroupElement<double> value(3);

  // Inverting yields expected value (negation, via addition).
  const R1GroupElement<double> value_inv = value.Inverse();
  EXPECT_EQ(value_inv.Value(), -3);

  // Inverting twice yields the original value.
  const R1GroupElement<double> value_inv2 = value_inv.Inverse();
  EXPECT_EQ(value_inv2.Value(), 3);
  EXPECT_EQ(value, value_inv2);

  // Inverting in place.
  value.Invert();
  EXPECT_EQ(value.Value(), -3);
  EXPECT_EQ(value, value_inv);

  // Inverting in place twice.
  value.Invert();
  EXPECT_EQ(value.Value(), 3);
  EXPECT_EQ(value, value_inv2);
}

TEST(LieGroupElement, Log) {
  R1GroupElement<double> value(3);

  // (Upper-case) Log yields a Lie algebra coordinate vector.
  R1AlgebraElement<double> log1 = value.Log();
  EXPECT_EQ(log1.Value(), 3);

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
    EXPECT_EQ(group.Value(), 3);
  }
  // (Lower-case) exp creates a group element from a Lie algebra element.
  {
    const R1AlgebraElement<double>::AlgebraElementStorage algebra(3);
    const R1GroupElement<double> group = R1GroupElement<double>::exp(algebra);
    EXPECT_EQ(group.Value(), 3);
  }
}

TEST(LieGroupElement, Compose) {
  // Check that we can compose two elements.
  const R1GroupElement<double> value1(3);
  const R1GroupElement<double> value2(5);
  const R1GroupElement<double> value3 = value1.Compose(value2);
  EXPECT_EQ(value3.Value(), 8);

  // R1 is Abelian --> same as value2 * value1.
  const R1GroupElement<double> value4 = value2.Compose(value1);
  EXPECT_EQ(value4.Value(), 8);
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
  EXPECT_EQ(group1.RightMinus(group2), R1AlgebraElement<double>(-2));
  EXPECT_EQ(group1.LeftMinus(group2), R1AlgebraElement<double>(2));
  EXPECT_EQ(group1.Minus(group2), R1AlgebraElement<double>(-2));
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
  R1GroupElement<double> group1(1);
  R1GroupElement<double> group2(1.1);
  R1GroupElement<double> group3(1.01);

  constexpr double kEpsilon = Constants<double>::kEpsilon;
  EXPECT_TRUE(group1.EqualTo(group2, /*tolerance=*/1.00 + kEpsilon));
  EXPECT_TRUE(group1.EqualTo(group2, /*tolerance=*/0.10 + kEpsilon));
  EXPECT_FALSE(group1.EqualTo(group2, /*tolerance=*/0.01 + kEpsilon));

  EXPECT_TRUE(group1.EqualTo(group3, /*tolerance=*/1.00 + kEpsilon));
  EXPECT_TRUE(group1.EqualTo(group3, /*tolerance=*/0.10 + kEpsilon));
  EXPECT_TRUE(group1.EqualTo(group3, /*tolerance=*/0.01 + kEpsilon));
}

TEST(LieGroupElement, Operators) {
  R1GroupElement<double> a(3);
  R1GroupElement<double> b(5);
  R1AlgebraElement<double> c(8);

  // Right-plus addition.
  R1GroupElement<double> d = a + c;
  EXPECT_EQ(d, R1GroupElement<double>(11));

  // In-place addition.
  a += c;
  EXPECT_EQ(a, d);

  // Right-minus subtraction.
  R1AlgebraElement<double> e = a - b;
  EXPECT_EQ(e, R1AlgebraElement<double>(-6));

  // Composition.
  R1GroupElement<double> f = a * b;
  EXPECT_EQ(f, R1GroupElement<double>(16));

  // In-place composition.
  a *= b.Inverse();
  EXPECT_EQ(a, R1GroupElement<double>(6));

  // Equality, inequality.
  EXPECT_FALSE(a == b);
  EXPECT_TRUE(a == R1GroupElement<double>(6));
  EXPECT_TRUE(a != b);
  EXPECT_FALSE(a != R1GroupElement<double>(6));

  // Multiplication with Lie-algebra elements. This should perform group
  // composition, not actual multiplication.
  R1GroupElement<double>::AlgebraElementStorage g(5);
  R1GroupElement<double>::AlgebraElementStorage h = a * g;
  EXPECT_EQ(h[0], 11);
  R1GroupElement<double>::AlgebraElementStorage i = g * a;
  EXPECT_EQ(i, h);

  // Ostreamable.
  std::stringstream ss;
  ss << a;
  EXPECT_FALSE(ss.str().empty());
}

TEST(LieAlgebraElement, Construction) {
  // Default constructor yields identity.
  {
    R1AlgebraElement<double> algebra;
    EXPECT_EQ(algebra.Value(), 0);
  }

  // Construct from scalar.
  {
    R1AlgebraElement<double> algebra(3);
    EXPECT_EQ(algebra.Value(), 3);
  }

  // Construct from underlying storage type.
  {
    R1AlgebraElement<double> algebra(Eigen::Matrix<double, 1, 1>::Constant(5));
    EXPECT_EQ(algebra.Value(), 5);
  }
}

TEST(LieAlgebraElement, Identity) {
  const R1AlgebraElement<double> identity =
      R1AlgebraElement<double>::Identity();
  EXPECT_EQ(identity.Value(), 0);
}

TEST(LieAlgebraElement, Hat) {
  const R1AlgebraElement<double> algebra(5);
  const Eigen::Matrix<double, 1, 1> hat = algebra.Hat();
  EXPECT_EQ(hat[0], 5);
}

TEST(LieAlgebraElement, Vee) {
  const Eigen::Matrix<double, 1, 1> element(5);
  const R1AlgebraElement<double> algebra =
      R1AlgebraElement<double>::Vee(element);
  EXPECT_EQ(algebra.Value(), 5);
}

TEST(LieAlgebraElement, Exp) {
  const R1AlgebraElement<double> algebra(5);
  const R1GroupElement<double> group = algebra.Exp();
  EXPECT_EQ(group.Value(), 5);
}

TEST(LieAlgebraElement, Jacobians) {
  const R1AlgebraElement<double> algebra(5);
  const Eigen::Matrix<double, 1, 1> J_left = algebra.LeftJacobian();
  EXPECT_EQ(J_left[0], 1);
  const Eigen::Matrix<double, 1, 1> J_right = algebra.RightJacobian();
  EXPECT_EQ(J_right[0], 1);
}

TEST(LieAlgebraElement, Storage) {
  R1AlgebraElement<double> algebra(5);
  EXPECT_EQ(algebra.Storage()[0], 5);
  algebra.Storage()[0] = 3;
  EXPECT_EQ(algebra.Storage()[0], 3);
  EXPECT_EQ(algebra.At(0), 3);
  algebra.At(0) = 5;
  EXPECT_EQ(algebra.At(0), 5);
}

TEST(LieAlgebraElement, Log) {
  const R1GroupElement<double> group(5);
  const R1AlgebraElement<double> algebra = R1AlgebraElement<double>::Log(group);
  EXPECT_EQ(algebra.Value(), 5);
}

TEST(LieAlgebraElement, Compose) {
  const R1AlgebraElement<double> algebra1(5);
  const R1AlgebraElement<double> algebra2(8);

  // R1 is abelian - Lie bracket is zero.
  const R1AlgebraElement<double> algebra3 = algebra1.Compose(algebra2);
  EXPECT_EQ(algebra3.Value(), 0);
}

TEST(LieAlgebraElement, Equal) {
  R1AlgebraElement<double> algebra1(1);
  R1AlgebraElement<double> algebra2(1.1);
  R1AlgebraElement<double> algebra3(1.01);

  constexpr double kEpsilon = Constants<double>::kEpsilon;
  EXPECT_TRUE(algebra1.EqualTo(algebra2, /*tolerance=*/1.00 + kEpsilon));
  EXPECT_TRUE(algebra1.EqualTo(algebra2, /*tolerance=*/0.10 + kEpsilon));
  EXPECT_FALSE(algebra1.EqualTo(algebra2, /*tolerance=*/0.01 + kEpsilon));

  EXPECT_TRUE(algebra1.EqualTo(algebra3, /*tolerance=*/1.00 + kEpsilon));
  EXPECT_TRUE(algebra1.EqualTo(algebra3, /*tolerance=*/0.10 + kEpsilon));
  EXPECT_TRUE(algebra1.EqualTo(algebra3, /*tolerance=*/0.01 + kEpsilon));
}

TEST(LieAlgebraElement, Operators) {
  R1AlgebraElement<double> a(5);
  R1AlgebraElement<double> b(3);

  // Negation.
  EXPECT_EQ(-a, R1AlgebraElement<double>(-5));

  // Coordinate vector addition.
  EXPECT_EQ(a + b, R1AlgebraElement<double>(8));

  // In-place coordinate vector addition.
  a += b;
  EXPECT_EQ(a, R1AlgebraElement<double>(8));

  // Coordinate vector subtraction.
  EXPECT_EQ(a - b, R1AlgebraElement<double>(5));

  // In-place coordinate vector subtraction.
  a -= b;
  EXPECT_EQ(a, R1AlgebraElement<double>(5));

  // Coordinate vector scaling.
  EXPECT_EQ(a * 2.0, R1AlgebraElement<double>(10));
  EXPECT_EQ(a / 2.0, R1AlgebraElement<double>(2.5));

  // In-place coordinate vector scaling.
  a *= 2.0;
  EXPECT_EQ(a, R1AlgebraElement<double>(10));
  a /= 2.0;
  EXPECT_EQ(a, R1AlgebraElement<double>(5));

  // Lie bracket.
  EXPECT_EQ(a * b, R1AlgebraElement<double>::Identity());

  // In-place Lie bracket.
  a *= b;
  EXPECT_EQ(a, R1AlgebraElement<double>::Identity());

  // Coordinate access.
  EXPECT_EQ(a[0], 0.0);
  a[0] = 5.0;
  EXPECT_EQ(a[0], 5.0);

  // Equality, inequality.
  EXPECT_FALSE(a == b);
  EXPECT_TRUE(a == R1AlgebraElement<double>(5));
  EXPECT_TRUE(a != b);
  EXPECT_FALSE(a != R1AlgebraElement<double>(5));

  // Ostreamable.
  std::stringstream ss;
  ss << a;
  EXPECT_FALSE(ss.str().empty());

  // Multiplication with adjoint / jacobian matrices.
  const typename Traits<R1AlgebraElement<double>>::Jacobian jacobian(5);
  const R1AlgebraElement<double> c = jacobian * b;
  EXPECT_EQ(c, R1AlgebraElement<double>(15.0));
}

}  // namespace mana