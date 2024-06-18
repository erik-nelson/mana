#include <utility>

#include "gtest/gtest.h"
#include "lie/algebra_element.h"

namespace mana {

// This test implements the algebra of 3D vectors - 3D vectors form an algebra
// using the cross product as their product operator.
class Vector3;

// Trait specialization for elements of the Vector3 algebra.
template <>
struct AlgebraTraits<Vector3> {
  using Element = Vector3;
  using Scalar = double;
  using Vector = Eigen::Vector<Scalar, 3>;
  static constexpr int Dimension = 3;
};

class Vector3 : public AlgebraElement<Vector3> {
 public:
  using Element = typename AlgebraTraits<Vector3>::Element;
  using Scalar = typename AlgebraTraits<Vector3>::Scalar;
  using Vector = typename AlgebraTraits<Vector3>::Vector;
  static constexpr int Dimension = AlgebraTraits<Vector3>::Dimension;

  // Convenience: construct from 3 scalars.
  Vector3(Scalar x, Scalar y, Scalar z) : Vector3(Vector(x, y, z)) {}

  // Implement `AlgebraElement` interface.
  explicit Vector3(Vector coordinates)
      : AlgebraElement<Vector3>(std::move(coordinates)) {}

  Element ComposeImpl(const Element& rhs) const {
    // The algebra product in R^3 is the cross product.
    return Element(coordinates_.cross(rhs.coordinates_));
  }
};

TEST(AlgebraElement, Coordinates) {
  Vector3 a(1, 2, 3);
  EXPECT_EQ(a.Coordinates(), Eigen::Vector3d(1, 2, 3));
  a.Coordinates()(0) = 3;
  EXPECT_EQ(a.Coordinates(), Eigen::Vector3d(3, 2, 3));
}

TEST(AlgebraElement, HatVee) {
  Vector3 a = Vector3::Hat(Eigen::Vector3d(1, 2, 3));
  EXPECT_EQ(a.Coordinates(), Eigen::Vector3d(1, 2, 3));
  EXPECT_EQ(a.Vee(), a.Coordinates());
}

TEST(AlgebraElement, Compose) {
  Vector3 a(1, 2, 3);
  Vector3 b(4, 5, 6);
  EXPECT_EQ(a.Compose(b).Coordinates(),
            Eigen::Vector3d(1, 2, 3).cross(Eigen::Vector3d(4, 5, 6)));
}

TEST(AlgebraElement, Operators) {
  Vector3 a(1, 2, 3);
  Vector3 b(4, 5, 6);

  // Vector addition: a + b.
  const Vector3 c = a + b;
  EXPECT_EQ(c.Coordinates(), Eigen::Vector3d(5, 7, 9));

  // Vector in-place addition: a += b.
  a += b;
  EXPECT_EQ(a.Coordinates(), Eigen::Vector3d(5, 7, 9));

  // Vector subtraction: a - b.
  const Vector3 d = a - b;
  EXPECT_EQ(d.Coordinates(), Eigen::Vector3d(1, 2, 3));

  // Vector in-place subtraction: a -= b.
  a -= b;
  EXPECT_EQ(a.Coordinates(), Eigen::Vector3d(1, 2, 3));

  // Vector product: a * b.
  const Vector3 e = a * b;
  EXPECT_EQ(e.Coordinates(),
            Eigen::Vector3d(1, 2, 3).cross(Eigen::Vector3d(4, 5, 6)));

  // Vector negation: -a;
  const Vector3 f = -a;
  EXPECT_EQ(f.Coordinates(), Eigen::Vector3d(-1, -2, -3));

  // Vector-scalar multiplication: v * s.
  const Vector3 g = a * 3.0;
  EXPECT_EQ(g.Coordinates(), Eigen::Vector3d(3, 6, 9));

  // Vector-scalar in-place multiplication: v *= s.
  a *= 3.0;
  EXPECT_EQ(a.Coordinates(), Eigen::Vector3d(3, 6, 9));

  // Vector-scalar division: v / s.
  const Vector3 h = a / 3.0;
  EXPECT_EQ(h.Coordinates(), Eigen::Vector3d(1, 2, 3));

  // Vector-scalar in-place division: v /= s.
  a /= 3.0;
  EXPECT_EQ(a.Coordinates(), Eigen::Vector3d(1, 2, 3));

  // Scalar-vector multiplication: s * v.
  const Vector3 i = 3.0 * a;
  EXPECT_EQ(i.Coordinates(), Eigen::Vector3d(3, 6, 9));

  // Equality comparison.
  const Vector3 j = a;
  EXPECT_EQ(a, j);
  EXPECT_TRUE(a == j);

  // Inequality comparison.
  EXPECT_NE(a, i);
  EXPECT_TRUE(a != i);
}

}  // namespace mana