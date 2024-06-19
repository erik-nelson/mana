#include <Eigen/Dense>
#include <cmath>
#include <iostream>

#include "gtest/gtest.h"
#include "lie/constants.h"
#include "lie/manifold_element.h"
#include "utils/angles.h"

namespace mana {

// This test implements the 2D sphere, S2, as an example manifold to test the
// interfaces of the `ManifoldElement<>`, `ManifoldChart<>`, and
// `ManifoldGeodesic<>` template types.
class S2Element;

// Trait specialization for elements of S2.
template <>
struct ManifoldTraits<S2Element> {
  using Element = S2Element;
  using Scalar = double;
  using Chart = ManifoldChart<S2Element>;
  using Geodesic = ManifoldGeodesic<S2Element>;
  // Tangent planes are copies of R^2.
  using TangentVector = Eigen::Vector<Scalar, 2>;
  // The 2D sphere is embedded in R^3.
  using EmbeddingPoint = Eigen::Vector<Scalar, 3>;
  // The manifold itself is 2 dimensional.
  static constexpr int Dimension = 2;
  // The embedding space is 3 dimensional.
  static constexpr int EmbeddingDimension = 3;
};

class S2Element : public ManifoldElement<S2Element> {
 public:
  using Element = typename ManifoldTraits<S2Element>::Element;
  using Scalar = typename ManifoldTraits<S2Element>::Scalar;
  using Chart = typename ManifoldTraits<S2Element>::Chart;
  using Geodesic = typename ManifoldTraits<S2Element>::Geodesic;
  using TangentVector = typename ManifoldTraits<S2Element>::TangentVector;
  using EmbeddingPoint = typename ManifoldTraits<S2Element>::EmbeddingPoint;
  static constexpr int Dimension = ManifoldTraits<S2Element>::Dimension;
  static constexpr int EmbeddingDimension =
      ManifoldTraits<S2Element>::EmbeddingDimension;

  // Construct from angles.
  S2Element(Scalar phi, Scalar theta)
      : phi_(Normalized(phi)), theta_(Normalized(theta)) {}

  // Helpers to get phi and theta from underlying storage.
  double& phi() { return phi_; }
  const double& phi() const { return phi_; }
  double& theta() { return theta_; }
  const double& theta() const { return theta_; }

  // Stream printing.
  friend std::ostream& operator<<(std::ostream& os, const Element& element) {
    return (os << "Phi, theta: (" << element.phi() << ", " << element.theta()
               << ")\n");
  }

  // Implement `ManifoldElement` interface.
  static S2Element ProjectImpl(const EmbeddingPoint& point) {
    const Scalar magnitude = point.norm();
    assert(magnitude != 0.0);
    const Scalar phi = std::acos(point.z() / magnitude);
    const Scalar theta = std::atan2(point.y(), point.x());
    return S2Element(phi, theta);
  }

  static bool IsValidImpl(const EmbeddingPoint& point) {
    return (point.squaredNorm() - 1) < Constants<Scalar>::kEpsilon;
  }

  EmbeddingPoint PointImpl() const {
    EmbeddingPoint result;
    SphericalToCartesian(/*rho=*/1.0, theta(), phi(),  //
                         result.x(), result.y(), result.z());
    return result;
  }

#if 0
  std::array<TangentVector, 1> TangentSpaceBasisImpl() const { return {}; }
#endif

  Scalar DistanceToImpl(const S2Element& rhs) const {
    return std::acos(Point().dot(rhs.Point()));
  }

  S2Element InterpolateImpl(const S2Element& rhs, Scalar fraction) const {
    // Slerp between *this and rhs.
    const EmbeddingPoint p0 = Point();
    const EmbeddingPoint p1 = rhs.Point();
    if (p0 == p1) return rhs;  // nothing to interpolate.
    const double angle = std::acos(p0.dot(p1));
    const double sin_angle_interp1 = std::sin((1.0 - fraction) * angle);
    const double sin_angle_interp2 = std::sin(fraction * angle);
    const double sin_angle = std::sin(angle);
    return S2Element::Project((1.0 / sin_angle) * sin_angle_interp1 * p0 +
                              sin_angle_interp2 * p1);
  }

 private:
  Scalar phi_, theta_;
};

TEST(ManifoldElement, Project) {
  constexpr double eps = Constants<double>::kEpsilon;
  S2Element element = S2Element::Project(Eigen::Vector3d(3, 3, 0));
  EXPECT_NEAR(element.phi(), DegToRad(90.0), eps);
  EXPECT_NEAR(element.theta(), DegToRad(45.0), eps);

  element = S2Element::Project(Eigen::Vector3d(0.1, 0, 0));
  EXPECT_NEAR(element.phi(), DegToRad(90.0), eps);
  EXPECT_NEAR(element.theta(), DegToRad(0.0), eps);

  element = S2Element::Project(Eigen::Vector3d(-1, 0, 0));
  EXPECT_NEAR(element.phi(), DegToRad(90.0), eps);
  EXPECT_NEAR(element.theta(), DegToRad(180.0), eps);

  element = S2Element::Project(Eigen::Vector3d(1, 0, 1));
  EXPECT_NEAR(element.phi(), DegToRad(45.0), eps);
  EXPECT_NEAR(element.theta(), DegToRad(0.0), eps);
}

TEST(ManifoldElement, IsValid) {
  constexpr double eps = Constants<double>::kEpsilon;
  EXPECT_FALSE(S2Element::IsValid(Eigen::Vector3d(1, 1, 1)));
  EXPECT_FALSE(S2Element::IsValid(Eigen::Vector3d(1 + eps, 0, 0)));

  double rho = 1.0;
  for (double phi = 0; phi < M_PI; phi += 0.01) {
    for (double theta = 0; theta < 2 * M_PI; theta += 0.01) {
      Eigen::Vector3d unit;
      SphericalToCartesian(rho, theta, phi,  //
                           unit.x(), unit.y(), unit.z());
      EXPECT_TRUE(S2Element::IsValid(unit));
    }
  }

  rho = 5.0;
  for (double phi = 0; phi < M_PI; phi += 0.01) {
    for (double theta = 0; theta < 2 * M_PI; theta += 0.01) {
      Eigen::Vector3d non_unit;
      SphericalToCartesian(rho, theta, phi,  //
                           non_unit.x(), non_unit.y(), non_unit.z());
      EXPECT_FALSE(S2Element::IsValid(non_unit));
    }
  }
}

TEST(ManifoldElement, Point) {
  S2Element element(/*phi=*/DegToRad(45.0), /*theta=*/0);

  // Get the element as a point in the manifold's embedding space.
  const Eigen::Vector3d point_act = element.Point();
  const Eigen::Vector3d point_exp = (0.5 * M_SQRT2) * Eigen::Vector3d(1, 0, 1);
  EXPECT_LT((point_act - point_exp).array().abs().maxCoeff(),
            Constants<double>::kEpsilon);
}

TEST(ManifoldElement, GeodesicTo) {
  S2Element a(/*phi=*/DegToRad(90), /*theta=*/0);
  S2Element b(/*phi=*/DegToRad(90), /*theta=*/DegToRad(90));

  // Build a geodesic between the two points.
  S2Element::Geodesic geodesic = a.GeodesicTo(b);
  EXPECT_EQ(geodesic.beg(), a);
  EXPECT_EQ(geodesic.end(), b);
  EXPECT_EQ(geodesic.Length(), DegToRad(90));
  EXPECT_EQ(geodesic.Interpolate(0), a);
  EXPECT_EQ(geodesic.Interpolate(1), b);
  EXPECT_EQ(geodesic.Interpolate(0.5),
            S2Element(/*phi=*/DegToRad(90), /*theta*/ DegToRad(45)));
}

TEST(ManifoldElement, DistanceTo) {
  S2Element a(/*phi=*/DegToRad(90), /*theta=*/0);
  S2Element b(/*phi=*/DegToRad(90), /*theta=*/DegToRad(90));

  // Check interface for the distance between the two points.

  // Symmetric.
  EXPECT_EQ(a.DistanceTo(b), DegToRad(90));
  EXPECT_EQ(b.DistanceTo(a), DegToRad(90));

  // Identity.
  EXPECT_EQ(a.DistanceTo(a), 0);
  EXPECT_EQ(b.DistanceTo(b), 0);

  // Triangle inequality.
  S2Element c(/*phi=*/DegToRad(45), /*theta=*/DegToRad(-45));
  EXPECT_GE(a.DistanceTo(b) + b.DistanceTo(c), a.DistanceTo(c));
  EXPECT_GE(a.DistanceTo(c) + c.DistanceTo(b), a.DistanceTo(b));
  EXPECT_GE(b.DistanceTo(a) + a.DistanceTo(c), b.DistanceTo(c));
}

TEST(ManifoldElement, EqualTo) {
  S2Element a(/*phi=*/DegToRad(90), /*theta=*/0);
  S2Element b = a;
  S2Element c(/*phi=*/DegToRad(90), /*theta=*/DegToRad(90));

  EXPECT_TRUE(a.EqualTo(b));
  EXPECT_TRUE(b.EqualTo(a));
  EXPECT_EQ(a, b);

  EXPECT_FALSE(a.EqualTo(c));
  EXPECT_FALSE(c.EqualTo(a));
  EXPECT_NE(a, c);

  const double tolerance = a.DistanceTo(c) + Constants<double>::kEpsilon;
  EXPECT_TRUE(a.EqualTo(c, tolerance));
}

}  // namespace mana