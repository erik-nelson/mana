#include <Eigen/Dense>

#include "gtest/gtest.h"
// #include "lie/lie_algebra_element.h"
// #include "lie/lie_group_element.h"

namespace mana {

#if 0

// This test implements the unit complex numbers S1 as an example of a (abelian)
// Lie group.
class S1Element;

// Manifold traits specialization.
template <>
struct ManifoldTraits<S1Element> {
  using Element = S1Element;
  using Scalar = double;
  using Chart = ManifoldChart<S1Element>;
  using Geodesic = ManifoldGeodesic<S1Element>;
  using TangentVector = Eigen::Vector<Scalar, 1>;
  using EmbeddingPoint = Eigen::Vector<Scalar, 2>;
  static constexpr int Dimension = 1;
  static constexpr int EmbeddingDimension = 2;
};

// Algebra traits specialization.
template <>
struct AlgebraTraits<S1Element> {
  using Element = S1Element;
  using Scalar = double;
  using Vector = Eigen::Vector<Scalar, 1>;
  static constexpr int Dimension = 1;
};

class S1GroupElement : public LieGroupElement<S1Element> {
  using Base = LieGroupElement<S1Element>;

 public:
  using Scalar = typename Base::Scalar;
  using TangentVector = typename Base::TangentVector;
  using EmbeddingPoint = typename Base::EmbeddingPoint;
  using Jacobian = typename Base::Jacobian;
  using Chart = typename Base::Chart;
  using Geodesic = typename Base::Geodesic;
  using GroupElement = typename Base::GroupElement;
  using AlgebraElement = typename Base::AlgebraElement;
  static constexpr int Dimension = Base::Dimension;
  static constexpr int EmbeddingDimension = Base::EmbeddingDimension;

  // Implement `GroupElement` interface.
  static S1GroupElement IdentityImpl();
  S1GroupElement Inverse() const;
  S1GroupElement Compose(const S1GroupElement& rhs) const;

  // Implement `ManifoldElement` interface.
  static S1GroupElement ProjectImpl(const EmbeddingPoint& point);
  static bool IsValidImpl(const EmbeddingPoint& point);
  EmbeddingPoint PointImpl() const;
  Scalar DistanceToImpl(const S1GroupElement& rhs) const;
  S1GroupElement InterpolateImpl(const S1GroupElement& rhs,
                                 Scalar fraction) const;

  // Implement `LieGroupElement` interface.
  AlgebraElement logImpl() const;
  TangentVector LogImpl() const;
  static S1GroupElement expImpl(const AlgebraElement& algebra);
  static S1GroupElement ExpImpl(const TangentVector& coordinate);
  Jacobian AdjointImpl() const;
};

class S1AlgebraElement : public LieAlgebraElement<S1Element> {
 public:
};

TEST(Foo, Bar) {}
#endif

}  // namespace mana