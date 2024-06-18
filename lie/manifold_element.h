#pragma once

#include <array>
#include <cassert>

#include "utils/crtp.h"

namespace mana {

// Traits template for a manifold.
template <typename Derived>
struct ManifoldTraits {};

// Base CRTP class for an element of a (differentiable) manifold. Manifolds are
// expressed in an extrinsic fashion, i.e. they are assumed to be embedded
// within a larger ambient (vector) space. Derived manifold element types should
// inherit from this like so:
//
//   class S2Element : public ManifoldElement<S2Element> {
//     ...
//   };
//
// Derived classes must implement these methods.
// - static Derived ProjectImpl(const EmbeddingPoint& point);
// - static bool IsValidImpl(const EmbeddingPoint& point);
// - EmbeddingPoint PointImpl() const;
// - StorageType& StorageImpl();
// - const StorageType& StorageImpl() const;
// - std::array<TangentVector, Dimension> TangentSpaceBasisImpl() const;
// - Scalar DistanceToImpl(const Element& rhs) const;
// - Element InterpolateImpl(const Element& rhs, Scalar fraction) const;
//
template <typename Derived>
class ManifoldElement : public Crtp<Derived> {
 public:
  // An element of the manifold (identical to `Derived`).
  using Element = typename ManifoldTraits<Derived>::Element;
  // The scalar type used to represent distances.
  using Scalar = typename ManifoldTraits<Derived>::Scalar;
  // The chart type (see `ManifoldChart` below).
  using Chart = typename ManifoldTraits<Derived>::Chart;
  // The geodesic type (see `ManifoldGeodesic` below).
  using Geodesic = typename ManifoldTraits<Derived>::Geodesic;
  // A vector in a tangent space of this manifold.
  using TangentVector = typename ManifoldTraits<Derived>::TangentVector;
  // A point in the embedding space that this manifold lies in. Not necessarily
  // a vector. For instance, points in the embedding space of 3D rotations SO(3)
  // are 3x3 matrices.
  using EmbeddingPoint = typename ManifoldTraits<Derived>::EmbeddingPoint;
  // The underlying storage type we use to implement this element. For example,
  // for 3D rotations this might be a 4x1 quaternion vector or a 3x1 axis-angle
  // vector to reduce storage size.
  using StorageType = typename ManifoldTraits<Derived>::StorageType;
  // The dimension of the manifold.
  static constexpr int Dimension = ManifoldTraits<Derived>::Dimension;
  // The dimension of the embedding space.
  static constexpr int EmbeddingDimension =
      ManifoldTraits<Derived>::EmbeddingDimension;

  // Construct by projecting from a point in the embedding space onto the
  // manifold. For example, an SO(3) derived manifold element type might
  // implement this by projecting a 3x3 matrix onto {R | R'R=I, det(R)=1}, the
  // space of valid rotation matrices.
  static Element Project(const EmbeddingPoint& point) {
    return Derived::ProjectImpl(point);
  }

  // Checks if this point in embedding space lies on the manifold. In general
  // this is true if the underlying `StorageType` is coordinates on the
  // manifold. However if the `StorageType` is the same as a point in embedding
  // space, elements of this class may not always be valid manifold elements.
  static bool IsValid(const EmbeddingPoint& point) {
    return Derived::IsValidImpl(point);
  }

  // Return this element's point in the underlying embedding space.
  EmbeddingPoint Point() const { return Crtp<Derived>::get().PointImpl(); }

  // Retrieve the underlying storage for this manifold element.
  StorageType& Storage() { return Crtp<Derived>::get().StorageImpl(); }
  const StorageType& Storage() const {
    return Crtp<Derived>::get().StorageImpl();
  }

#if 0
  // Build a chart at this point on the manifold.
  Chart LocalChart() const { return Chart(*this); }

  // Returns a basis for the tangent space at this point on the manifold.
  std::array<TangentVector, Dimension> TangentSpaceBasis() const {
    return Crtp<Derived>::get().TangentSpaceBasisImpl();
  }
#endif

  // Builds a geodesic curve parameterized between two points on the manifold.
  Geodesic GeodesicTo(const Element& rhs) const {
    return Geodesic(Crtp<Derived>::get(), rhs);
  }

  // Compute the distance between two points on the manifold.
  Scalar DistanceTo(const Element& rhs) const {
    return Crtp<Derived>::get().DistanceToImpl(rhs);
  }

  // Interpolate along the path from this manifold element to `rhs`, following
  // the geodesic. The provided fraction should be in [0, 1] for points along
  // the geodesic, and outside of that range to perform extrapolation instead.
  Element Interpolate(const Element& rhs, Scalar fraction) const {
    return Crtp<Derived>::get().InterpolateImpl(rhs, fraction);
  }

  // Equality checking.
  bool operator==(const Element& rhs) const {
    return Storage() == rhs.Storage();
  }

  // Inequality checking.
  bool operator!=(const Element& rhs) const { return !(*this == rhs); }
};

// Class representing a chart on a manifold, mapping from the manifold to its
// tangent space.
template <typename Derived>
class ManifoldChart {
 public:
  using Element = typename ManifoldTraits<Derived>::Element;
  using TangentVector = typename ManifoldTraits<Derived>::TangentVector;

  // Construct from origin point on the manifold.
  explicit ManifoldChart(Element origin) : origin_(std::move(origin)) {}

  // Map an element on the manifold to a tangent vector using the chart's
  // forward map.
  TangentVector ToTangent(const Element& rhs) const {
    // TODO(erik): implement.
    return {};
  }

  // Map a tangent vector to an element on the manifold using the chart's
  // reverse map.
  Element ToManifold(const TangentVector& rhs) const {
    // TODO(erik): implement.
    return {};
  }

 private:
  // The point forming the origin of this chart. The zero tangent vector is
  // mapped to this point on the manifold.
  Element origin_;
};

template <typename Derived>
class ManifoldGeodesic {
 public:
  using Scalar = typename ManifoldTraits<Derived>::Scalar;
  using Element = typename ManifoldTraits<Derived>::Element;

  // Construct from start and end points.
  ManifoldGeodesic(Element beg, Element end)
      : beg_(std::move(beg)), end_(std::move(end)) {}

  // Return the geodesic's begin point.
  const Element& beg() const { return beg_; }

  // Return the geodesic's end point.
  const Element& end() const { return end_; }

  // Interpolate along the geodesic at the provided fraction. Values in [0, 1]
  // perform true interpolation, values outside of this range perform
  // extrapolation.
  Element Interpolate(Scalar fraction) const {
    return beg_.Interpolate(end_, fraction);
  }

  // Return the length of this geodesic.
  Scalar Length() const { return beg_.DistanceTo(end_); }

 private:
  // The geodesic's start and end points.
  Element beg_, end_;
};

}  // namespace mana