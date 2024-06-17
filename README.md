Just screwing around... we will see where this gets.

I wanted to make it easier to build and run least squares problems, and also make a library that can easily represent 
continuous time state estimation problems via cubic hermite splines. Right now most libraries I know force you define 
factor and variable types. I want something more like the pytorch interface, where the computation graph is defined 
implicitly, and the user doesn't need to define factors with their Jacobians unless it's a highly nontrivial computation.

For example, I want something like the following interface for a point registration problem:

```cpp
// Initialize a spline representing a trajectory. Add a couple of knots.
CubicHermiteSpline<SE3GroupElement> spline;
spline.AddKnot(timestamp1, pose1, velocity1);
spline.AddKnot(timestamp2, pose2, velocity2);

// Get an expression for the transform at a given timestamp. Note that this transform
// has not been evaluated yet.
Variable<SE3GroupElement> Tx_world_body = spline.At(timestamp);

// Constrain the transform via a point matching measurement.
// argmin_{Tx_wb} || p_i - (Rx_w_b * q_i + tx_w_b) ||^2_{sigma}
Expression<double> cost = 0;
for (size_t i = 0; i < num_points; ++i) {
  Vector3 pi = points_p[i];
  Vector3 qi = points_q[i];
  cost += (pi - Tx_world_body * qi).norm2();
}

GaussNewtonOptimizer optimizer(cost);
optimizer.Optimize();

// Updates propagated back to knot points of the original spline.
```

We will see how much of this interface is possible to implement, but I'd like to basically not need to ever implement 
a factor or Jacobian for simple problems like this, and instead model everything as an expression tree (where Jacobians 
are all implemented in the library already).
