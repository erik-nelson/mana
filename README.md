# Mana 
(idk, it'll have a manifold class and that name sounds cool)

Just screwing around... we will see where this gets.

I wanted to make it easier to build and run least squares problems, and also make a library that can easily represent 
continuous time state estimation problems via cubic hermite splines. Right now most libraries I know force you define 
factor and variable types. I want something more like the pytorch interface, where the computation graph is defined 
implicitly, and the user doesn't need to define factors with their Jacobians unless it's a highly nontrivial computation.

For example, I want something like the following interface for a point registration problem:

```cpp
// Initialize a spline representing a trajectory. Add a couple of knots.
CubicHermiteSpline<SE3GroupElement> spline;

// First knot point created from a transform and velocity variable.
Variable<SE3GroupElement> pose1 = ...;
Variable<Vector6> velocity1 = ...;
spline.AddKnot(timestamp1, pose1, velocity1);

// Second knot point created from a transform and velocity variable.
Variable<SE3GroupElement> pose2 = ...;
Variable<Vector6> velocity2 = ...;
spline.AddKnot(timestamp2, pose2, velocity2);

// Get an expression for the transform at a given timestamp. This transform
// won't be evaluated until later, but is a function of `pose{1,2}`, and 
// `velocity{1,2}`.
Expression<SE3GroupElement> Tx_world_body = spline.At(timestamp);

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

I suppose this is pretty similar to pytorch, except it's in C++ and is mainly created to handle least squares problems: the loss/cost function is always of the form $$\rho(||r(x(t))||_2)$$, where $$x : R --> X$$ is some continuous function mapping a time interval to elements of some Lie group $$X$$. I have no delusions about not being able to implement a huge library like pytorch myself as a hobby. 
This is also pretty similar to symforce, except that the computation graph can be built and run dynamically, and does not need to be compiled.
That means it'll be slow. But it'll look really slick too 8-)
