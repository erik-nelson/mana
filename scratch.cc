// Want something like Pytorch, that allows "backpropping" (in a quadratic /
// least-squares fashion) loss through to poses. Example: Point cloud
// registration cost function.
//   sum_i || p_i - (R q_i + t) ||^2

Variable<SE3Groupd> Tx_p_q = SE3Groupd::Identity();
Expression<double> cost = 0;

for (int i = 0; i < num_points; ++i) {
  const Vector3d pi = p[i];
  const Vector3d qi = q[i];
  cost += (pi - Tx_p_q * qi).norm2();
}

cost.Optimize(/*options=...*/);

// `Tx_p_q` automagically updated to local minimum. Optimization by default is
// quadratic. All the various Jacobians related to Lie groups are implemented
// already.