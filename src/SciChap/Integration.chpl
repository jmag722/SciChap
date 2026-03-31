/*
  Integration: Basic numerical integration techniques.
*/
module Integration {
  import Math;
  import SciChap.Array;

  /*
  Performs trapezoidal rule, accounting for nonuniform domain spacing

  *References:*

  "Trapezoidal Rule: Numerical Implementation: Non-uniform grid". Wikipedia.

  :arg y: function values
  :type y: [?D] real

  :arg x: function input points
  :type x: [D] real

  :returns: Integral of y for x domain
  :rtype: real

  :throws Error: Domain `x` not strictly increasing or decreasing
  */
  proc trapezoid(const ref y: [?D] real, const ref x: [D] real): real throws
                 where D.rank == 1 {
    if x.size < 2 {
      return 0.0;
    }
    if !Array.isMonotonic(x) {
      throw new owned Error("Domain not strictly increasing or decreasing.");
    }
    const N: int = D.size;
    const ref x0 = x.reindex(0..<N);
    const ref y0 = y.reindex(0..<N);
    var sum: real = 0.0;
    sum += y0[0] * (x0[1] - x0[0]);
    sum += y0[N - 1] * (x0[N - 1] - x0[N - 2]);
    forall i in 1..<N-1 with (+ reduce sum) {
      sum += y0[i] * (x0[i + 1] - x0[i - 1]);
    }
    return 0.5 * sum;
  }

  /*
  Performs simpson rule, accounting for nonuniform domain spacing.
  If number of input points are less than 3, uses the trapezoidal rule.

  *References:*

  Cartwright, Kenneth V. (September 2017). "Simpson's Rule Cumulative
  Integration with MS Excel and Irregularly-spaced Data". Journal of
  Mathematical Sciences and Mathematics Education. 12 (2): 1–9.

  "Simpson's Rule: Composite Simpson's rule for irregularly spaced data".
  Wikipedia.

  :arg y: function values
  :type y: [?D] real

  :arg x: function input points
  :type x: [D] real

  :returns: Integral of y for x domain
  :rtype: real

  :throws Error: Domain `x` not strictly increasing or decreasing
  */
  proc simpson(const ref y: [?D] real, const ref x: [D] real): real throws
               where D.rank == 1 {
    if x.size < 3 {
      return trapezoid(y=y, x=x);
    } else {
      if !Array.isMonotonic(x) {
        throw new owned Error("Domain not strictly increasing or decreasing.");
      }
    }

    // algorithm below depends on 0-based indexing
    const ref x0 = x.reindex(0..<D.size);
    const ref y0 = y.reindex(0..<D.size);

    const h: [0..<D.size - 1] real = Array.diff(x0);
    const N: int = h.size;  // number of intervals
    var sum: real = 0.0;
    forall i in 0..<N/2 with (+ reduce sum) {
      sum += (h[2*i] + h[2*i+1]) / 6.0 * (
             (2.0 - h[2*i+1] / h[2*i]) * y0[2*i]
           + (h[2*i] + h[2*i+1])**2 / (h[2*i] * h[2*i+1]) * y0[2*i+1]
           + (2.0 - h[2*i] / h[2*i+1]) * y0[2*i+1 + 1]);
    }

    if Math.mod(N, 2) == 1 {
      const alpha: real = (2 * h[N - 1]**2 + 3 * h[N - 1] * h[N - 2])
                        / (6.0 * (h[N - 2] + h[N - 1]));
      const beta: real = (h[N - 1]**2 + 3 * h[N - 1] * h[N - 2])
                       / (6.0 * h[N - 2]);
      const eta: real = h[N - 1]**3 / (6.0 * h[N - 2] * (h[N - 2] + h[N - 1]));
      sum += alpha * y0[N] + beta * y0[N - 1] - eta * y0[N - 2];
    }
    return sum;
  }

}