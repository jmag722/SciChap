/*
  Linear algebra methods beyond those found in the standard LinearAlgebra
  library.

*/
module LinAlg {

  import Math.fma;

  /*
    Serial tri-diagonal matrix algorithm (TDMA), or the Thomas algorithm.
    Matrix should be diagonally dominant (or close enough).

    The solution is computed in-place and stored at ``rhs``.

    Reference:
      Laszlo, E., Giles, M., & Appleyard, J. (2016).
      Manycore algorithms for batch scalar and block tridiagonal solvers.
      `ACM Transactions on Mathematical Software (TOMS)`, `42`(4), 1-36.

    :arg lower: lower diagonal (subdiagonal), length ``N-1``

    :arg diagonal: main diagonal, length ``N``

    :arg upper: upper diagonal (superdiagonal), length ``N-1``

    :arg rhs: right hand side of system, length ``N``

    :returns: void (solution variable stored at ``rhs``)
  */
  proc tdma(const lower: [?Dm1] real, const diagonal: [?D] real,
            ref upper: [Dm1] real, ref rhs: [D] real): void
            where Dm1.rank == 1 && D.rank == 1 {
    const N = rhs.size;
    // index lo starting at 1 to avoid "-1" everywhere below
    const ref lo = lower.reindex(1..N-1);
    const ref diag = diagonal.reindex(0..N-1);
    ref up = upper.reindex(0..N-2);
    ref x = rhs.reindex(0..N-1);

    // forward substitution
    x[0] /= diag[0];
    up[0] /= diag[0];
    var r: real;
    for idx in 1..N-2 {
      r = 1.0 / fma(-lo[idx], up[idx-1], diag[idx]);
      x[idx] = r * fma(-lo[idx], x[idx-1], x[idx]);
      up[idx] *= r;
    }
    x[N-1] = fma(-lo[N-1], x[N-2], x[N-1])
           / fma(-lo[N-1], up[N-2], diag[N-1]);

    // back substitution
    for idx in 0..N-2 by -1 {
      x[idx] = fma(-up[idx], x[idx+1], x[idx]);
    }
  }

}