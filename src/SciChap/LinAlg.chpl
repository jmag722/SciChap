/*
  Linear algebra methods beyond those found in the standard LinearAlgebra
  library.

*/
module LinAlg {

  /*
    Serial tri-diagonal matrix algorithm (TDMA), or the Thomas algorithm.
    Matrix should be diagonally dominant (or close enough).

    The solution is computed in-place at the variable ``rhs``.

    Reference:
      Laszlo, E., Giles, M., & Appleyard, J. (2016).
      Manycore algorithms for batch scalar and block tridiagonal solvers.
      `ACM Transactions on Mathematical Software (TOMS)`, `42`(4), 1-36.

    :arg lower: lower diagonal, length ``N-1``

    :arg diagonal: diagonal, length ``N``

    :arg upper: upper diagonal, length ``N-1``

    :arg rhs: right hand side of system, length ``N``

    :returns: void (solution variable is ``rhs``)
  */
  proc tdma(const lower: [?Dm1] real, const diagonal: [?D] real,
            ref upper: [Dm1] real, ref rhs: [D] real): void
            where Dm1.rank == 1 && D.rank == 1 {
    const N = rhs.size;
    // index lo starting at 1 to avoid "-1" everywhere below
    const ref lo = lower.reindex(1..N-1);
    const ref diag = diagonal.reindex(0..N-1);
    ref up = upper.reindex(0..N-2);
    ref rh = rhs.reindex(0..N-1);

    // forward substitution
    rh[0] /= diag[0];
    up[0] /= diag[0];
    for idx in 1..N-2 {
      rh[idx] = (rh[idx] - lo[idx] * rh[idx-1])
             / (diag[idx] - lo[idx] * up[idx-1]);
      up[idx] /= diag[idx] - lo[idx] * up[idx-1];
    }
    rh[N-1] = (rh[N-1] - lo[N-1] * rh[N-2])
            / (diag[N-1] - lo[N-1] * up[N-2]);

    // back substitution
    for idx in 0..N-2 by -1 {
      rh[idx] -= up[idx] * rh[idx+1];
    }
  }

}