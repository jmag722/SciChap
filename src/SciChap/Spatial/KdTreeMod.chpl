module KdTreeMod {
  use Math;

  class KdTree {
    const dataDom: domain(2) = {1..0, 1..0};
    const data: [dataDom] real;
    const leafSize: int;

    const npoints: int;
    const ndim: int;

    var nodesDom: domain(1) = {1..0};
    var nodes: [nodesDom] real;
    var axes: [nodesDom] int;

    proc init(const in data: [?D] real, in leafSize: int=1): void
              where D.rank == 2 {
      this.dataDom = {0..#D.shape[0], 0..#D.shape[1]};
      this.data = data;
      this.leafSize = leafSize;

      this.npoints = data.shape[0];
      this.ndim = data.shape[1];

      nodesDom = {0..#this.npoints};
      nodes = nan;
      axes = -1;
      init this;

      constructTree(nodesDom.first, this.nodes, this.data, this.axes);
    }
    proc constructTree(in level: int, const ref _mid: [?D] real,
                       const ref _dat: [] real,
                       const ref _axes: [D] int): void {
      // for each axis in data, compute spread (max - min)
      // choose axis with highest spread
      // cut axis w/highest spread in half orthogonal to that direction (midpt)
      // find all points to left of plane
      // find all points to right of plane
      // if left points empty, get min of right points as new node value
      // if right points empty, get min of left points as new node value
    }

    proc _findPivot(subDom: domain(1)) {
      var dimRng: range = 0..#ndim;
      var maxVals: [dimRng] real;
      var minVals: [dimRng] real;
      forall axis in dimRng {
        maxVals[axis] = max reduce data[subDom, axis];
        minVals[axis] = min reduce data[subDom, axis];
      }
      var (maxSpread, maxSpreadAxis) = maxloc reduce zip(maxVals - minVals,
                                                         dimRng);
      return (0.5 * maxSpread, maxSpreadAxis);
    }

    proc query(const ref pt: [?D] real): int where D.rank == 1 {

    }

  }

}