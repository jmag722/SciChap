module KdTreeMod {
  use Map;
  use Math;
  import SciChap.Array;

  record bucket {
    var dom: domain(1) = {1..0};
    var arr: [dom] int;

    proc init(in arr: [?D] int) {
      this.dom = D;
      this.arr = arr;
    }
  }

  class KdTree {
    param ptsAxis: int = 0;
    param dimAxis: int = 1;
    param emptyNodeVal: real = nan;
    param emptyAxisVal: int = -1;
    const dataDom: domain(2) = {1..0, 1..0};
    const data: [dataDom] real;
    const leafSize: int;

    var nodesDom: domain(1) = {1..0};
    var nodes: [nodesDom] real;
    var axes: [nodesDom] int;
    var leaves: map(int, bucket);

    proc init(const in data: [?D] real, in leafSize: int=1): void
              where D.rank == 2 {
      this.dataDom = {0..#D.shape[ptsAxis], 0..#D.shape[dimAxis]};
      this.data = data;
      this.leafSize = leafSize;

      // can't access proc npoints before init this
      var npoints: int = data.shape[ptsAxis];
      nodesDom = {0..#4*npoints};
      nodes = emptyNodeVal;
      axes = emptyAxisVal;
      leaves = new map(int, bucket);
      init this;

      var pointIndicesRange: range = dataDom.dim(ptsAxis);
      var pointIndices: [pointIndicesRange] int = pointIndicesRange;
      constructTree(pointIndices);
    }

    inline proc npoints: int {
      return data.shape[ptsAxis];
    }
    inline proc ndim: int {
      return data.shape[dimAxis];
    }

    proc constructTree(ref pointIndices: [] int): void {
      var level: int = 0;
      leaves.add(level, new bucket(pointIndices));
      while anyBucketAboveMinSize(leaves) {
        forall levelIdx in KdTree.levelIdxs(level) {
          var defaultBucket = new bucket(Array.empty(int));
          var pointIdxs = leaves.get(levelIdx, defaultBucket).arr;

          // skip leaf nodes and nonexistent nodes
          if pointIdxs.size <= leafSize then continue;

          var (splitVal, splitAxis): (real, int) = findSplit(pointIdxs);
          var leftMask = data[pointIdxs, splitAxis] <= splitVal;
          var rightMask = data[pointIdxs, splitAxis] > splitVal;
          var nLeftNodes: int = leftMask.count(true);
          var nRightNodes: int = rightMask.count(true);
          // TODO: if leftMask.count(true) < 1, slide midpoint, adjust mask
          var leftArr: [0..#nLeftNodes] int = Array.findAll(leftMask, true);
          var rightArr: [0..#nRightNodes] int = Array.findAll(rightMask, true);

          leaves.add(KdTree.childIdxLeft(levelIdx), new bucket(leftArr));
          leaves.add(KdTree.childIdxRight(levelIdx), new bucket(rightArr));
          nodes[levelIdx] = splitVal;
          axes[levelIdx] = splitAxis;

          leaves.remove(levelIdx);
        }
        level += 1;
      }
    }

    /*
      Find the appropriate splitting value and axis for a subdomain of the
      data points.

      :arg unvisited: continuous function on which the root will be computed
      :type unvisited: [] int

      :returns: splitting value (separator) of the hyperplane, and the index
      of the hyperplane axis where the split occurs
      :rtype: tuple(real, int)
     */
    proc findSplit(const ref unvisited: [] int): (real, int) {
      var dimRng: range = 0..#ndim;
      var maxVals: [dimRng] real;
      var minVals: [dimRng] real;
      forall axis in dimRng {
        // TODO: when chapel handles slicing sparse domains, `unvisited` can
        // be a sparse domain, and reduction is data[unvisited][.., axis]
        maxVals[axis] = max reduce data[unvisited, axis];
        minVals[axis] = min reduce data[unvisited, axis];
      }
      var (maxSpread, maxSpreadAxis) = maxloc reduce zip(maxVals - minVals,
                                                         dimRng);
      var midpoint = 0.5 * (maxVals[maxSpreadAxis] + minVals[maxSpreadAxis]);
      return (midpoint, maxSpreadAxis);
    }

    proc anyBucketAboveMinSize(const ref bucketMap: map(int, bucket)): bool {
      for bucket in bucketMap.values() {
        if bucket.arr.size > leafSize then return true;
      }
      return false;
    }

    /*
      Get the child indices from the parent index.

      :arg parentIdx: parent node index
      :type parentIdx: int

      :returns: left and right child indices
      :rtype: tuple(int, int)
     */
    inline proc type childIdxs(const in parentIdx: int): (int, int) {
      return (childIdxLeft(parentIdx), childIdxRight(parentIdx));
    }

    inline proc type childIdxLeft(const in parentIdx: int): int {
      return 2 * parentIdx + 1;
    }
    inline proc type childIdxRight(const in parentIdx: int): int {
      return 2 * parentIdx + 2;
    }

    /*
      Get parent index from a given child index

      :arg childIdx: child index
      :type childIdx: int

      :returns: parent index
      :rtype: int
     */
    inline proc type parentIdx(const in childIdx: int): int {
      return floor(0.5 * childIdx - 0.5): int;
    }

    /*
      Get the number of node indices for a given tree level

      :arg level: tree level (0-based)
      :type level: int

      :returns: number of nodes at the given tree level
      :rtype: int
     */
    inline proc type nIdxs(const in level: int): int {
      return 2**level;
    }

    /*
      Get the node indices for all indices of a given level

      :arg level: tree level (0-based)
      :type level: int

      :returns: all node indices at the given tree level
      :rtype: range(int)
     */
    proc type levelIdxs(const in level: int): range(int) {
      return nIdxs(level) - 1..nIdxs(level+1) - 2;
    }

  }

}