module KdTreeMod {
  use Map;
  use Math;
  import SciChap.Array;

  record bucket {
    var dom: domain(1) = {1..0};
    var arr: [dom] int;

    proc init() {}

    proc init(in arr: [?D] int) where D.rank == 1 {
      this.dom = {0..#D.size};
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
          const pointIdxs = leaves.get(levelIdx, new bucket()).arr;

          // skip leaf nodes and nonexistent nodes
          if pointIdxs.size <= leafSize then continue;

          const (splitVal, splitAxis): (real, int) = findSplit(pointIdxs);
          const leftMask = data[pointIdxs, splitAxis] <= splitVal;
          const rightMask = data[pointIdxs, splitAxis] > splitVal;
          const nLeftNodes: int = leftMask.count(true);
          const nRightNodes: int = rightMask.count(true);
          // TODO: if leftMask.count(true) < 1, slide midpoint, adjust mask
          const leftArr: [0..#nLeftNodes] int = Array.boolSlice(pointIdxs,
                                                                leftMask);
          const rightArr: [0..#nRightNodes] int = Array.boolSlice(pointIdxs,
                                                                  rightMask);
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
     Find the index of the data point closest to the query point.

     :arg queryPoint: point being queried
     :type queryPoint: [] real

     :returns: index of the closest data point
     :rtype: int

     :throws IllegalArgumentError: queryPoint domain does not match KdTree
      domain
     */
    proc query(const queryPoint: [?queryD] real): int throws
               where queryD.rank == 1 {
      if queryD.size != dataDom.shape[dimAxis] {
        throw new owned IllegalArgumentError(
          "queryPoint domain does not match KdTree data dimensionality");
      }
      var closestIdx: int = -1;
      var closestDistSq: real = 1e99;
      queryRecurse(queryPoint, nodeIdx=0, closestIdx, closestDistSq);
      return closestIdx;
    }

    proc queryRecurse(const queryPoint: [] real, const nodeIdx: int,
                      ref closestIdx: int, ref closestDistSq: real): void {
      if isEmptyNode(nodeIdx) then return;
      if isLeafNode(nodeIdx) {
        const bucketPtIdxs = leaves.get(nodeIdx, new bucket()).arr;
        const nBucketPts = bucketPtIdxs.size;
        if nBucketPts == 0 then halt(); // should never get here

        const dimRng: range = data.dim(dimAxis);
        var bucketPoints: [{0..#nBucketPts, dimRng}] real;
        forall i in bucketPtIdxs.domain {
          bucketPoints[i, dimRng] = data[bucketPtIdxs[i], dimRng];
        }

        const distSq = [i in bucketPoints.dim(ptsAxis)]
                        + reduce (bucketPoints[i, ..] - queryPoint)**2;
        const (newClosestDistSq,
               newClosestIdx) = minloc reduce zip(distSq, bucketPtIdxs);
        if newClosestDistSq < closestDistSq {
          closestDistSq = newClosestDistSq;
          closestIdx = newClosestIdx;
        }
        return;
      }

      const currentAxis = axes[nodeIdx];
      const currentSplit = nodes[nodeIdx];
      const dist2planeSq: real = (currentSplit - queryPoint[currentAxis])**2;
      if queryPoint[currentAxis] <= currentSplit {
        queryRecurse(queryPoint, KdTree.childIdxLeft(nodeIdx),
                     closestIdx, closestDistSq);
        if closestDistSq >= dist2planeSq {
          queryRecurse(queryPoint, KdTree.childIdxRight(nodeIdx),
                       closestIdx, closestDistSq);
        }
      }
      else {
        queryRecurse(queryPoint, KdTree.childIdxRight(nodeIdx),
                     closestIdx, closestDistSq);
        if closestDistSq >= dist2planeSq {
          queryRecurse(queryPoint, KdTree.childIdxLeft(nodeIdx),
                       closestIdx, closestDistSq);
        }
      }
    }

    inline proc isLeafNode(in nodeIdx: int): bool {
      return leaves.contains(nodeIdx);
    }
    inline proc isEmptyNode(in nodeIdx: int): bool {
      return axes[nodeIdx] == emptyAxisVal && !isLeafNode(nodeIdx);
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