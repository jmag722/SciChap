module KdTreeMod {
  use Heap;
  use Map;
  use Math;
  use Sort;

  import SciChap.Array;

  record leafBucket {
    var dom: domain(1) = {1..0};
    var arr: [dom] int;

    proc init() {}

    proc init(in arr: [?D] int) where D.rank == 1 {
      this.dom = {0..#D.size};
      this.arr = arr;
    }
  }

  record kdTreeHeap {

    // index of the point index in each tuple on the heap
    param heapIdxIdx: int = 0;
    // index of the distance in each tuple on the heap
    param heapDistIdx: int = 1;
    type heapType = (int, real);
    record tupleComparator: keyComparator {
      proc key(elt) do return elt[kdTreeHeap.heapDistIdx];
    }

    const maxSize: int;
    var distQueue: heap(heapType, parSafe=true, tupleComparator);

    proc init(in maxSize) {
      this.maxSize = maxSize;
      this.distQueue = new heap(heapType, parSafe=true, new tupleComparator());
      init this;
    }
    proc size: int do return distQueue.size;
    proc isEmpty(): bool do return distQueue.isEmpty();
    proc isFull(): bool do return this.size >= maxSize;
    proc last: real do return distQueue.top()[heapDistIdx];
    proc ref push(in newItem: heapType): void {
      distQueue.push(newItem);
      if this.size > maxSize then distQueue.pop();
    }
    proc toArray(): ([0..#size] int, [0..#size] real) {
      var heapArr: [0..#size] heapType = distQueue.toArray(); // unsorted
      sort(heapArr, comparator=new tupleComparator());
      var indices: [0..#size] int = [tup in heapArr] tup[heapIdxIdx];
      var distances: [0..#size] real = [tup in heapArr] tup[heapDistIdx];
      return (indices, distances);
    }
  }

  class KdTree {
    param ptsAxis: int = 0;
    param dimAxis: int = 1;
    param emptyNodeVal: real = nan;
    param emptyAxisVal: int = -1;
    type leafBucketMap = map(int, leafBucket);
    const dataDom: domain(2) = {1..0, 1..0};
    const data: [dataDom] real;
    const leafSize: int;


    var nodesDom: domain(1) = {1..0};
    var nodes: [nodesDom] real;
    var axes: [nodesDom] int;
    var leaves: leafBucketMap;

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
      leaves = new leafBucketMap();
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

    proc constructTree(const ref pointIndices: [] int): void {
      var level: int = 0;
      leaves.add(level, new leafBucket(pointIndices));
      while anyBucketAboveMinSize(leaves) {
        forall levelIdx in KdTree.levelIdxs(level) {
          const pointIdxs = leaves.get(levelIdx, new leafBucket()).arr;

          // skip leaf nodes and nonexistent nodes
          if pointIdxs.size <= leafSize then continue;

          const (splitVal, splitAxis): (real, int) = findSplit(pointIdxs);

          // TODO: when chapel supports passing promoted expressions into
          // procs with array arguments, merge leftMask and leftArr logic
          const leftMask = data[pointIdxs, splitAxis] <= splitVal;
          const rightMask = data[pointIdxs, splitAxis] > splitVal;
          const nLeftNodes: int = leftMask.count(true);
          const nRightNodes: int = rightMask.count(true);

          // TODO: if leftMask.count(true) < 1, slide midpoint, adjust mask
          const leftArr: [0..#nLeftNodes] int = pointIdxs[
                                                Array.trueIdxs(leftMask)];
          const rightArr: [0..#nRightNodes] int = pointIdxs[
                                                  Array.trueIdxs(rightMask)];
          leaves.add(KdTree.childIdxLeft(levelIdx), new leafBucket(leftArr));
          leaves.add(KdTree.childIdxRight(levelIdx), new leafBucket(rightArr));
          nodes[levelIdx] = splitVal;
          axes[levelIdx] = splitAxis;

          leaves.remove(levelIdx);
        }
        level += 1;
      }
    }

    /*
     Find the indices of the points closest to the query point.

     :arg queryPoint: point being queried
     :type queryPoint: [] real

     :arg nnearest: number of nearest data points to obtain, defaults to `1`.
     If there are fewer points in the tree than specified by `nnearest`, the
     query only output the max size of the tree.
     :type nnearest: int

     :returns: indices and distances of the N nearest data points, sorted
     from closest to farthest
     :rtype: ([] int, [] real)

     :throws IllegalArgumentError: queryPoint domain does not match KdTree
      domain
     */
    proc query(const queryPoint: [?queryD] real, in nnearest:int=1):
               ([0..#nnearest] int, [0..#nnearest] real) throws
               where queryD.rank == 1 {
      if queryD.size != dataDom.shape[dimAxis] {
        throw new owned IllegalArgumentError(
          "queryPoint domain does not match KdTree data dimensionality");
      }
      nnearest = min(nnearest, this.npoints);
      var search: kdTreeHeap = new kdTreeHeap(nnearest);
      queryRecurse(queryPoint, nodeIdx=0, search);
      var (indices, distances) = search.toArray();
      [idx in distances.domain] distances[idx] **= 0.5;
      return (indices, distances);
    }

    proc queryRecurse(const queryPoint: [] real, const nodeIdx: int,
                      ref search: kdTreeHeap): void {
      if isEmptyNode(nodeIdx) then return;
      if isLeafNode(nodeIdx) {
        const leafIdxs = leaves.get(nodeIdx, new leafBucket()).arr;
        const nleaves = leafIdxs.size;
        if nleaves == 0 then halt(); // should never get here

        const dimRng: range = data.dim(dimAxis);
        var leafPoints: [{0..#nleaves, dimRng}] real;
        forall i in leafIdxs.domain {
          leafPoints[i, dimRng] = data[leafIdxs[i], dimRng];
        }

        const distSq = [i in leafPoints.dim(ptsAxis)]
                        + reduce (leafPoints[i, ..] - queryPoint)**2;
        const (newClosestIdx,
               newClosestDistSq) = minloc reduce zip(leafIdxs, distSq);
        if !search.isFull() || newClosestDistSq < search.last {
          search.push((newClosestIdx, newClosestDistSq));
        }
        return;
      }

      const currentAxis = axes[nodeIdx];
      const currentSplit = nodes[nodeIdx];
      // recurse, and if alternate branch plane is closer than the farthest
      // accumulated point we have, search that branch too
      const dist2planeSq: real = (currentSplit - queryPoint[currentAxis])**2;
      if queryPoint[currentAxis] <= currentSplit {
        queryRecurse(queryPoint, KdTree.childIdxLeft(nodeIdx), search);
        if search.last >= dist2planeSq {
          queryRecurse(queryPoint, KdTree.childIdxRight(nodeIdx), search);
        }
      } else {
        queryRecurse(queryPoint, KdTree.childIdxRight(nodeIdx), search);
        if search.last >= dist2planeSq {
          queryRecurse(queryPoint, KdTree.childIdxLeft(nodeIdx), search);
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

    proc anyBucketAboveMinSize(const ref leaves: leafBucketMap): bool {
      for leafBucket in leaves.values() {
        if leafBucket.arr.size > leafSize then return true;
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