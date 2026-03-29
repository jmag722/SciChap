/*
  **NOTE**

  This module is not intended to be imported directly. Instead, use

  .. code-block:: chapel

       import Spatial.KdTree
*/
module KdTreeMod {
  import Heap.heap;
  import Map.map;
  import Sort.{keyComparator, sort};

  import SciChap.Array;

  @chpldoc.nodoc
  record leafBucket {
    var dom: domain(1) = {1..0};
    // point indices in the bucket
    var pointIdxs: [dom] int;

    proc init() {}

    proc init(in pointIdxs: [?D] int) where D.rank == 1 {
      this.dom = {0..#D.size};
      this.pointIdxs = pointIdxs;
    }
  }

  @chpldoc.nodoc
  record nearestPtsQueue {

    // index of the point indices in each tuple on the heap
    param pointIdxIdx: int = 0;
    // index of the distances in each tuple on the heap
    param distIdx: int = 1;
    type heapType = (int, real);

    // heaps need a way to sort their tuples
    record tupleComparator: keyComparator {
      proc key(elt) do return elt[nearestPtsQueue.distIdx];
    }

    const maxSize: int;
    // a 'queue' where tuples of point index and distance to the queried
    // point are stored.
    var distQueue: heap(heapType, parSafe=true, tupleComparator);

    proc init(in maxSize) {
      this.maxSize = maxSize;
      this.distQueue = new heap(heapType, parSafe=true, new tupleComparator());
      init this;
    }
    proc size: int do return distQueue.size;
    proc isEmpty(): bool do return distQueue.isEmpty();
    proc isFull(): bool do return this.size >= maxSize;
    proc last: real do return distQueue.top()[distIdx];
    proc ref push(in newItem: heapType): void {
      distQueue.push(newItem);
      if this.size > maxSize then distQueue.pop();
    }
    proc toArray(): ([0..#size] int, [0..#size] real) {
      var heapArr: [0..#size] heapType = distQueue.toArray(); // unsorted
      sort(heapArr, comparator=new tupleComparator());
      var indices: [0..#size] int = [tup in heapArr] tup[pointIdxIdx];
      var distances: [0..#size] real = [tup in heapArr] tup[distIdx];
      return (indices, distances);
    }
  }

  /*
    A KdTree implementation that builds the tree in parallel. It currently uses
    the midpoint of max spread to determine the partitions.
  */
  class KdTree {
    @chpldoc.nodoc
    param ptsAxis: int = 0;

    @chpldoc.nodoc
    param dimAxis: int = 1;
    @chpldoc.nodoc
    const ptsDom: domain(2) = {1..0, 1..0};
    @chpldoc.nodoc
    const points: [ptsDom] real;

    @chpldoc.nodoc
    param emptyNodeVal: real = nan;
    @chpldoc.nodoc
    var nodesDom: domain(1) = {1..0};
    // the planar value at which the current node is split
    @chpldoc.nodoc
    var nodes: [nodesDom] real;

    @chpldoc.nodoc
    param emptyAxisVal: int = -1;
    // the axis of the planar value used in the split
    @chpldoc.nodoc
    var axes: [nodesDom] int;

    @chpldoc.nodoc
    const leafSize: int;
    @chpldoc.nodoc
    type leafBucketMap = map(int, leafBucket);
    @chpldoc.nodoc
    var leaves: leafBucketMap;

    /*
      Initialize the KdTree

      :arg points: n-dimensional set of points The size of the first
                    dimension is the number of points, and the size of the
                    second dimension is the number of dimensions

      :arg leafSize: minimum number of points needed before simple distance
                     method is used, defaults to unity

    */
    proc init(const in points: [?D] real, in leafSize: int=1): void
              where D.rank == 2 {
      this.ptsDom = {0..#D.shape[ptsAxis], 0..#D.shape[dimAxis]};
      this.points = points;


      // can't access proc npoints before init this
      var npoints: int = points.shape[ptsAxis];
      nodesDom = {0..#10*npoints};
      nodes = emptyNodeVal;
      axes = emptyAxisVal;

      this.leafSize = leafSize;
      leaves = new leafBucketMap();

      init this;

      var pointIdxsRange: range = ptsDom.dim(ptsAxis);
      var allPointIdxs: [pointIdxsRange] int = pointIdxsRange;
      constructTree(allPointIdxs);
    }

    @chpldoc.nodoc
    inline proc npoints: int {
      return points.shape[ptsAxis];
    }
    @chpldoc.nodoc
    inline proc ndim: int {
      return points.shape[dimAxis];
    }


    /*
      Construct the KdTree using a breadth-first approach, allowing for
      a parallel build process.

      Note this function does not need to be called by the user.

      :arg allPointIdxs: indices of all the points input by the user
    */
    @chpldoc.nodoc
    proc constructTree(const ref allPointIdxs: [] int): void {
      // the level index of the binary tree
      var level: int = 0;
      // all points start out as leaves in the initial leaf node, then are
      // subdivided out
      leaves.add(level, new leafBucket(allPointIdxs));
      while anyBucketAboveMinSize(leaves) {
        forall nodeIdx in KdTree.nodeIdxs(level) {
          const pointIdxs = leaves.get(nodeIdx, new leafBucket()).pointIdxs;

          // skip leaf nodes and nonexistent nodes
          if pointIdxs.size <= leafSize then continue;

          const (splitVal, splitAxis): (real, int) =
            splitMidpointMaxSpread(pointIdxs);

          // TODO: when chapel supports passing promoted expressions into
          // procs with array arguments, merge leftMask and leftPtIdxs logic
          const leftMask = points[pointIdxs, splitAxis] <= splitVal;
          const rightMask = points[pointIdxs, splitAxis] > splitVal;
          const nLeftPts: int = leftMask.count(true);
          const nRightPts: int = rightMask.count(true);

          // TODO: if leftMask.count(true) < 1, slide midpoint, adjust mask
          const leftPtIdxs: [0..#nLeftPts] int = pointIdxs[
                                                Array.trueIdxs(leftMask)];
          const rightPtIdxs: [0..#nRightPts] int = pointIdxs[
                                                  Array.trueIdxs(rightMask)];
          leaves.add(KdTree.childIdxLeft(nodeIdx), new leafBucket(leftPtIdxs));
          leaves.add(KdTree.childIdxRight(
            nodeIdx), new leafBucket(rightPtIdxs)
          );
          nodes[nodeIdx] = splitVal;
          axes[nodeIdx] = splitAxis;

          leaves.remove(nodeIdx);
        }
        level += 1;
      }
    }

    /*
     Find the indices of the points closest to the query point.

     :arg queryPoint: point being queried

     :arg nnearest: number of nearest data points to obtain, defaults to `1`.
                    If there are fewer points in the tree than specified by
                    `nnearest`, the query only output the max size of the tree.

     :returns: indices and distances of the N nearest data points, sorted
               from closest to farthest

     :throws IllegalArgumentError: queryPoint domain does not match KdTree
                                   domain
     */
    proc query(const queryPoint: [?queryD] real, in nnearest:int=1):
               ([0..#nnearest] int, [0..#nnearest] real) throws
               where queryD.rank == 1 {
      if queryD.size != ptsDom.shape[dimAxis] {
        throw new owned IllegalArgumentError(
          "queryPoint domain does not match KdTree points dimensionality");
      }
      nnearest = min(nnearest, this.npoints);
      var search: nearestPtsQueue = new nearestPtsQueue(nnearest);
      queryRecurse(queryPoint, nodeIdx=0, search);
      var (indices, distances) = search.toArray();
      [idx in distances.domain] distances[idx] **= 0.5;
      return (indices, distances);
    }

    @chpldoc.nodoc
    proc queryRecurse(const queryPoint: [] real, const nodeIdx: int,
                      ref search: nearestPtsQueue): void {
      if isEmptyNode(nodeIdx) then return;
      if isLeafNode(nodeIdx) {
        const leafIdxs = leaves.get(nodeIdx, new leafBucket()).pointIdxs;
        const nleaves = leafIdxs.size;
        if nleaves == 0 then halt(); // should never get here

        const dimRng: range = points.dim(dimAxis);
        var leafPoints: [{0..#nleaves, dimRng}] real;
        forall i in leafIdxs.domain {
          leafPoints[i, dimRng] = points[leafIdxs[i], dimRng];
        }

        const distsSq = [i in leafPoints.dim(ptsAxis)]
                         + reduce (leafPoints[i, ..] - queryPoint)**2;
        for (leafIdx, distSqr) in zip(leafIdxs, distsSq) {
          if !search.isFull() || distSqr < search.last {
            search.push((leafIdx, distSqr));
          }
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
        if !search.isFull() || search.last >= dist2planeSq {
          queryRecurse(queryPoint, KdTree.childIdxRight(nodeIdx), search);
        }
      } else {
        queryRecurse(queryPoint, KdTree.childIdxRight(nodeIdx), search);
        if !search.isFull() || search.last >= dist2planeSq {
          queryRecurse(queryPoint, KdTree.childIdxLeft(nodeIdx), search);
        }
      }
    }


    /*
     Find the indices of the points within a spherical radius of the query

     :arg queryPoint: point being queried

     :arg radius: Radius of spherical ball around queryPoint to search

     :returns: indices and distances of the points within the ball relative to
               the query, sorted from closest to farthest

     :rtype: ([] int, [] real)

     :throws IllegalArgumentError: queryPoint domain does not match KdTree
                                   domain
     */
    proc queryBallPoint(const queryPoint: [?queryD] real, in radius:real)
                        throws where queryD.rank == 1 {
      if queryD.size != ptsDom.shape[dimAxis] {
        throw new owned IllegalArgumentError(
          "queryPoint domain does not match KdTree points dimensionality");
      }
      // use nearestPtsQueue, but never pop any off because no limit for qty
      var search: nearestPtsQueue = new nearestPtsQueue(this.npoints);
      queryBallPointRecurse(queryPoint, nodeIdx=0, search=search,
                            radiusSqr=radius**2);
      var (indices, distances) = search.toArray();
      [idx in distances.domain] distances[idx] **= 0.5;
      return (indices, distances);
    }

    @chpldoc.nodoc
    proc queryBallPointRecurse(const queryPoint: [] real, const nodeIdx: int,
                               ref search: nearestPtsQueue,
                               in radiusSqr:real): void {
      if isEmptyNode(nodeIdx) then return;
      if isLeafNode(nodeIdx) {
        const leafIdxs = leaves.get(nodeIdx, new leafBucket()).pointIdxs;
        const nleaves = leafIdxs.size;
        if nleaves == 0 then halt(); // should never get here

        const dimRng: range = points.dim(dimAxis);
        var leafPoints: [{0..#nleaves, dimRng}] real;
        forall i in leafIdxs.domain {
          leafPoints[i, dimRng] = points[leafIdxs[i], dimRng];
        }

        const distsSq = [i in leafPoints.dim(ptsAxis)]
                         + reduce (leafPoints[i, ..] - queryPoint)**2;
        for (leafIdx, distSqr) in zip(leafIdxs, distsSq) {
          if distSqr <= radiusSqr {
            search.push((leafIdx, distSqr));
          }
        }
        return;
      }

      const currentAxis = axes[nodeIdx];
      const currentSplit = nodes[nodeIdx];
      // recurse, and if alternate branch plane is closer than the ball radius,
      // then search that branch too
      const dist2planeSq: real = (currentSplit - queryPoint[currentAxis])**2;
      if queryPoint[currentAxis] <= currentSplit {
        queryBallPointRecurse(queryPoint, KdTree.childIdxLeft(nodeIdx),
                              search, radiusSqr);
        if !search.isFull() || radiusSqr >= dist2planeSq {
          queryBallPointRecurse(queryPoint, KdTree.childIdxRight(nodeIdx),
                                search, radiusSqr);
        }
      } else {
        queryBallPointRecurse(queryPoint, KdTree.childIdxRight(nodeIdx),
                              search, radiusSqr);
        if !search.isFull() || radiusSqr >= dist2planeSq {
          queryBallPointRecurse(queryPoint, KdTree.childIdxLeft(nodeIdx),
                                search, radiusSqr);
        }
      }
    }

    @chpldoc.nodoc
    inline proc isLeafNode(in nodeIdx: int): bool {
      return leaves.contains(nodeIdx);
    }
    @chpldoc.nodoc
    inline proc isEmptyNode(in nodeIdx: int): bool {
      return axes[nodeIdx] == emptyAxisVal && !isLeafNode(nodeIdx);
    }

    /*
      Find the appropriate splitting value and axis for a subdomain of the
      data points.

      :arg unvisited: unvisited point indices yet to be partitioned

      :returns: splitting value (separator) of the hyperplane, and the index
                of the hyperplane axis where the partition occurs
     */
    @chpldoc.nodoc
    proc splitMidpointMaxSpread(const ref unvisited: [] int): (real, int) {
      const (minVals, maxVals) = findHyperRectangleDims(unvisited);
      var (maxSpread, maxSpreadAxis) = maxloc reduce zip(maxVals - minVals,
                                                         minVals.domain);
      var midpoint = 0.5 * (maxVals[maxSpreadAxis] + minVals[maxSpreadAxis]);
      return (midpoint, maxSpreadAxis);
    }

    @chpldoc.nodoc
    proc findHyperRectangleDims(const unvisited: [] int): ([0..#ndim] real,
                                                           [0..#ndim] real) {
      const dimRng: range = 0..#ndim;
      var minVals: [dimRng] real;
      var maxVals: [dimRng] real;
      forall axis in dimRng {
        minVals[axis] = min reduce points[unvisited, axis];
        maxVals[axis] = max reduce points[unvisited, axis];
      }
      return (minVals, maxVals);
    }

    @chpldoc.nodoc
    proc anyBucketAboveMinSize(const ref leaves: leafBucketMap): bool {
      for leafBucket in leaves.values() {
        if leafBucket.pointIdxs.size > leafSize then return true;
      }
      return false;
    }

    /*
      Get the left child node index from the parent node index.

      :arg parentNodeIdx: parent node index

      :returns: left child index
    */
    @chpldoc.nodoc
    inline proc type childIdxLeft(const in parentNodeIdx: int): int {
      return 2 * parentNodeIdx + 1;
    }
    /*
      Get the right child node index from the parent node index.

      :arg parentNodeIdx: parent node index

      :returns: right child index
    */
    @chpldoc.nodoc
    inline proc type childIdxRight(const in parentNodeIdx: int): int {
      return 2 * parentNodeIdx + 2;
    }

    /*
      Get parent node index from a given child node index

      :arg childNodeIdx: child node index

      :returns: parent node index
     */
    @chpldoc.nodoc
    inline proc type parentIdx(const in childNodeIdx: int): int {
      return floor(0.5 * childNodeIdx - 0.5): int;
    }

    /*
      Get the number of node indices for a given tree level

      :arg level: tree level (0-based)

      :returns: number of nodes at the given tree level
     */
    @chpldoc.nodoc
    inline proc type nIdxs(const in level: int): int {
      return 2**level;
    }

    /*
      Get the node indices for a given tree level

      :arg level: tree level (0-based)

      :returns: all node indices at the given tree level
     */
    @chpldoc.nodoc
    proc type nodeIdxs(const in level: int): range(int) {
      return nIdxs(level) - 1..nIdxs(level+1) - 2;
    }

  }

}