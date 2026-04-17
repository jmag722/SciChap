/*

Private kd-tree module

**NOTE**

This module is not intended to be imported directly. Instead, use

.. code-block:: chapel

      import Spatial.KdTree

..
  START_HERE
*/
module KdTreeMod {
  import Heap.heap;
  import Map.map;
  import Math;
  import Sort.{keyComparator, sort};

  import SciChap.Array;
  import SciChap.Statistics;

  @chpldoc.nodoc
  /* A "bucket" to hold the indices of the points for a given leaf node
  */
  record leafBucket {
    var dom: domain(1) = {1..0};
    /* point indices in the bucket */
    var pointIdxs: [dom] int;

    proc init() {}

    proc init(in pointIdxs: [?D] int) where D.rank == 1 {
      this.dom = {0..#D.size};
      this.pointIdxs = pointIdxs;
    }
  }

  /* This queue stores and sorts the list of point indices while the tree
  is being searched
  */
  @chpldoc.nodoc
  record nearestPtsQueue {
    /* the data type stored in the queue/heap */
    type heapType = (int, real);

    /* index of the point indices in each tuple on the heap */
    param pointIdxIdx: int = 0;
    /* index of the distances in each tuple on the heap */
    param distIdx: int = 1;

    // heaps need a way to sort their tuples
    record tupleComparator: keyComparator {
      proc key(elt) do return elt[nearestPtsQueue.distIdx];
    }

    /* the the max size this queue will store before removing points
      to add a new point */
    const maxSize: int;

    /* the actual 'queue' where tuples of point index and distance to the
      queried point are stored. */
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

    /*
      Add a new item to the queue

      :arg newItem: New heap element to add
    */
    proc ref push(in newItem: heapType): void {
      distQueue.push(newItem);
      if this.size > maxSize then distQueue.pop();
    }

    /*
      Convert the queue to a sorted array

      :returns: A tuple of sorted arrays of point indices and their
                corresponding distances to the queried point
     */
    proc toArray(): ([0..#size] int, [0..#size] real) {
      var heapArr: [0..#size] heapType = distQueue.toArray(); // unsorted
      sort(heapArr, comparator=new tupleComparator());
      var indices: [0..#size] int = [tup in heapArr] tup[pointIdxIdx];
      var distances: [0..#size] real = [tup in heapArr] tup[distIdx];
      return (indices, distances);
    }
  }

  enum Splitter {
    /* Apply a median split along the data dimension with the largest spread */
    MedianMaxSpread,
    /* Apply a sliding-midpoint split along the hyperrectangle dimension
      of greatest length */
    MidpointRectangle,
    /* Apply a midpoint split along the data dimension with the largest spread*/
    MidpointMaxSpread,
  }

  /*
    A KdTree implementation that builds the tree in parallel.
  */
  class KdTree {
    @chpldoc.nodoc
    /* axis along which point differ */
    param ptsAxis: int = 0;

    @chpldoc.nodoc
    /* axis along which the dimension differs */
    param dimAxis: int = 1;

    @chpldoc.nodoc
    const ptsDom: domain(2) = {1..0, 1..0};
    @chpldoc.nodoc
    /* input data points */
    const points: [ptsDom] real;

    @chpldoc.nodoc
    param emptyNodeVal: real = nan;
    @chpldoc.nodoc
    var nodesDom: domain(1) = {1..0};
    @chpldoc.nodoc
    /* the nodes of the Kdtree, which are planar values at which the current
      node of the binary tree are divided */
    var nodes: [nodesDom] real;

    @chpldoc.nodoc
    param emptyAxisVal: int = -1;
    /* the dimensional axis of the planar node value used in the split */
    @chpldoc.nodoc
    var axes: [nodesDom] int;

    @chpldoc.nodoc
    /* max number of leaves stored in a single bucket of a leaf node */
    const leafSize: int;

    @chpldoc.nodoc
    type leafBucketMap = map(int, leafBucket);
    /* a map for relating which nodes of the tree are leaf nodes */
    @chpldoc.nodoc
    var leaves: leafBucketMap;

    /* Method for partitioning the tree */
    @chpldoc.nodoc
    var splitter: Splitter;

    /*
      Initialize the KdTree

      :arg points: n-dimensional set of points The size of the first
                   dimension is the number of points, and the size of the
                   second dimension is the number of dimensions

      :arg leafSize: minimum number of points needed before simple distance
                     method is used

      :arg splitter: The desired method for partitioning the KdTree. Median-
                     based rules ensure log2(n) tree height. Sliding midpoint
                     rules prevent skinny rectangles and can lead to more
                     efficient searches.

      :arg memFactor: Multiplier for the tree size (increases memory). With
                      the midpoint and even median based splitting, when leaf
                      sizes are small, often the tree size can be larger than
                      the actual number of points (in degenerate cases).

    */
    proc init(const in points: [?D] real, in leafSize: int=10,
              in splitter: Splitter=Splitter.MedianMaxSpread,
              in memFactor: real=1.0): void
              where D.rank == 2 {
      this.ptsDom = {0..#D.shape[ptsAxis], 0..#D.shape[dimAxis]};
      this.points = points;

      // can't access `proc npoints` before init this
      var npoints: int = points.shape[ptsAxis];

      this.nodesDom = {0..#allocatedTreeSize(npoints, memFactor)};
      this.nodes = emptyNodeVal;
      this.axes = emptyAxisVal;

      this.leafSize = leafSize;
      this.leaves = new leafBucketMap();

      this.splitter = splitter;

      init this;

      var pointIdxsRange: range = ptsDom.dim(ptsAxis);
      var allPointIdxs: [pointIdxsRange] int = pointIdxsRange;
      constructTree(allPointIdxs);
    }

    /*
      Number of points in tree
    */
    inline proc npoints: int {
      return points.shape[ptsAxis];
    }

    /*
      Dimensionality of the tree
    */
    inline proc ndim: int {
      return points.shape[dimAxis];
    }

    /*
      Dimensions of the tree
    */
    inline proc dimRng: range {
      return points.dim(dimAxis);
    }

    @chpldoc.nodoc
    // Domain of the hyperrectangle min-max bounds
    inline proc rectDom: domain(1) do return {0..#2*ndim};


    /*
      Construct the KdTree using a breadth-first approach, allowing for
      a parallel build process.

      Note this function does not need to be called by the user.

      :arg allPointIdxs: indices of all the points input by the user
    */
    @chpldoc.nodoc
    proc constructTree(const ref allPointIdxs: [] int): void {
      // the level index of the binary tree, starting at the root
      var level: int = 0;

      // all points start out as leaves to the root node, then are split up
      leaves.add(level, new leafBucket(allPointIdxs));

      // initialize the hyperrectangle for the entire domain
      var rectMinMax: [rectDom] real = findHyperRectangleBounds(allPointIdxs);

      /*
        Get boolean masks for the points on the left and right side of the tree
      */
      proc getLRMasks(const _ptIdxs: [?D] int, const splitVal: real,
                      const splitAxis: int): ([0..#D.size] bool,
                                              [0..#D.size] bool) {
        var left: [0..#D.size] bool = points[_ptIdxs, splitAxis] <= splitVal;
        var right: [0..#D.size] bool = points[_ptIdxs, splitAxis] > splitVal;
        return (left, right);
      }

      /*
       Slide the current partition of either the left or right group of points
       is empty, such that there is at least 1 point in each
      */
      proc slidePartition(ref leftMask: [?D] bool, ref rightMask: [D] bool,
                          ref splitVal: real, const splitAxis: int,
                          const _ptIdxs: [D] int): void {
        if leftMask.count(true) == 0 {
          var (_splitVal, minIdx) = minloc reduce zip(
            points[_ptIdxs, splitAxis], _ptIdxs.domain
          );
          // set masks explicitly like this to avoid float comparison errors
          // that we could get by recreating a mask
          leftMask[minIdx] = true;
          rightMask[minIdx] = false;
          splitVal = _splitVal;
        } else if rightMask.count(true) == 0 {
          var (_splitVal, maxIdx) = maxloc reduce zip(
            points[_ptIdxs, splitAxis], _ptIdxs.domain
          );
          leftMask[maxIdx] = false;
          rightMask[maxIdx] = true;
          splitVal = _splitVal;
        }
      }

      while anyBucketAboveMinSize(leaves) {
        forall nodeIdx in KdTree.nodeIdxs(level) {
          const pointIdxs = leaves.get(nodeIdx, new leafBucket()).pointIdxs;

          // skip leaf nodes and nonexistent nodes
          if pointIdxs.size <= leafSize then continue;

          var splitVal:real;
          var splitAxis: int;
          var leftMask: [pointIdxs.domain] bool;
          var rightMask: [pointIdxs.domain] bool;
          select splitter {
            when Splitter.MedianMaxSpread {
              (splitVal, splitAxis) = splitMedianMaxSpread(pointIdxs);
              (leftMask, rightMask) = getLRMasks(pointIdxs, splitVal,
                splitAxis);
              slidePartition(leftMask, rightMask, splitVal, splitAxis,
                             pointIdxs);
            }
            when Splitter.MidpointRectangle {
              (splitVal, splitAxis) = splitMidpointRectangle(nodeIdx,
                rectMinMax, level);
              (leftMask, rightMask) = getLRMasks(pointIdxs, splitVal,
                splitAxis);
              slidePartition(leftMask, rightMask, splitVal, splitAxis,
                pointIdxs);
            }
            when Splitter.MidpointMaxSpread {
              (splitVal, splitAxis) = splitMidpointMaxSpread(pointIdxs);
              (leftMask, rightMask) = getLRMasks(pointIdxs, splitVal,
                splitAxis);
            }
            otherwise do halt();
          }

          const leftPtIdxs: [0..#leftMask.count(true)] int = pointIdxs[
                                                Array.trueIdxs(leftMask)];
          const rightPtIdxs: [0..#rightMask.count(true)] int = pointIdxs[
                                                  Array.trueIdxs(rightMask)];
          leaves.add(KdTree.childIdxLeft(nodeIdx), new leafBucket(leftPtIdxs));
          leaves.add(KdTree.childIdxRight(nodeIdx),
                     new leafBucket(rightPtIdxs));

          nodes[nodeIdx] = splitVal;
          axes[nodeIdx] = splitAxis;
          // current node is no longer a leaf node
          leaves.remove(nodeIdx);
        }
        level += 1;
      }
    }

    /*
     Find the indices of the points closest to the query point.

     :arg queryPoint: point being queried

     :arg nnearest: number of nearest data points to obtain.
                    If there are fewer points in the tree than specified by
                    `nnearest`, the query only outputs the max size of the tree.

     :returns: indices and distances of the N nearest data points, sorted
               from closest to farthest

     :throws IllegalArgumentError: queryPoint domain does not match KdTree
                                   domain
     */
    proc query(const queryPoint: [?queryD] real, in nnearest:int=1):
               ([0..#nnearest] int, [0..#nnearest] real) throws
               where queryD.rank == 1 {
      if queryD.size != ndim {
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
      if queryD.size != ndim {
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
                               const ref radiusSqr:real): void {
      if isEmptyNode(nodeIdx) then return;
      if isLeafNode(nodeIdx) {
        const leafIdxs = leaves.get(nodeIdx, new leafBucket()).pointIdxs;
        const nleaves = leafIdxs.size;
        if nleaves == 0 then halt(); // should never get here

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
        if radiusSqr >= dist2planeSq {
          queryBallPointRecurse(queryPoint, KdTree.childIdxRight(nodeIdx),
                                search, radiusSqr);
        }
      } else {
        queryBallPointRecurse(queryPoint, KdTree.childIdxRight(nodeIdx),
                              search, radiusSqr);
        if radiusSqr >= dist2planeSq {
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
      Split the unvisited data points at the midpoint, along the axis with
      the greatest spread. Probably the most efficient splitting algorithm.

      :arg unvisited: unvisited point indices yet to be partitioned

      :returns: splitting value (separator) of the hyperplane, and the index
                of the hyperplane axis where the partition occurs
     */
    @chpldoc.nodoc
    proc splitMidpointMaxSpread(const ref unvisited: [] int): (real, int) {
      // the spread of points becomes the hyperrectangle to split
      const bounds = findHyperRectangleBounds(unvisited);
      return splitMidpoint(bounds);
    }

    /*
     Split the unvisited data points at the midpoint, along the axis of the
     rectangle of greatest length.

     :arg nodeIdx: The node at which we wish to compute the hyperrectangle

     :rectMinMax: the initial hyperrectangle to split

     :level: The level of the tree where nodeIdx is. This could easily be
             derived from nodeIdx, but it's already computed and cheap
             to pass in
    */
    @chpldoc.nodoc
    proc splitMidpointRectangle(const nodeIdx:int,
                                in rectMinMax: [rectDom] real,
                                const level: int): (real, int) {
      const parentIdxs = allParentIdxs(nodeIdx, level);
      // iterate from root node down to current node, splitting hyperrectangle
      for idx in parentIdxs.domain[1..] by -1 {
        var splitVal: real = this.nodes[parentIdxs[idx]];
        var splitAxis: int = this.axes[parentIdxs[idx]];
        var nextIdx = idx - 1;
        // whether the child "nextIdx" node index is even or odd determines if
        // we went left or right, thus deciding which index of rectMinMax to
        // update
        rectMinMax[ndim * (parentIdxs[nextIdx] % 2) + splitAxis] = splitVal;
      }

      return splitMidpoint(rectMinMax);
    }

    /*
      Split the unvisited data points at the median, along the axis with
      the greatest spread.

      :arg unvisited: unvisited point indices yet to be partitioned

      :returns: splitting value (separator) of the hyperplane, and the index
                of the hyperplane axis where the partition occurs
     */
    @chpldoc.nodoc
    proc splitMedianMaxSpread(const ref unvisited: [] int): (real, int) {
      const bounds = findHyperRectangleBounds(unvisited);
      var (_, maxSpreadAxis) = boundsMaxSpread(bounds);
      var slice: [0..#unvisited.size] real = points[unvisited, maxSpreadAxis];
      return (Statistics.median(slice), maxSpreadAxis);
    }

    /*
      Split a hyperrectangle at the midpoint, along the side with greatest
      length.

      :arg minVals: min hyperrectangle values along each dimension

      :arg maxVals: max hyperrectangle values along each dimension

      :returns: splitting value (separator) of the hyperplane, and the index
                of the hyperplane axis where the partition occurs
     */
    @chpldoc.nodoc
    proc splitMidpoint(const bounds: [rectDom] real): (real, int) {
      var (_, maxSpreadAxis) = boundsMaxSpread(bounds);
      var midpoint = 0.5 * (bounds[maxSpreadAxis] + bounds[ndim+maxSpreadAxis]);
      return (midpoint, maxSpreadAxis);
    }

    @chpldoc.nodoc
    proc boundsMaxSpread(const bounds: [rectDom] real): (real, int) {
      return maxloc reduce zip(bounds[ndim..#ndim] - bounds[0..#ndim],
                               bounds.domain[0..#ndim]);
    }

    @chpldoc.nodoc
    proc findHyperRectangleBounds(const unvisited: [] int): [rectDom] real {
      var bounds: [rectDom] real;
      forall axis in dimRng {
        bounds[axis] = min reduce points[unvisited, axis];
        bounds[ndim + axis] = max reduce points[unvisited, axis];
      }
      return bounds;
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

  @chpldoc.nodoc
  proc allocatedTreeSize(in npoints: int, in memScale:real): int {
    return Math.ceil(memScale*npoints):int;
  }

  /*
    Compute the ancestry of a given node index - all node indices above it
    to the root node

    :arg nodeIdx: current node index

    :arg level: level at which the `nodeIdx` resides. could be computed without
                but we already know it

    :returns: starting from and including `nodeIdx`, all the parent nodes of
              `nodeIdx`
  */
  @chpldoc.nodoc
  proc allParentIdxs(const nodeIdx: int, const level: int): [0..#level+1] int {
    const nAncestors = level+1;
    var ancestors: [0..#nAncestors] int;
    ancestors[0] = nodeIdx;
    for idx in ancestors.domain[1..] {
      ancestors[idx] = KdTree.parentIdx(ancestors[idx-1]);
    }
    return ancestors;
  }

}