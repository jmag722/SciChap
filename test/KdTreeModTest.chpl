module KdTreeModTest {
  use Map;
  use Math;
  use Sort;
  use UnitTest;
  import SciChap.Spatial;

  proc init_simple(test: borrowed Test) throws {
    var x: [{5..6, 11..12}] real = [1.0, 2.0; 3.0, 4.0];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x);
    test.assertEqual(tree.ptsDom, {0..1, 0..1});
    test.assertEqual(tree.ptsDom, tree.points.domain);
    test.assertEqual(tree.points, [1.0, 2.0; 3.0, 4.0]);
    test.assertEqual(x, [1.0, 2.0; 3.0, 4.0]);

    var expectedDom = {0..#10*x.shape[tree.ptsAxis]};
    test.assertEqual(tree.nodesDom, expectedDom);

    var en: real = tree.emptyNodeVal;
    var expectedNodeArr: [expectedDom] real = en;
    expectedNodeArr[0] = 2.0;
    test.assertClose(tree.nodes, expectedNodeArr, relTol=1e-15);

    var ea: int = tree.emptyAxisVal;
    var expectedAxisArr: [expectedDom] int = ea;
    expectedAxisArr[0] = 0;
    test.assertEqual(tree.axes, expectedAxisArr);

    var expectedLeaves = new map(int, Spatial.leafBucket);
    expectedLeaves.add(1, new Spatial.leafBucket([0]));
    expectedLeaves.add(2, new Spatial.leafBucket([1]));
    test.assertEqual(tree.leaves, expectedLeaves);
  }

  proc splitMidpointMaxSpread_entireDomain(test: borrowed Test) throws {
    var x: [{0..4, 0..2}] real = [1.0, 2.0, 13.0;
                                  3.0, 4.0, 14.0;
                                  5.0, 6.0, -3.0;
                                  7.0, -9.0, 0.0;
                                  2.0, 12.0, -6.0;];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x);
    var (value, axis) = tree.splitMidpointMaxSpread([0,1,2,3,4]);
    test.assertEqual(value, 1.5);
    test.assertEqual(axis, 1);
  }

  proc dist2query(const pts: [] real, const query: [] real): [] real
       where pts.rank == 2 && query.rank == 1 {
    var expectDists: [pts.dim(0)] real;
    forall queryIdx in expectDists.domain {
      expectDists[queryIdx] = dist2query(pts[queryIdx, ..], query);
    }
    return expectDists;
  }
  proc dist2query(const pts: [] real, const query: [] real): real
       where pts.rank == 1 && query.rank == 1 {
    return (+ reduce (pts - query)**2)**0.5;
  }

  proc query_1nearest2D(test: borrowed Test) throws {
    var x: [{0..4, 0..2}] real = [1.0, 2.0, 13.0;
                                  3.0, 4.0, 14.0;
                                  5.0, 6.0, -3.0;
                                  7.0, -9.0, 0.0;
                                  2.0, 12.0, -6.0;];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x);
    forall queryIdx in x.dim(0) {
      var queryPoint = x[queryIdx, ..] + 0.1 * x[queryIdx, ..];
      var (indices, distances) = tree.query(queryPoint);
      var expectedDist = dist2query(x[queryIdx, ..], queryPoint);
      test.assertEqual(indices[0], queryIdx);
      test.assertClose(distances[0], expectedDist, relTol=1e-15);
    }
  }

  proc query_allnearest2D(test: borrowed Test) throws {
    var x: [{0..4, 0..2}] real = [1.0, 2.0, 13.0;
                                  3.0, 4.0, 14.0;
                                  5.0, 6.0, -3.0;
                                  7.0, -9.0, 0.0;
                                  2.0, 12.0, -6.0;];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x, leafSize=3);
    var queryPoint = [0.0, -1.0, 0.25];
    var (indices, distances) = tree.query(queryPoint, nnearest=5);
    test.assertEqual(indices, [2, 3, 0, 4, 1]);
    var expectDists = dist2query(x, queryPoint);
    sort(expectDists);
    test.assertClose(distances, expectDists, relTol=1e-15);
  }

  proc query_nnearestTooBig(test: borrowed Test) throws {
    var x: [{0..0, 0..2}] real = [1.0, 2.0, 13.0;];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x);
    var queryPoint = x[0, ..] + 0.1 * x[0, ..];
    var (indices, distances) = tree.query(queryPoint, 11);
    test.assertEqual(indices.size, 1);
    test.assertEqual(indices[0], 0);
    var expectedDist = dist2query(x[0, ..], queryPoint);
    test.assertClose(distances[0], expectedDist, relTol=1e-15);
  }

  proc query2D(test: borrowed Test) throws {
    var x = [
      -1.1, 0.0;
      -1.2, 0.0;
       0.0, 0.1;
       0.0, 0.2;
       0.0, 0.0;
       1.0, 0.0;
       0.3, 0.4;
       5.0, 3.6;
    ];
    var tree: Spatial.KdTree = new owned Spatial.KdTree(x, leafSize=2);
    var queryPoint = [0.0, 0.0];
    var (indices, distances) = tree.query(queryPoint, nnearest=8);
    test.assertEqual(indices, [4, 2, 3, 6, 5, 0, 1, 7]);
    var expectDists = dist2query(x, queryPoint);
    sort(expectDists);
    test.assertClose(distances, expectDists, relTol=1e-15);
  }

  proc query1D(test: borrowed Test) throws {
    var x = [
      -1.1;
      -1.2;
       0.1;
       0.0;
       1.0;
       0.3;
       5.0;
    ];
    var tree: Spatial.KdTree = new owned Spatial.KdTree(x, leafSize=2);
    var queryPoint = [0.0];
    var (indices, distances) = tree.query(queryPoint, nnearest=7);
    test.assertEqual(indices, [3, 2, 5, 4, 0, 1, 6]);
    var expectDists = dist2query(x, queryPoint);
    sort(expectDists);
    test.assertClose(distances, expectDists, relTol=1e-15);
  }

  proc queryBallPoint(test: borrowed Test) throws {
    var x = [
      -1.1, 0.0;
      -1.2, 0.0;
       0.0, 1.3;
       0.0, 1.35;
       0.0, 0.1;
       0.0, 0.2;
       0.0, 0.0;
       1.0, 0.0; // right on the edge, should be included
       0.3, 0.4;
       5.0, 3.6;
    ];
    var tree: Spatial.KdTree = new owned Spatial.KdTree(x, leafSize=2);
    var queryPoint = [0.0, 0.0];
    var (indices, distances) = tree.queryBallPoint(queryPoint, radius=1);
    test.assertEqual(indices, [6, 4, 5, 8, 7]);
    var allDists = dist2query(x, queryPoint);
    sort(allDists);
    var expectedDists = allDists[0..#5];
    test.assertClose(distances, expectedDists, relTol=1e-15);
  }

  proc splitMidpointMaxSpread_subdomain(test: borrowed Test) throws {
    var x: [{0..4, 0..2}] real = [1.0, 2.0, 13.0;
                                  3.0, 4.0, 14.0;
                                  5.0, 6.0, -3.0;
                                  7.0, -9.0, 0.0;
                                  2.0, 12.0, -6.0;];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x);
    var (value, axis) = tree.splitMidpointMaxSpread([0,2,3]);
    test.assertEqual(value, 5.0);
    test.assertEqual(axis, 2);
  }

  proc childIdxLeft(test:borrowed Test) throws {
    test.assertEqual(Spatial.KdTree.childIdxLeft(0), 1);
    test.assertEqual(Spatial.KdTree.childIdxLeft(3), 7);
    test.assertEqual(Spatial.KdTree.childIdxLeft(4), 9);
  }
  proc childIdxRight(test:borrowed Test) throws {
    test.assertEqual(Spatial.KdTree.childIdxRight(0), 2);
    test.assertEqual(Spatial.KdTree.childIdxRight(3), 8);
    test.assertEqual(Spatial.KdTree.childIdxRight(4), 10);
  }
  proc parentIdx(test: borrowed Test) throws {
    test.assertEqual(Spatial.KdTree.parentIdx(1), 0);
    test.assertEqual(Spatial.KdTree.parentIdx(2), 0);
    test.assertEqual(Spatial.KdTree.parentIdx(3), 1);
    test.assertEqual(Spatial.KdTree.parentIdx(11), 5);
    test.assertEqual(Spatial.KdTree.parentIdx(14), 6);
  }

  proc nIdxs(test: borrowed Test) throws {
    test.assertEqual(Spatial.KdTree.nIdxs(0), 1);
    test.assertEqual(Spatial.KdTree.nIdxs(1), 2);
    test.assertEqual(Spatial.KdTree.nIdxs(2), 4);
    test.assertEqual(Spatial.KdTree.nIdxs(3), 8);
    test.assertEqual(Spatial.KdTree.nIdxs(10), 1024);
  }

  proc nodeIdxs(test: borrowed Test) throws {
    test.assertEqual(Spatial.KdTree.nodeIdxs(0), 0..0);
    test.assertEqual(Spatial.KdTree.nodeIdxs(1), 1..2);
    test.assertEqual(Spatial.KdTree.nodeIdxs(2), 3..6);
    test.assertEqual(Spatial.KdTree.nodeIdxs(3), 7..14);
  }

  proc main() throws {
    UnitTest.main();
  }
}