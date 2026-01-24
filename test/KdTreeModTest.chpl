module KdTreeModTest {
  use Map;
  use Math;
  use Sort;
  use UnitTest;
  import SciChap.Spatial;

  /*
   Determine if 2 1D arrays are equal, with NaN meaning they are.
   */
  proc assertEqualsNanArray(const ref actual: [?D] real,
                            const ref expected: [D] real): bool
                            where D.rank == 1 {
    var equals: bool = true;
    for i in D {
      if isNan(expected[i]) {
        equals &&= isNan(actual[i]);
      }
      else {
        equals &&= actual[i] == expected[i];
      }
    }
    return equals;
  }

  proc init_simple(test: borrowed Test) throws {
    var x: [{5..6, 11..12}] real = [1.0, 2.0; 3.0, 4.0];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x);
    test.assertEqual(tree.dataDom, {0..1, 0..1});
    test.assertEqual(tree.dataDom, tree.data.domain);
    test.assertEqual(tree.data, [1.0, 2.0; 3.0, 4.0]);
    test.assertEqual(x, [1.0, 2.0; 3.0, 4.0]);

    test.assertEqual(tree.nodesDom, {0..#4*x.shape[tree.ptsAxis]});
    var en: real = tree.emptyNodeVal;
    test.assertTrue(assertEqualsNanArray(tree.nodes,
                    [2.0, en, en, en, en, en, en, en]));
    var ea: int = tree.emptyAxisVal;
    test.assertEqual(tree.axes, [0, ea, ea, ea, ea, ea, ea, ea]);
    var expectedLeaves = new map(int, Spatial.leafBucket);
    expectedLeaves.add(1, new Spatial.leafBucket([0]));
    expectedLeaves.add(2, new Spatial.leafBucket([1]));
    test.assertEqual(tree.leaves, expectedLeaves);
  }

  proc findSplit_entireDomain(test: borrowed Test) throws {
    var x: [{0..4, 0..2}] real = [1.0, 2.0, 13.0;
                                  3.0, 4.0, 14.0;
                                  5.0, 6.0, -3.0;
                                  7.0, -9.0, 0.0;
                                  2.0, 12.0, -6.0;];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x);
    var (value, axis) = tree.findSplit([0,1,2,3,4]);
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
      test.assertEqual(distances[0], expectedDist);
    }
  }

  proc query_allnearest2D(test: borrowed Test) throws {
    var x: [{0..4, 0..2}] real = [1.0, 2.0, 13.0;
                                  3.0, 4.0, 14.0;
                                  5.0, 6.0, -3.0;
                                  7.0, -9.0, 0.0;
                                  2.0, 12.0, -6.0;];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x);
    var queryPoint = [0.0, -1.0, 0.25];
    var (indices, distances) = tree.query(queryPoint, nnearest=5);
    test.assertEqual(indices, [2, 3, 0, 4, 1]);
    var expectDists = dist2query(x, queryPoint);
    sort(expectDists);
    test.assertEqual(distances, expectDists);
  }

  proc query_nnearestTooBig(test: borrowed Test) throws {
    var x: [{0..0, 0..2}] real = [1.0, 2.0, 13.0;];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x);
    var queryPoint = x[0, ..] + 0.1 * x[0, ..];
    var (indices, distances) = tree.query(queryPoint, 11);
    test.assertEqual(indices.size, 1);
    test.assertEqual(indices[0], 0);
    var expectedDist = dist2query(x[0, ..], queryPoint);
    test.assertEqual(distances[0], expectedDist);
  }

  proc findSplit_subdomain(test: borrowed Test) throws {
    var x: [{0..4, 0..2}] real = [1.0, 2.0, 13.0;
                                  3.0, 4.0, 14.0;
                                  5.0, 6.0, -3.0;
                                  7.0, -9.0, 0.0;
                                  2.0, 12.0, -6.0;];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x);
    var (value, axis) = tree.findSplit([0,2,3]);
    test.assertEqual(value, 5.0);
    test.assertEqual(axis, 2);
  }

  proc childIdxs(test: borrowed Test) throws {
    test.assertEqual(Spatial.KdTree.childIdxs(0), (1, 2));
    test.assertEqual(Spatial.KdTree.childIdxs(4), (9, 10));
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

  proc levelIdxs(test: borrowed Test) throws {
    test.assertEqual(Spatial.KdTree.levelIdxs(0), 0..0);
    test.assertEqual(Spatial.KdTree.levelIdxs(1), 1..2);
    test.assertEqual(Spatial.KdTree.levelIdxs(2), 3..6);
    test.assertEqual(Spatial.KdTree.levelIdxs(3), 7..14);
  }

  proc main() throws {
    UnitTest.main();
  }
}