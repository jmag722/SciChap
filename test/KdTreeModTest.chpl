module KdTreeModTest {
  use UnitTest;
  import SciChap.Spatial;

  proc init_simple(test: borrowed Test) throws {
    var x: [{5..6, 11..12}] real = [1.0, 2.0; 3.0, 4.0];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x);
    test.assertEqual(tree.dataDom, {0..1, 0..1});
    test.assertEqual(tree.dataDom, tree.data.domain);
    test.assertEqual(tree.data, [1.0, 2.0; 3.0, 4.0]);
    test.assertEqual(x, [1.0, 2.0; 3.0, 4.0]);
  }

  proc _findPivot(test: borrowed Test) throws {
    var x: [{1..5, 1..3}] real = [1.0, 2.0, 13.0;
                                  3.0, 4.0, 14.0;
                                  5.0, 6.0, -3.0;
                                  7.0, -9.0, 0.0;
                                  2.0, 12.0, -6.0;];
    var tree : Spatial.KdTree = new owned Spatial.KdTree(x);
    var (pivot, axis) = tree._findPivot({x.domain.dim(0) - 1});
    test.assertEqual(pivot, 10.5);
    test.assertEqual(axis, 1);
  }

  proc main() throws {
    UnitTest.main();
  }
}