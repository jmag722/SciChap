module KdTreeTest {
  use UnitTest;
  import SciChap.Spatial.KdTree;

  proc init_simple(test: borrowed Test) throws {
    var x: [{0..1, 0..1}] real = [1.0, 2.0; 3.0, 4.0];
    var tree : KdTree.kdTree = new KdTree.kdTree(x);
    test.assertEqual(tree.data, [1.0, 2.0; 3.0, 4.0]);
  }

  proc main() throws {
    UnitTest.main();
  }
}