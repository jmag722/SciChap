module StatisticsTest {
  use UnitTest;

  import SciChap.Statistics;

  proc median_1d(test: borrowed Test) throws {
    var arr = [5.0, 3.0, -1.0, 0.5];
    test.assertEqual(Statistics.median(arr), 1.75);
    test.assertEqual(Statistics.median([11, 3, 5]), 5.0);
  }
  proc median_2dslice(test: borrowed Test) throws {
    var arr = [
      11, 0, 5;
      12, -3, -19;
      44, 0, 0;
    ];
    test.assertEqual(Statistics.median(arr[.., 0]), 12.0);
    test.assertEqual(Statistics.median(arr[.., 1]), 0.0);
    test.assertEqual(Statistics.median(arr[.., 2]), 0.0);
    test.assertEqual(Statistics.median(arr[0, ..]), 5.0);
    test.assertEqual(Statistics.median(arr[1, ..]), -3.0);
    test.assertEqual(Statistics.median(arr[2, ..]), 0.0);
  }

  proc main() throws {
    UnitTest.main();
  }

}