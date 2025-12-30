module ArrayTest {
  use UnitTest;
  import SciChap.Array;

  proc diff_uniform(test: borrowed Test) throws {
    test.assertEqual(Array.diff([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]),
                                [1.0, 1.0, 1.0, 1.0, 1.0]);
  }
  proc diff_nonuniform(test: borrowed Test) throws {
    test.assertEqual(Array.diff([-12.0, 3.0, 7.0, 5.0, -60.1]),
                                [15.0, 4.0, -2.0, -65.1]);
  }
  proc diff_int(test: borrowed Test) throws{
    test.assertEqual(Array.diff([1, 3, -9]), [2, -12]);
  }
  proc diff_complex(test: borrowed Test) throws{
    test.assertEqual(Array.diff([1+2i, 3+5i]), [2+3i]);
  }
  proc diff_len2(test: borrowed Test) throws{
    test.assertEqual(Array.diff([1.0, 3.0]), [2.0]);
  }
  proc diff_len1(test: borrowed Test) throws{
    var empty: [1..0] real;
    test.assertEqual(Array.diff([1.0]), empty);
  }

  proc isMonotonic_len1(test: borrowed Test) throws {
    test.assertTrue(Array.isMonotonic([1]));
    test.assertTrue(Array.isMonotonic([-5.0]));
  }
  proc isMonotonic_len2(test: borrowed Test) throws {
    test.assertTrue(Array.isMonotonic([15, 25]));
    test.assertTrue(Array.isMonotonic([-5.1, 3.6]));
  }
  proc isMonotonic_increasing(test: borrowed Test) throws {
    test.assertTrue(Array.isMonotonic([15.0, 25.0, 65.0, 65.0001]));
  }
  proc isMonotonic_decreasing(test: borrowed Test) throws {
    test.assertTrue(Array.isMonotonic([-15.0, -25.0, -65.0, -65.0001]));
  }
  proc isMonotonic_false(test: borrowed Test) throws {
    test.assertFalse(Array.isMonotonic([-15.0, 25.0, -65.0, -65.0001]));
  }


  proc linspace_increasing(test: borrowed Test) throws {
    var start = 0.0;
    var end = 5.0;
    var num = 6;
    test.assertEqual(Array.linspace(start, end, num),
                                    [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]);
    test.assertEqual(Array.linspace(start, end, num, endpoint=false),
                                    [0.0, 5.0/6, 5.0/3, 2.5, 10.0/3, 25.0/6]);
  }
  proc linspace_decreasing(test: borrowed Test) throws {
    var atol: real = 0.0;
    var rtol: real = 1e-15;
    var start = 6.0;
    var end = -2.0;
    var num = 8;
    var expectNoEnd = [6.0,34.0/7, 26.0/7, 18.0/7, 10.0/7, 2.0/7, -6.0/7, -2.0];
    test.assertLessThan(abs(Array.linspace(start, end, num) - expectNoEnd),
                        atol + rtol*(expectNoEnd));
    test.assertEqual(Array.linspace(start, end, num, endpoint=false),
                     [6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0, -1.0]);
  }
  proc linspace_complex(test: borrowed Test) throws {
    test.assertEqual(Array.linspace(-3.0+0i, 3.0+2i, 3),
                     [-3.0+0i, 0.0+1i, 3.0+2i]);
    test.assertEqual(Array.linspace(-3.0+0i, 3.0+2i, 3, endpoint=false),
                     [-3.0+0i, -1.0+2.0i/3, 1.0+4.0i/3]);
  }

  proc linspace_len2(test: borrowed Test) throws {
    test.assertEqual(Array.linspace(-16.0, -32.0, 2), [-16.0, -32.0]);
    test.assertEqual(Array.linspace(-16.0, -32.0, 2, endpoint=false),
                     [-16.0, -24.0]);
  }

  proc linspace_len1(test: borrowed Test) throws {
    test.assertEqual(Array.linspace(-16.0, -32.0, 1), [-16.0]);
    test.assertEqual(Array.linspace(-16.0, -32.0, 1, endpoint=false), [-16.0]);
  }

  proc empty_1(test: borrowed Test) throws {
    var b: [1..0] real;
    test.assertEqual(Array.empty(real), b);
    var c: [1..0] complex;
    test.assertEqual(Array.empty(complex), c);
    var d: [1..0] int(32);
    test.assertEqual(Array.empty(int(32)), d);
  }


  proc arange_range1D(test: borrowed Test) throws {
    test.assertEqual(Array.arange(-1..6), [-1, 0, 1, 2, 3, 4, 5, 6]);
    test.assertEqual(Array.arange(-1..6, int), [-1, 0, 1, 2, 3, 4, 5, 6]);
    test.assertEqual(Array.arange(1..5, real), [1.0, 2.0, 3.0, 4.0, 5.0]);
    test.assertEqual(Array.arange(7..10, complex), [7+0i, 8+0i, 9+0i, 10+0i]);
  }

  proc findAll_real(test: borrowed Test) throws {
    test.assertEqual(Array.findAll([1.0, 3.0, 5.0, 0.0, 5.0], 5.0), [2, 4]);
  }
  proc findAll_bool(test: borrowed Test) throws {
    test.assertEqual(Array.findAll([true, true, false], true), [0, 1]);
  }
  proc findAll_noneMatch(test: borrowed Test) throws {
    var empty: [1..0] int;
    test.assertEqual(Array.findAll([true, true, true], false), empty);
  }

  proc boolSlice_1d(test: borrowed Test) throws {
    test.assertEqual(Array.boolSlice([1, 2, 3], [true, false, false]), [1]);
    test.assertEqual(Array.boolSlice([1, 2, 3], [true, true, false]), [1, 2]);
    test.assertEqual(Array.boolSlice([1.3, 2.6, 3.5], [false, true, true]),
                     [2.6, 3.5]);
    test.assertEqual(Array.boolSlice(["foo", "bar"], [true, false]), ["foo"]);
  }
  proc boolSlice_noneMatch(test: borrowed Test) throws {
    var empty: [1..0] int;
    test.assertEqual(Array.boolSlice([1, 2], [false, false]), empty);
  }

  proc main() throws {
    UnitTest.main();
  }
}