module IntegrationTest {
  use Math;
  use UnitTest;

  import SciChap.Integration;

  proc trapezoid_len1(test: borrowed Test) throws {
    test.assertEqual(Integration.trapezoid([1.0], [0.0]), 0.0);
  }

  proc trapezoid_len2(test: borrowed Test) throws {
    var dom: domain(1) = {6..7};
    var x: [dom] real = [1, 2];
    var y: [dom] real = [2, 4];
    test.assertEqual(Integration.trapezoid(y, x), 3.0);
  }

  proc trapezoid_linear(test: borrowed Test) throws {
    var dom: domain(1) = {0..<5};
    var x: [dom] real = [0, 1, 2, 3, 4];
    var y: [dom] real = [0, 3, 6, 9, 12];
    test.assertClose(Integration.trapezoid(y, x), 24.0, relTol=1e-15);
  }

  proc trapezoid_cos(test: borrowed Test) throws {
    var dom: domain(1) = {0..<4};
    var x: [dom] real = [0.0, 0.3, 0.6, 0.95];
    var y: [dom] real = cos(x);
    test.assertClose(Integration.trapezoid(y, x), 0.8066295622395069,
                     relTol=1e-15);
  }

  proc simpson_linearOdd(test: borrowed Test) throws {
    var dom: domain(1) = {0..5};
    var x: [dom] real = [1, 2, 3, 4, 5, 6];
    var y: [dom] real = [2, 4, 6, 8, 10, 12];
    test.assertClose(Integration.simpson(y=y, x=x), 35.0, relTol=1e-15);
  }
  proc simpson_linearEven(test: borrowed Test) throws {
    var dom: domain(1) = {0..6};
    var x: [dom] real = [1, 2, 3, 4, 5, 6, 8];
    var y: [dom] real = [2, 4, 6, 8, 10, 12, 16];
    test.assertClose(Integration.simpson(y=y, x=x), 63.0, relTol=1e-15);
  }
  proc simpson_quadraticEven(test: borrowed Test) throws {
    var n: int = 400;
    var dom: domain(1) = {0..#n};
    var x: [dom] real;
    var y: [dom] real;
    for i in dom {
      x[i] = i: real;
      y[i] = x[i]**2 + 3.0;
    }
    var actual: real = Integration.simpson(y=y, x=x);
    var expected: real = 1.0/3.0 * dom.last**3 + 3 * dom.last;
    test.assertClose(actual, expected, relTol=1e-15);
  }
  proc simpson_cubicOdd(test: borrowed Test) throws {
    var n: int = 525;
    var dom: domain(1) = {0..#n};
    var x: [dom] real;
    var y: [dom] real;
    for i in dom {
      x[i] = i: real;
      y[i] = x[i]**3 + 5*x[i]**2 + 6.0;
    }
    var actual: real = Integration.simpson(y=y, x=x);
    var xN = dom.last;
    var expected: real = 1.0/4.0 * xN**4 + 5.0/3.0 * xN**3 + 6.0 * xN;
    var atol: real = 0.0;
    var rtol: real = 1e-15;
    test.assertClose(actual, expected, relTol=1e-15);
  }

  proc main() throws {
    UnitTest.main();
  }

}


