module RootTest {
  use UnitTest;
  import SciChap.Root;

  proc bisect_maxIters(test: borrowed Test) throws {
    test.assertEqual(Root.bisect_maxIters(-5.0, 10.0, 3.3e-5), 19);
    test.assertEqual(Root.bisect_maxIters(25.0, 0.0, 1e-15), 55);
    test.assertEqual(Root.bisect_maxIters(5.0, 10.3, 1e-30), 103);
    test.assertEqual(Root.bisect_maxIters(5.0, 10.3, 1e-300), 999);
  }

  proc bisect_polyOdd(test: borrowed Test) throws {
    var func = proc(in x: real): real {return x**3 - x**5;};
    test.assertClose(Root.bisect(func, -0.5, 0.5), 0.0);
    test.assertClose(Root.bisect(func, 2.5, 0.8), 1.0);
    test.assertClose(Root.bisect(func, -2.5, -0.8), -1.0);
  }
  proc bisect_DoubleRoot(test: borrowed Test) throws {
    var func = proc(in x: real): real {return (x-2)**2 * (-2*x+5);};
    test.assertClose(Root.bisect(func, 1.8, 2.51), 2.5);
    var funcPrime = proc(in x: real): real {
      return 2*(x-2)*(-2*x+5) - 2*(x-2)**2;
    };
    // to find double root, run bisection on derivative
    test.assertClose(Root.bisect(funcPrime, 0.0, 2.2), 2.0);
  }

  proc main() throws {
    UnitTest.main();
  }
}