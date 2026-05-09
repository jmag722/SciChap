module LinAlgTest {
  use UnitTest;

  import SciChap.LinAlg;

  proc tdma_2x2(test: borrowed Test) throws {
    var a = [1.3];
    var b = [3.0, 5.9];
    var c = [1.0];
    var d = [11.0, 12.0];
    LinAlg.tdma(lower=a, diagonal=b, upper=c, rhs=d);
    var expected = [3.2256097560975605, 1.3231707317073171];
    test.assertClose(d, expected, relTol=1e-15);
  }
  proc tdma_3x3(test: borrowed Test) throws {
    var a = [1.3, -0.5];
    var b = [3.0, 5.9, 7.9];
    var c = [1.0, 2.1];
    var d = [11.0, 12.0, 13.0];
    LinAlg.tdma(lower=a, diagonal=b, upper=c, rhs=d);
    var expected = [3.441790369979655, 0.6746288900610353, 1.688267651269686];
    test.assertClose(d, expected, relTol=1e-15);
  }

  proc main() throws {
    UnitTest.main();
  }


}