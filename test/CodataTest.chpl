module CodataTest {
  use UnitTest;

  import SciChap.Constants;

  proc speedOfLightInVacuum(test: borrowed Test) throws {
    test.assertClose(Constants.speedOfLightInVacuum.value, 3e8, relTol=1e-2);
  }

  proc main() throws {
    UnitTest.main();
  }


}