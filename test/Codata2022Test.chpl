module Codata2022Test {
  use UnitTest;

  import SciChap.Constants;

  proc planckEv(test: borrowed Test) throws {
    test.assertClose(Constants.planckEv.value, 4.135667696e-15,
                     absTol=1e-10);
  }
  proc reducedPlanckC(test: borrowed Test) throws {
    test.assertClose(Constants.reducedPlanckC.value, 197.3269804,
                     absTol=1e-8);
  }

  proc main() throws {
    UnitTest.main();
  }


}