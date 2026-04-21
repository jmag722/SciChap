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

  proc elementaryChargeHbar(test: borrowed Test) throws {
    test.assertClose(Constants.elementaryChargeHbar.value, 1.519267447e15,
                     absTol=1e-10);
  }
  proc magneticFluxQuantum(test: borrowed Test) throws {
    test.assertClose(Constants.magneticFluxQuantum.value, 2.067833848e-15,
                     absTol=1e-10);
  }
  proc conductanceQuantumInverse(test: borrowed Test) throws {
    test.assertClose(Constants.conductanceQuantumInverse.value, 12906.40372,
                     absTol=1e-6);
  }
  proc josephson(test: borrowed Test) throws {
    test.assertClose(Constants.josephson.value, 483597.8484e9, absTol=1e-5);
  }
  proc vonKlitzing(test: borrowed Test) throws {
    test.assertClose(Constants.vonKlitzing.value, 25812.80745, absTol=1e-6);
  }

  proc main() throws {
    UnitTest.main();
  }


}