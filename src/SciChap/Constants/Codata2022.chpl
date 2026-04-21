module Codata2022 {
  import SciChap.Constants.ConstantsMod.{constant};
  import SciChap.Constants.ConstantsMod;

  /* electron volt */
  const electronVolt = new constant(
    value=1.602_176_634e-19,
    unit="J",
    symbol="eV"
  );

  // Universal constants

  /* Speed of light in vacuum */
  const speedOfLight = new constant(
    value=299_792_458,
    unit="m s^-1",
    symbol="c"
  );

  /* Vacuum magnetic permeability */
  const permeabilityVacuum = new constant(
    value=1.256_637_061_27e-6,
    uncertainty=0.000_000_000_20e-6,
    unit="N A^-2",
    symbol="\\mu_0"
  );

  /* Vacuum electric permittivity */
  const permittivityVacuum = new constant(
    value=8.854_187_8188e-12,
    uncertainty=0.000_000_0014e-12,
    unit="F m^-1",
    symbol="\\varepsilon_0"
  );

  /* Characteristic impedance of vacuum */
  const characteristicImpedance = new constant(
    value=376.730_313_412,
    uncertainty=0.000_000_059,
    unit="ohm",
    symbol="Z_0"
  );

  /* Newtonian constant of gravitation */
  const gravitationalConstant = new constant(
    value=6.674_30e-11,
    uncertainty=0.000_15e-11,
    unit="m^3 kg^-1 s^-2",
    symbol="G"
  );

  /* Newtonian constant of gravitation over h-bar c */
  const gravitationalConstantHbar = new constant(
    value=6.708_83e-39,
    uncertainty=0.000_15e-39,
    unit="(GeV/c^2)^-2",
    symbol="G/(\\hbar c)"
  );

  /* Planck constant */
  const planck = new constant(
    value=6.626_070_15e-34,
    unit="J Hz^-1",
    symbol="h"
  );

  /* Planck constant in eV/Hz */
  const planckEv = new constant(
    value=planck.value / electronVolt.value,
    unit="eV Hz^-1",
    symbol="h"
  );

  /* Reduced Planck constant */
  const reducedPlanck = new constant(
    value=planck.value / (2 * ConstantsMod.pi),
    unit="J s",
    symbol="\\hbar"
  );

  /* Reduced Planck constant in eV s */
  const reducedPlanckEv = new constant(
    value=reducedPlanck.value / electronVolt.value,
    unit="eV s",
    symbol="\\hbar"
  );

  /* Reduced Planck constant times c in MeV fm */
  const reducedPlanckC = new constant(
    value=reducedPlanckEv.value * speedOfLight.value
         / (ConstantsMod.mega * ConstantsMod.femto),
    unit="MeV fm",
    symbol="\\hbar c"
  );

  /* Planck mass */
  const planckMass = new constant(
    value=2.176_434e-8,
    uncertainty=0.000_024e-8,
    unit="kg",
    symbol="m_P"
  );

  /* Planck mass energy equivalent in GeV */
  const planckMassEnergy = new constant(
    value=1.220_890e19,
    uncertainty=0.000_014e19,
    unit="GeV",
    symbol="m_P c^2"
  );

  /* Planck temperature */
  const planckTemperature = new constant(
    value=1.416_784e32,
    uncertainty=0.000_016e32,
    unit="K",
    symbol="T_P"
  );

  /* Planck length */
  const planckLength = new constant(
    value=1.616_255e-35,
    uncertainty=0.000_018e-35,
    unit="m",
    symbol="l_P"
  );

  /* Planck time */
  const planckTime = new constant(
    value=5.391_247e-44,
    uncertainty=0.000_060e-44,
    unit="s",
    symbol="t_P"
  );

}