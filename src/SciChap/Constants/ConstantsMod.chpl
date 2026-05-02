/*
  Private Constants module
*/
module ConstantsMod {
  import Math;

  record constant {
    /* full name of constant */
    param name: string;
    /* numerical value */
    param value: real;
    /* standard uncertainty */
    param uncertainty: real=0.0;
    /* unit */
    param unit: string="";

    /* relative standard uncertainty */
    proc relUncertainty do return value / uncertainty;
  }

  /* Pi (~3.14159) */
  param pi = Math.pi;
  // TODO: change to param when ** operator works for param
  /* Golden ratio (~1.618) */
  const goldenRatio: real = 0.5 + 0.5 * 5**0.5;

  /* quetta (Q) */
  param quetta: real = 1e30;
  /* ronna (R) */
  param ronna: real = 1e27;
  /* yotta (Y) */
  param yotta: real = 1e24;
  /* zetta (Z) */
  param zetta: real = 1e21;
  /* exa (E) */
  param exa: real = 1e18;
  /* peta (P) */
  param peta: real = 1e15;
  /* tera (T) */
  param tera: real = 1e12;
  /* giga (G) */
  param giga: real = 1e9;
  /* mega (M) */
  param mega: real = 1e6;
  /* kilo (k) */
  param kilo: real = 1e3;
  /* hecto (h) */
  param hecto: real = 1e2;
  /* deca (da) */
  param deca: real = 1e1;
  /* deci (d) */
  param deci: real = 1e-1;
  /* centi (c) */
  param centi: real = 1e-2;
  /* milli (m) */
  param milli: real = 1e-3;
  /* micro (\\mu) */
  param micro: real = 1e-6;
  /* nano (n) */
  param nano: real = 1e-9;
  /* pico (p) */
  param pico: real = 1e-12;
  /* femto (f) */
  param femto: real = 1e-15;
  /* atto (a) */
  param atto: real = 1e-18;
  /* zepto (z) */
  param zepto: real = 1e-21;
  /* yocto (y) */
  param yocto: real = 1e-24;
  /* ronto (r) */
  param ronto: real = 1e-27;
  /* quecto (q) */
  param quecto: real = 1e-30;


  // Binary prefixes
  param kibi: real = 2e10;
  param mebi: real = 2e20;
  param gibi: real = 2e30;
  param tebi: real = 2e40;
  param pebi: real = 2e50;
  param exbi: real = 2e60;
  param zebi: real = 2e70;
  param yobi: real = 2e80;


  // Length
  param inch: real = 0.0254;
  param foot: real = 12 * inch;
  param yard: real = 3 * foot;
  param mile: real = 5280 * foot;
  param mil: real = 0.001 * inch;
  param fermi: real = femto;
  param angstrom: real = 1e-10;
  param micron: real = micro;
  /* 1 nautical mile (M, NM, nmi) in meter */
  param nauticalMile: real = 1852;
  /* 1 astronomical unit (au) in meter */
  param astronomicalUnit: real = 149_597_870_700;
  /* 1 light-year (ly) in meter (~9.5e12)*/
  param lightYear: real = 9_460_730_472_580.8 * kilo;
  /* 1 parsec (pc) in meter (~31 petameter) */
  param parsec: real = 648_000 / pi * astronomicalUnit;


  // Area
  /* 1 acre (ac) in square meter */
  param acre: real = 4840 * yard;
  /* 1 hectare (ha) in square meter */
  param hectare: real = 1e4;


  // Volume
  /* 1 liter (L) to cubic meter */
  param liter: real = 0.001;
  /* 1 US gallon in cubic meter */
  param gallonUS: real = 231 * inch * inch * inch;
  /* 1 imperial gallon in cubic meter */
  param gallonImperial: real = 4.54609 * liter;
  param fluidOunceUS: real = gallonUS / 128;
  param fluidOunceImperial: real = gallonImperial / 160;


  // Mass
  /* 1 gram in kg */
  param gram: real = milli;
  /* 1 metric ton (tonne) in kg */
  param metricTon: real = kilo;
  /* 1 ounce in kg */
  param ounce: real = 28.349523125 * gram;
  /* 1 pound-mass (lbm) in kg (0.45359237 kg) */
  param poundMass: real = 16.0 * ounce;
  /* 1 grain in kg */
  param grain: real = 64.798_91 * milli;
  /* 1 slug in kg (~14.59 kg) */
  param slug: real = poundMass * 9.806_65 / foot;
  /* 1 stone in kg */
  param stone: real = 14 * poundMass;
  /* 1 short ton (US ton) in kg */
  param shortTon: real = 2000 * poundMass;
  /* 1 long ton in kg (1016.047 kg) */
  param longTon: real = 2240*poundMass;
  /* 1 Troy ounce (oz t) in kg */
  param troyOunce: real = 480 * grain;
  /* 1 Troy pound in kg */
  param troyPound: real = 12 * troyOunce;
  /* 1 carat (ct) in kg */
  param carat: real = 200 * milli;


  // Force
  /* 1 dyne in Newtons */
  param dyne: real = 1e-5;
  /* 1 kilogram-force (kgf) in Newton */
  param kilogramForce: real = 9.806_65;
  /* 1 pound-force (lbf) in Newton */
  param poundForce = poundMass * kilogramForce;


  // Pressure
  /* 1 atmosphere (atm) in Pascal */
  param atmosphere: real = 101_325;
  /* 1 bar in Pascal */
  param bar: real = 100_000;
  /* 1 Torr (torr, mmHg) in Pascal */
  param torr: real = atmosphere / 760;
  /* 1 psi in Pascal (~6894.76 Pa) */
  param psi: real = poundForce / (inch * inch);


  // Temperature
  /* 1 degree Fahrenheit in degrees Celsius */
  param fahrenheit: real = 9/5.0;
  /* 1 Rankine in Kelvin */
  param rankine: real = 9/5.0;
  /* 0 degrees Celsius in Kelvin (offset) */
  param zeroCelsius: real = 273.15;
  /* 0 degrees Fahrenheit in Rankine (offset) */
  param zeroFahrenheit: real = 459.67;


  // Energy
  /* calorie(Th) in joules */
  param calorieThermochemical: real = 4.184;
  /* calorie(IT) (International Steam Table, 1956) in joules */
  param calorieIT: real = 4.1868;
  /* erg in joules */
  param erg: real = 1e-7;
  /* BTU(Th) in joules (~1054.35 J) */
  param btuThermochemical: real =
    calorieThermochemical * kilo * poundMass / fahrenheit;
  /* BTU(IT) in joules (~1055.06 J) */
  param btuIT: real = calorieIT * kilo * poundMass / fahrenheit;


  // Power
  /* 1 imperial horsepower (hp, bhp) in Watts (~745.7 W) */
  param horsepowerImperial: real = 550 * poundForce * foot;
  /* 1 metric horsepower (cv, PS) in Watts (~735.5 W)*/
  param horsepowerMetric: real = 75 * kilogramForce;
  /* 1 electric horsepower (hp_E) in Watts */
  param horsepowerElectric: real = 746;


  // Angle
  /* 1 degree in radians */
  param degree: real = pi / 180.0;
  /* 1 arcminute in radians */
  param arcminute: real = degree / 60.0;
  /* 1 arcsecond in radians */
  param arcsecond: real = degree / 3600.0;


  // Time
  /* 1 minute in seconds */
  param minute: real = 60;
  /* 1 hour in seconds */
  param hour: real = 60 * minute;
  /* 1 day in seconds */
  param day: real = 24 * hour;
  /* 1 week in seconds */
  param week: real = 7 * day;
  /* 1 year in seconds */
  param year: real = 365 * day;
  /* 1 Julian year in seconds */
  param julianYear: real = 365.25 * day;


  // Speed
  /* 1 kmh in m/s */
  param kmh = kilo / hour;
  /* 1 mph in m/s */
  param mph = mile / hour;
  /* 1 knot in m/s */
  param knot = 1.852 * kmh;

}