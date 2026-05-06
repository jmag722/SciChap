/*

Private Constants module

TODO: figure out how to add sections so it's not a mass of constants in doc

..
  START_HERE
*/
module ConstantsMod {
  import Math;

  /* Constant type used for CODATA values */
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


  // SI prefixes

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

  /* one inch in meters */
  param inch: real = 0.0254;

  /* one foot in meters */
  param foot: real = 12 * inch;

  /* one yard in meters */
  param yard: real = 3 * foot;

  /* one mile in meters */
  param mile: real = 5280 * foot;

  /* one mil in meters */
  param mil: real = 0.001 * inch;

  /* one fermi in meters */
  param fermi: real = femto;

  /* one angstrom in meters */
  param angstrom: real = 1e-10;

  /* one micron in meters */
  param micron: real = micro;

  /* one nautical mile (M, NM, nmi) in meters */
  param nauticalMile: real = 1852;

  /* one astronomical unit (au) in meters */
  param astronomicalUnit: real = 149_597_870_700;

  /* one light-year (ly) in meters (~9.5e12 m)*/
  param lightYear: real = 9_460_730_472_580.8 * kilo;

  /* one parsec (pc) in meters (~31 petameter) */
  param parsec: real = 648_000 / pi * astronomicalUnit;


  // Area

  /* one acre (ac) in meters squared */
  param acre: real = 4840 * yard;

  /* one hectare (ha) in meters squared */
  param hectare: real = 1e4;


  // Volume

  /* one liter (L) in cubic meters */
  param liter: real = 0.001;

  /* one US gallon in cubic meters */
  param gallonUS: real = 231 * inch * inch * inch;

  /* one imperial gallon in cubic meters */
  param gallonImperial: real = 4.54609 * liter;

  /* one US fluid ounce in cubic meters */
  param fluidOunceUS: real = gallonUS / 128;

  /* one imperial fluid ounce in cubic meters */
  param fluidOunceImperial: real = gallonImperial / 160;


  // Mass

  /* one gram in kilograms */
  param gram: real = milli;

  /* one metric ton (tonne) in kilograms  */
  param metricTon: real = kilo;

  /* one ounce in kilograms  */
  param ounce: real = 28.349523125 * gram;

  /* one pound-mass (lbm) in kilograms  (0.45359237 kg) */
  param poundMass: real = 16.0 * ounce;

  /* one grain in kilograms  */
  param grain: real = 64.798_91 * milli;

  /* one slug in kilograms  (~14.59 kg) */
  param slug: real = poundMass * 9.806_65 / foot;

  /* one stone in kilograms  */
  param stone: real = 14 * poundMass;

  /* one short ton (US ton) in kilograms  */
  param shortTon: real = 2000 * poundMass;

  /* one long ton in kilograms  (1016.047 kg) */
  param longTon: real = 2240*poundMass;

  /* one Troy ounce (oz t) in kilograms  */
  param troyOunce: real = 480 * grain;

  /* one Troy pound in kilograms  */
  param troyPound: real = 12 * troyOunce;

  /* one carat (ct) in kilograms  */
  param carat: real = 200 * milli;


  // Force

  /* one dyne in Newtons */
  param dyne: real = 1e-5;

  /* one kilogram-force (kgf) in Newtons */
  param kilogramForce: real = 9.806_65;

  /* one pound-force (lbf) in Newtons */
  param poundForce = poundMass * kilogramForce;


  // Pressure

  /* one atmosphere (atm) in Pascals */
  param atmosphere: real = 101_325;

  /* one bar in Pascals */
  param bar: real = 100_000;

  /* one Torr (torr, mmHg) in Pascals (~133.3 Pa) */
  param torr: real = atmosphere / 760;

  /* one psi in Pascals (~6894.76 Pa) */
  param psi: real = poundForce / (inch * inch);


  // Temperature

  /* one degree Fahrenheit in degrees Celsius */
  param fahrenheit: real = 9/5.0;

  /* one Rankine in Kelvin */
  param rankine: real = 9/5.0;

  /* zero degrees Celsius in Kelvin (offset) */
  param zeroCelsius: real = 273.15;

  /* zero degrees Fahrenheit in Rankine (offset) */
  param zeroFahrenheit: real = 459.67;


  // Energy

  /* one calorie(Th) in joules */
  param calorieThermochemical: real = 4.184;

  /* one calorie(IT) (International Steam Table, 1956) in joules */
  param calorieIT: real = 4.1868;

  /* one erg in joules */
  param erg: real = 1e-7;

  /* one BTU(Th) in joules (~1054.35 J) */
  param btuThermochemical: real =
    calorieThermochemical * kilo * poundMass / fahrenheit;

  /* one BTU(IT) in joules (~1055.06 J) */
  param btuIT: real = calorieIT * kilo * poundMass / fahrenheit;


  // Power

  /* one imperial horsepower (hp, bhp) in Watts (~745.7 W) */
  param horsepowerImperial: real = 550 * poundForce * foot;

  /* one metric horsepower (cv, PS) in Watts (~735.5 W)*/
  param horsepowerMetric: real = 75 * kilogramForce;

  /* one electric horsepower (hp_E) in Watts */
  param horsepowerElectric: real = 746;


  // Angle

  /* one degree in radians */
  param degree: real = pi / 180.0;

  /* one arcminute in radians */
  param arcminute: real = degree / 60.0;

  /* one arcsecond in radians */
  param arcsecond: real = degree / 3600.0;


  // Time

  /* one minute in seconds */
  param minute: real = 60;

  /* one hour in seconds */
  param hour: real = 60 * minute;

  /* one day in seconds */
  param day: real = 24 * hour;

  /* one week in seconds */
  param week: real = 7 * day;

  /* one year in seconds */
  param year: real = 365 * day;

  /* one Julian year in seconds */
  param julianYear: real = 365.25 * day;


  // Speed
  /* one kilometer per hour in meters per second */
  param kmh = kilo / hour;

  /* one mile per hour in meters per second */
  param mph = mile / hour;

  /* one knot in meters per second */
  param knot = 1.852 * kmh;

}