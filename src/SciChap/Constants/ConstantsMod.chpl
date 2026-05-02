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

  param pi = Math.pi;

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

}