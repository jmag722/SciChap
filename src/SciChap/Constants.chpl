/*
  Useful constants for unit conversions.

  The latest CODATA values will be available directly from this module. Older
  values are available only through directly importing the relevant submodule.

  .. include:: Constants/ConstantsMod.rst
    :start-after: START_HERE

*/
module Constants {
  include module ConstantsMod;
  include module Codata2022;
  include module Codata2018;
  include module Codata2014;
  include module Codata2010;

  public use this.ConstantsMod;
  public use this.Codata2022;
  public import this.Codata2018;
  public import this.Codata2014;
  public import this.Codata2010;
}