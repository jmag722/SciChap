module GenerateCodata {
  import IO;
  import IO.FormattedIO.format;
  import Map;
  import Math.{isClose, pi};
  import Regex;
  use URL;

  /* map from CODATA year to the URL of the file */
  var codataMap = new Map.map(int, string);
  codataMap.add(2022, "https://pml.nist.gov/cuu/Constants/Table/allascii.txt");

  /* User-specified CODATA year */
  config const year:int = 2022;

  const sep = "  ";

  class DeveloperError: Error {
    proc init(msg:string) { super.init(msg); }
  }


  /* get information from CODATA file online

  :arg reader: open URL reader

  :returns: the number of CODATA constants, the number of header lines in the
            file, and the year of the file
  */
  proc getCodataInfo(const filename): (int, int, int) throws {
    var nlines = 0;
    var header = true;
    var nheader = 0;
    var yearRegex = new Regex.regex("\\d{4}");
    var year: int = -1;
    var reader = openUrlReader(filename);

    for line in reader.lines() {
      nlines += 1;

      if header {
        nheader += 1;
        var potentialMatch = yearRegex.search(line);
        if potentialMatch.matched {
          year = line[potentialMatch]:int;
        }
        if line.startsWith("----------") then header = false;
      }
    }
    var nconstants = nlines - nheader;
    return (nconstants, nheader, year);
  }

  /* record for storing CODATA values while processing
  */
  record codata {
    var n: int;
    var year: int;
    var file: string;  // filename or url
    var names: [0..#n] string;  // raw CODATA name
    var varnames: [0..#n] string;  // variable names for chapel module
    var values: [0..#n] string;
    var uncertainties: [0..#n] string;
    var units: [0..#n] string;
    var derived: [0..#n] bool;  // value should be computed for inc. precision
    var unitless: [0..#n] bool;  // no units
    var exact: [0..#n] bool;  // is value exact
  }

  /* Read online CODATA file

    :arg filename: name of file (url)

    :arg nconstants: number of constants to read

    :arg nheader: number of header lines to skip

    :arg yr: year of file

    :returns: CODATA object
  */
  proc readData(const filename, in nconstants: int, in nheader: int,
                in yr:int): codata throws {
    var reader = openUrlReader(filename);
    var cd = new codata(n=nconstants, year=yr, file=filename);

    for 1..#nheader do reader.readLine();

    for idx in 0..#nconstants {
      var line = reader.readLine();
      cd.names[idx] = line[0..<60].strip();
      cd.units[idx] = line[110..].strip();
      cd.unitless[idx] = cd.units[idx] == "";

      // create variable name, remove symbols
      cd.varnames[idx] = cd.names[idx];
      cd.varnames[idx] = cd.varnames[idx].replace(
        new Regex.regex("[ .()/,-]"), "_");
      // make varname camel case
      var varlist = cd.varnames[idx].split("_");
      cd.varnames[idx] = varlist[0].toLower();
      for v in varlist[1..] {
        cd.varnames[idx] += v.toTitle();
      }

      var val_tmp = line[60..<85].strip();
      val_tmp = val_tmp.replace(" ", "_");
      val_tmp = val_tmp.replace("_e", "e");
      cd.derived[idx] = val_tmp.contains("...");
      // don't replace ""..."" yet, if it's not caught downstream it'll
      // result in a build error
      cd.values[idx] = val_tmp;

      var u_tmp = line[85..<110].strip();
      u_tmp = u_tmp.replace(" ", "_");
      u_tmp = u_tmp.replace("_e", "e");
      cd.exact[idx] = u_tmp.contains("(exact)");
      cd.uncertainties[idx] = u_tmp.replace("(exact)", "0.0");
    }
    return cd;
  }

  /* Write CODATA into chapel module

     :arg cd: codata object
  */
  proc writeData(const cd: codata): void throws {
    const outFilename = "Codata" + cd.year:string + ".chpl";
    const outFile = IO.openWriter(outFilename);

    outFile.write("/* CODATA ", cd.year:string, "\n");
    outFile.write(sep, "from NIST, ", cd.file, "\n");
    outFile.write("*/\n\n");

    outFile.write("module Codata", cd.year:string, " {\n");
    outFile.write(sep, "import ConstantsMod.constant;\n\n");

    for idx in 0..#cd.n {
      // write comment describing variable, if exact or not
      outFile.write(sep, "/* ", cd.names[idx]);
      if cd.exact[idx] {
        outFile.write(" (exact)");
      }
      outFile.write(" */\n");

      outFile.write(sep, "const ", cd.varnames[idx], " = new constant(\n");
      outFile.write(sep, sep, "name=", "\"", cd.names[idx], "\"", ",\n");
      outFile.write(sep, sep, "value=", cd.values[idx]);
      // if something after `value`
      if !cd.exact[idx] || !cd.unitless[idx] {
        outFile.write(",\n");
      }

      if !cd.exact[idx] {
        outFile.write(sep, sep, "uncertainty=", cd.uncertainties[idx]);
      }
      if !cd.unitless[idx] {
        if !cd.exact[idx] then outFile.write(",\n");
        outFile.write(sep, sep, "unit=", "\"", cd.units[idx], "\"", "\n");
      } else {
        outFile.write("\n");
      }

      outFile.write(sep, ");\n\n");
    }

    outFile.write("}\n");
  }

  @chplcheck.ignore("LineLength")
  /* Compute derived CODATA quantities

    :arg cd: codata object
  */
  proc computeDerived(ref cd: codata) throws {
    var derivedMap = new Map.map(string, string);
    for (n, v) in zip(cd.names, cd.values) {
      derivedMap.add(n, v);
    }
    const h = derivedMap["Planck constant"]:real;
    const k = derivedMap["Boltzmann constant"]:real;
    const c = derivedMap["speed of light in vacuum"]:real;
    var hbar = h / (2 * pi);
    var e = derivedMap["elementary charge"]:real;
    var eV = derivedMap["electron volt"]:real;
    var kj90 = derivedMap["conventional value of Josephson constant"]:real;
    var rk90 = derivedMap["conventional value of von Klitzing constant"]:real;
    var rk = 2*pi*hbar / e**2;  // R_K
    var kj = 2 * e / h;  // K_J
    var na = derivedMap["Avogadro constant"]:real;  // N_A
    var gasRu = na * k;  // R (universal)
    var vm100 = gasRu*273.15/1e5;  // V_m
    var vm101325 = gasRu*273.15/101325; // V_m
    // Wien displacement constants - see "Wien's displacement law" on wikipedia
    var xpeak_lam = 4.965114231744276303;  // xpeak (wavelength)
    var xpeak_nu = 2.821439372122078893;  // xpeak (frequency)

    var fmtStr = "%.17r";
    for idx in 0..#cd.n {
      if cd.derived[idx] {
        select cd.names[idx] {
          when "atomic unit of action" do cd.values[idx] = fmtStr.format(hbar);
          when "Boltzmann constant in eV/K" do cd.values[idx] = fmtStr.format(k / e);
          when "Boltzmann constant in Hz/K" do cd.values[idx] = fmtStr.format(k / h);
          when "Boltzmann constant in inverse meter per kelvin" do cd.values[idx] = fmtStr.format(k / (h * c));
          when "conductance quantum" do cd.values[idx] = fmtStr.format(2*e**2 / (2*pi*hbar));
          when "conventional value of ampere-90" do cd.values[idx] = fmtStr.format(kj90 * rk90 / (kj * rk));
          when "conventional value of coulomb-90" do cd.values[idx] = fmtStr.format(kj90 * rk90 / (kj * rk));
          when "conventional value of farad-90" do cd.values[idx] = fmtStr.format(rk90 / rk);
          when "conventional value of henry-90" do cd.values[idx] = fmtStr.format(rk / rk90);
          when "conventional value of ohm-90" do cd.values[idx] = fmtStr.format(rk / rk90);
          when "conventional value of volt-90" do cd.values[idx] = fmtStr.format(kj90 / kj);
          when "conventional value of watt-90" do cd.values[idx] = fmtStr.format(kj90**2 * rk90 / (kj**2 * rk));
          when "electron volt-hertz relationship" do cd.values[idx] = fmtStr.format(eV / h);
          when "electron volt-inverse meter relationship" do cd.values[idx] = fmtStr.format(eV / (h * c));
          when "electron volt-kelvin relationship" do cd.values[idx] = fmtStr.format(eV / k);
          when "electron volt-kilogram relationship" do cd.values[idx] = fmtStr.format(eV / c**2);
          when "elementary charge over h-bar" do cd.values[idx] = fmtStr.format(e / hbar);
          when "Faraday constant" do cd.values[idx] = fmtStr.format(na*e);
          when "first radiation constant" do cd.values[idx] = fmtStr.format(2*pi*h*c**2);
          when "first radiation constant for spectral radiance" do cd.values[idx] = fmtStr.format(2*h*c**2);
          when "hertz-electron volt relationship" do cd.values[idx] = fmtStr.format(h / e);
          when "hertz-inverse meter relationship" do cd.values[idx] = fmtStr.format(1 / c);
          when "hertz-kelvin relationship" do cd.values[idx] = fmtStr.format(h / k);
          when "hertz-kilogram relationship" do cd.values[idx] = fmtStr.format(h / c**2);
          when "inverse meter-electron volt relationship" do cd.values[idx] = fmtStr.format(h*c/e);
          when "inverse meter-joule relationship" do cd.values[idx] = fmtStr.format(h*c);
          when "inverse meter-kelvin relationship" do cd.values[idx] = fmtStr.format(h*c/k);
          when "inverse meter-kilogram relationship" do cd.values[idx] = fmtStr.format(h / c);
          when "inverse of conductance quantum" do cd.values[idx] = fmtStr.format(2*pi*hbar / (2*e**2));
          when "Josephson constant" do cd.values[idx] = fmtStr.format(kj);
          when "joule-electron volt relationship" do cd.values[idx] = fmtStr.format(1/e);
          when "joule-hertz relationship" do cd.values[idx] = fmtStr.format(1/h);
          when "joule-inverse meter relationship" do cd.values[idx] = fmtStr.format(1 / (h*c));
          when "joule-kelvin relationship" do cd.values[idx] = fmtStr.format(1/k);
          when "joule-kilogram relationship" do cd.values[idx] = fmtStr.format(1/c**2);
          when "kelvin-electron volt relationship" do cd.values[idx] = fmtStr.format(k / e);
          when "kelvin-hertz relationship" do cd.values[idx] = fmtStr.format(k / h);
          when "kelvin-inverse meter relationship" do cd.values[idx] = fmtStr.format(k / (h*c));
          when "kelvin-kilogram relationship" do cd.values[idx] = fmtStr.format(k / c**2);
          when "kilogram-electron volt relationship" do cd.values[idx] = fmtStr.format(c**2/e);
          when "kilogram-hertz relationship" do cd.values[idx] = fmtStr.format(c**2 / h);
          when "kilogram-inverse meter relationship" do cd.values[idx] = fmtStr.format(c / h);
          when "kilogram-joule relationship" do cd.values[idx] = fmtStr.format(c**2);
          when "kilogram-kelvin relationship" do cd.values[idx] = fmtStr.format(c**2 / k);
          when "Loschmidt constant (273.15 K, 100 kPa)" do cd.values[idx] = fmtStr.format(na / vm100);
          when "Loschmidt constant (273.15 K, 101.325 kPa)" do cd.values[idx] = fmtStr.format(na / vm101325);
          when "mag. flux quantum" do cd.values[idx] = fmtStr.format(2*pi*hbar / (2*e));
          when "molar gas constant" do cd.values[idx] = fmtStr.format(gasRu);
          when "molar Planck constant" do cd.values[idx] = fmtStr.format(na * h);
          when "molar volume of ideal gas (273.15 K, 100 kPa)" do cd.values[idx] = fmtStr.format(vm100);
          when "molar volume of ideal gas (273.15 K, 101.325 kPa)" do cd.values[idx] = fmtStr.format(vm101325);
          when "natural unit of action" do cd.values[idx] = fmtStr.format(hbar);
          when "natural unit of action in eV s" do cd.values[idx] = fmtStr.format(hbar / e);
          when "Planck constant in eV/Hz" do cd.values[idx] = fmtStr.format(h / e);
          when "reduced Planck constant" do cd.values[idx] = fmtStr.format(hbar);
          when "reduced Planck constant in eV s" do cd.values[idx] = fmtStr.format(hbar / e);
          when "reduced Planck constant times c in MeV fm" do cd.values[idx] = fmtStr.format(hbar * c / (e * 1e6 * 1e-15));
          when "second radiation constant" do cd.values[idx] = fmtStr.format(h*c/k);
          when "Stefan-Boltzmann constant" do cd.values[idx] = fmtStr.format(pi**2/60 * k**4/(hbar**3 * c**2));
          when "von Klitzing constant" do cd.values[idx] = fmtStr.format(rk);
          when "Wien frequency displacement law constant" do cd.values[idx] = fmtStr.format(xpeak_nu * k / h);
          when "Wien wavelength displacement law constant" do cd.values[idx] = fmtStr.format(h*c/(k * xpeak_lam));
          otherwise halt("derived constant is not expected");
        }
        // check that this computed value matches what was read from the file
        var fileValue = derivedMap[cd.names[idx]].replace("...", ""):real;
        if !isClose(cd.values[idx]:real, fileValue, relTol=8e-10) {
          var errStr = "Developer, you derived this constant (%s) incorrectly.\n".format(cd.names[idx]);
          errStr += "Computed value: %.17r\n".format(cd.values[idx]:real);
          errStr += "Expected value: %.17r...".format(fileValue);
          throw new DeveloperError(errStr);
        }
      }
    }
  }

  proc createCodataFile(const yr: int): void throws {
    var filename = codataMap[yr];
    var (nconstants, nheader, year_check) = getCodataInfo(filename);
    if yr != year_check then throw new DeveloperError(
      "The year in the CODATA file header does not match the input year"
    );
    var cd = readData(filename, nconstants, nheader, yr);
    computeDerived(cd);
    writeData(cd);
  }

  proc main() throws {
    createCodataFile(year);
  }

}