module GenerateCodata {
  import IO;
  import Map;
  import Regex;

  var codataMap = new Map.map(int, string);
  codataMap.add(2022, "https://pml.nist.gov/cuu/Constants/Table/allascii.txt");
  config const year:int = 2022;
  const sep = "  ";

  proc getCodataInfo(reader): (int, int, int) throws {
    var nlines = 0;
    var header = true;
    var nheader = 0;
    var yearRegex = new Regex.regex("\\d{4}");
    var year: int = -1;

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

  record codata {
    var n: int;
    var names: [0..#n] string;
    var varnames: [0..#n] string;
    var values: [0..#n] string;
    var uncertainties: [0..#n] string;
    var units: [0..#n] string;
    var derived: [0..#n] bool;
    var unitless: [0..#n] bool;
    var exact: [0..#n] bool;
  }

  proc readData(filereader, nconstants: int, nheader: int): codata throws {
    var cd = new codata(n=nconstants);

    for 1..#nheader do filereader.readLine();

    for idx in 0..#nconstants {
      var line = filereader.readLine();
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
      cd.values[idx] = val_tmp.replace("...", "");

      var u_tmp = line[85..<110].strip();
      u_tmp = u_tmp.replace(" ", "_");
      u_tmp = u_tmp.replace("_e", "e");
      cd.exact[idx] = u_tmp.contains("(exact)");
      cd.uncertainties[idx] = u_tmp.replace("(exact)", "0.0");
    }
    return cd;
  }

  proc writeData(cd: codata, yr: int): void throws {
    const outFilename = "Codata" + yr:string + ".chpl";
    const outFile = IO.openWriter(outFilename);

    outFile.write("module Codata", yr:string, " {\n");
    outFile.write(sep, "import ConstantsMod.constant;\n\n");

    for idx in 0..#cd.n {
      outFile.write(sep, "/* ", cd.names[idx], " */\n");
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

  proc computeDerived(cd: codata) throws {
    for idx in 0..#cd.n {
      if cd.derived[idx] {
        // TODO: cd must be a map to acceses vars by name
        // select cd.names[idx] {
        // when "atomic unit of action" do cd.values[idx] = ...;
        // otherwise halt("derived constant not yet supported");
      }
    }
  }

  proc createCodataFile(const yr: int): void throws {
    var file = codataMap[yr];
    use URL;
    var (nconstants, nheader, year_check) = getCodataInfo(openUrlReader(file));
    if yr != year_check then throw new Error(
      "Loaded file year does not match the codataMap (dev must update)"
    );
    var cd = readData(openUrlReader(file), nconstants, nheader);
    computeDerived(cd);
    writeData(cd, year);
  }

  proc main() throws {
    createCodataFile(year);
  }

}