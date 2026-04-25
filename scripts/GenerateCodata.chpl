module GenerateCodata {
  import IO;
  import Regex;

  proc main() {
    var file = "codata2022.txt";
    var filereader = try! IO.openReader(file);

    var header = true;
    const sep = "  ";

    const outFile = try! IO.openWriter(
      "../src/SciChap/Constants/Codata2022.chpl"
    );
    try! outFile.write("module Codata2022 {\n");
    try! outFile.write(sep, "import ConstantsMod.constant;\n\n");

    for line in filereader.lines() {
      if header {
        if line.startsWith("----------") {
          header = false;
        }
        continue;
      }

      // split columns
      const name = line[0..<60].strip();
      var value = line[60..<85].strip();
      var uncertainty = line[85..<110].strip();
      const unit = line[110..].strip();

      // create variable name - remove symbols, make camel case
      var varname = name;
      varname = try! varname.replace(new Regex.regex("[ .()/,-]"), "_");
      var varlist = varname.split("_");
      varname = varlist[0].toLower();
      for v in varlist[1..] {
        varname += v.toTitle();
      }

      value = value.replace(" ", "_");
      value = value.replace("_e", "e");
      var isDerivedConst = false;
      if value.contains("...") then isDerivedConst = true;
      value = value.replace("...", "");

      uncertainty = uncertainty.replace(" ", "_");
      uncertainty = uncertainty.replace("_e", "e");
      const isExact = uncertainty == "(exact)";

      const unitless = unit == "";

      try! outFile.write(sep, "/* ", name, " */\n");
      try! outFile.write(sep, "const ", varname, " = new constant(\n");
      try! outFile.write(sep, sep, "name=", "\"", name, "\"", ",\n");

      try! outFile.write(sep, sep, "value=", value);
      if !isExact || !unitless {
        try! outFile.write(",\n");
      }

      if !isExact {
        try! outFile.write(sep, sep, "uncertainty=", uncertainty);
      }
      if !unitless {
        if !isExact then try! outFile.write(",\n");
        try! outFile.write(sep, sep, "unit=", "\"", unit, "\"", "\n");
      } else {
        try! outFile.write("\n");
      }

      try! outFile.write(sep, ");\n\n");
    }

    try! outFile.write("}\n");
  }

}