# SciChap

SciChap is a numerical methods library aimed at facilitating scientific computations with standardized, tested modules. It is akin and similar in design scope to Scipy, but written in Chapel.

## Why would I start this library

Chapel has [built-in support for calling Python libraries](https://chapel-lang.org/docs/modules/packages/Python.html), thus providing access to things like Scipy. Regardless, many other C/C++/Fortran libraries exist to provide these same capabilities I'm (re)implementing.

So, why?

Comes down to:

1. I love to code

2. I like Chapel and want to learn it better

3. I want to learn more numerical methods and parallel computing, and coding is a way I feel I've internalized and understood something

Writing a library like this with relatively discrete modules, rather than a large monolithic app, allows me the flexibility to pick up and come back to where I left off even after weeks/months of being busy with other things.

## Roadmap

You'd be right to be skeptical about a programmer trying to write a library even close to the quality of Scipy. The repository now is pretty modest, but I envision it being added to over time, by me and anyone else who wants to contribute (stats module anyone?). Some things I'm planning on (because I want to learn more about them):

- [X] constants
- [ ] integration
  - [ ] Gaussian quadrature in multiple dimensions
  - [x] trapezoidal/simpson
- [ ] interpolation (nearest neighbor, linear, spline, RBF ...)
- [ ] linear algebra
  - [x] TDMA
  - [ ] PCR or some hybrid suitable for GPU
- [ ] root-finding
  - [ ] bracketed methods (bisection, regula falsi, ...)
  - [ ] open and other hybrid approaches (Newton, Brent, etc)
- [ ] spatial
  - [X] K-D Tree
  - [ ] 2D Delaunay triangulation
  - [ ] 3D Delaunay triangulation
- [ ] statistics

Other items:
- [ ] C/Fortran interfaces for library - would be amazing to have an easy-to-maintain, parallel numerical methods library available for other languages
- [ ] CMake support

## Installation

Currently using Mason for the build system.

To add to your application, put the following in your Mason.toml:

```toml
[dependencies]
SciChap = { git = "https://github.com/jmag722/SciChap", version = "0.1.0" }
```

To build chapel from source, see their [excellent docs](https://chapel-lang.org/docs/usingchapel/QUICKSTART.html), or my cheat sheet below.

### Chapel installation cheat sheet

Configure Chapel with something like the following:

```bash
# user must provide location of chapel source (CHPL_HOME) and compiler (LLVM here, see below)
CHPL_HOME="..."
export CHPL_HOME=${CHPL_HOME}

CHPL_BIN_SUBDIR=`"$CHPL_HOME"/util/chplenv/chpl_bin_subdir.py`
export PATH="$PATH":"$CHPL_HOME/bin/$CHPL_BIN_SUBDIR"
export MANPATH="$MANPATH":"$CHPL_HOME"/man

export CHPL_LLVM=system
LLVM_DIR="..."
LLVM_BIN="${LLVM_DIR}/bin"
export CHPL_HOST_CC="${LLVM_BIN}/clang"
export CHPL_HOST_CXX="${LLVM_BIN}/clang++"
export CHPL_TARGET_CC="${LLVM_BIN}/clang"
export CHPL_TARGET_CXX="${LLVM_BIN}/clang++"
export CHPL_LLVM_CONFIG="${LLVM_BIN}/llvm-config"
# if you want regex
export CHPL_RE2=bundled

INSTALL_DIR="${CHPL_HOME}"

./configure --prefix=${INSTALL_DIR}
```

Then to build:

```bash
n=12
make -j${n}
# `make install` doesn't work when CHPL_HOME is the source dir
make -j${n} chpldoc
make -j${n} mason
make -j${n} chplcheck
make -j${n} chpl-language-server
```

# License

SciChap is licensed under BSD 3-Clause License. See the [LICENSE](LICENSE) file for more information.
