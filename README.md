# SciChap

SciChap is a numerical methods library aimed at facilitating scientific computations with standardized, tested modules. It is akin and similar in design scope to Scipy, but written in Chapel.

## Why would I start this library

Chapel has [built-in support for calling Python libraries](https://chapel-lang.org/docs/modules/packages/Python.html), thus providing access to things like Scipy. Regardless, many other C/C++/Fortran libraries exist to provide these same capabilities I'm (re)implementing.

So, why?

Comes down to:

1. I love to code

2. I like Chapel and want to learn it better

3. I want to learn more numerical methods and parallel computing, and coding is a way I feel I've internalized and understood something

4. Writing a library like this with relatively discrete modules, rather than a large monolithic apps, allows me the flexibility to pick up and come back to where I left off even after weeks/months of being busy with other things.

## Roadmap

You'd be right to be skeptical about a programmer trying to write a library even close to the quality of Scipy. The repository now is pretty modest, but I envision it being added to over time, by me and anyone else who wants to contribute (stats module anyone?). Some things I'm planning on (because I want to learn more about them):

- [ ] TDMA/PCR and any other relevant algorithms for solving tridiagonal systems
- [ ] Spatial library
  - [X] K-D Tree
  - [ ] 2D and (someday) 3D Delaunay triangulation
- [ ] Integration
  - [ ] Gaussian quadrature
  - [x] trapezoidal/simpson
- [ ] root-finding
  - [ ] bracketed methods (midpoint, regula falsi, ...)
  - [ ] open and other hybrid approaches (newton, Brent, etc)
- [X] scientific constants
- [ ] statistics

# License

SciChap is licensed under BSD 3-Clause License. See the [LICENSE](LICENSE) file for more information.
