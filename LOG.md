# LOG.md - Development diary

## [2026-06-18] v0.0.2: Bux fixes
- The integrators pointers was pointing to the same object, so using the same method as coarse and as fine wasn't working. Now it is.

## [2026-06-15] v0.0.1: _lerigo_!

First implementation of the program!
- Sequential and Parareal methods available.
- Only Symplectic Euler and Velocity-Verlet integrators avaiable for now.
- A bit of user-friendly CLI.
- JSON-Fortran support.
- Simple, double and quadruple precision support.
- Output support.
- Some nice examples!