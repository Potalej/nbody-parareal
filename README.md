# nbody-parareal

<p align="center">
    <img src="https://img.shields.io/github/languages/top/potalej/nbody-parareal?label=Fortran&logo=fortran&labelColor=%236f4c91&color=gray" alt="fortran">
    <img src="https://img.shields.io/github/languages/code-size/potalej/nbody-parareal?label=Size&labelColor=%236f4c91&color=gray" alt="Size">
    <img src="https://img.shields.io/github/last-commit/potalej/nbody-parareal?label=Modified&labelColor=%236f4c91&color=gray" alt="Last Commit">
    <img src="https://img.shields.io/github/issues/potalej/nbody-parareal?label=issues&labelColor=%236f4c91&color=gray" alt="issues">
</p>

A gravitational N-body solver using the parallel-in-time Parareal method!

---

## Requirements

- [gfortran](https://gcc.gnu.org/fortran/)
- [Make](https://www.gnu.org/software/make/)
- [CMake](https://cmake.org/)
- [OpenMP](https://www.openmp.org/)

This software also uses the [JSON-Fortran](https://github.com/jacobwilliams/json-fortran) API. If you are in an Anaconda environment you can [install locally](https://github.com/jacobwilliams/json-fortran#conda); in other cases our CMake will automatically compile it locally.

## ⚙️ Building

After cloning the repository, build using CMake and Make (or another build system like Ninja):

```
cmake -B build
cd build
make
```

### Debugging

To see the debug messages, just compile with the flag `-DPRINT_LEVEL`:
```
cmake -B build -DPRINT_LEVEL=2
```

### Precision

The default precision is double (64). To use single (32) or quadruple precision (128), use the flag `-DPRECISION`
```
cmake -B build -DPRECISION=32  # single
cmake -B build -DPRECISION=128 # quadruple
```

> Notice that using single or quad, the program will need to compile the JSON-Fortran library locally, and it may take a while.

---

## 🧮 Examples

There are four examples available, being two sequential and two parareal.

```
./nbody_parareal -r examples/sequential_random.json  # generate random initial values
./nbody_parareal -r examples/sequential_iv.json      # use defined initial values

./nbody_parareal -r examples/parareal_random.json    # generate random initial values
./nbody_parareal -r examples/parareal_iv.json        # use defined initial values
```
