# Phonon Dispersion Generator for ShengBTE

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20275233.svg)](https://doi.org/10.5281/zenodo.20275233)


This Python script generates phonon dispersion curves directly from [ShengBTE](https://bitbucket.org/sousaw/shengbte) calculation files. The script reads the input `CONTROL` file together with the ShengBTE output files `BTE.qpoints` and `BTE.omega`, then produces both a phonon dispersion figure and a clean tabulated data file.

Crystallographic analysis is handled through [Seekpath](https://github.com/giovannipizzi/seekpath) and [spglib](https://spglib.readthedocs.io/), ensuring that the high-symmetry coordinates, labels, and reciprocal-space path are derived from the detected Bravais lattice of the input structure. 

---

## Inputs

Three files from a completed ShengBTE run are required:

| File          | Content                                                     |
| ------------- | ----------------------------------------------------------- |
| `CONTROL`     | lattice vectors, atomic positions, q-mesh                   |
| `BTE.qpoints` | irreducible-BZ q-point list                                 |
| `BTE.omega`   | phonon angular frequencies in rad/ps                        |

---

## Installation

Python ≥ 3.8 together with five packages:

```bash
pip install numpy matplotlib ase seekpath spglib
```

* [numpy](https://numpy.org/) and [matplotlib](https://matplotlib.org/) — numerical routines and plotting.
* [ase](https://wiki.fysik.dtu.dk/ase/) — element symbol to atomic number mapping.
* [seekpath](https://github.com/giovannipizzi/seekpath) — Bravais lattice detection and standard high-symmetry paths.
* [spglib](https://spglib.readthedocs.io/) — full-BZ to irreducible-BZ symmetry mapping.

---

## Usage

The script is invoked by passing the three input files and an output
directory:

```bash
python dispersion_bte.py \
    --control     CONTROL \
    --qpoints     BTE.qpoints \
    --omega       BTE.omega \
    --output      ./output/
```

After parsing the structure, the script prints the high-symmetry
points expressed in the primitive reciprocal basis of the input cell
and prompts for a path:

```
Type the path you want, separated by spaces (e.g.  G X K G L):
>  G X K G L
```

The prompt can be bypassed with `--path`:

```bash
python dispersion_bte.py --path "G X K G L"
```

Two output files are produced:

| File                              | Content                                                          |
| --------------------------------- | ---------------------------------------------------------------- |
| `phonon_dispersion_<PATH>.png`    | dispersion figure                                                |
| `dispersion_<PATH>.out`           | distance, qx, qy, qz, one column per band (units in the header)  |

The path is encoded in the filename, so successive runs along
different paths do not overwrite previous results.

---

## Examples

Two complete demonstrations are shipped in `examples/`. Running the
commands below should reproduce the figures stored in
`examples/*/output/`.

### Silicon (FCC, mesh: 48 × 48 × 48)

```bash
python ./dispersion_bte.py --control ./examples/Si/CONTROL --qpoints ./examples/Si/BTE.qpoints --omega ./examples/Si/BTE.omega --output ./examples/Si/output/
```

### Graphene (hexagonal, mesh: 250 × 250 × 1)

```bash
python ./dispersion_bte.py --control ./examples/Graphene/CONTROL --qpoints ./examples/Graphene/BTE.qpoints --omega ./examples/Graphene/BTE.omega --output ./examples/Graphene/output/
```

---

## Command-line options

```
python dispersion_bte.py [options]
```

| Option       | Default               | Description                                                  |
| ------------ | --------------------- | ------------------------------------------------------------ |
| `--control`  | `CONTROL`             | ShengBTE input file                                          |
| `--qpoints`  | `BTE.qpoints`         | ShengBTE q-point file                                        |
| `--omega`    | `BTE.omega`           | ShengBTE frequency file (rad/ps)                             |
| `--path`     | *interactive prompt*  | path labels, for example `--path "G X K G L"`                |
| `--unit`     | `THz`                 | `THz` or cm<sup>−1</sup> = `cm-1`                          |
| `--npts`     | `200`                 | sample points per segment                                    |
| `--output`   | script directory      | directory for the figure and `.out` file                     |
| `--no-show`  | (off)                 | save files without opening a plot window                     |

---

## Layout

```
.
├── dispersion_bte.py       
├── README.md               
└── examples/
    ├── Si/{CONTROL, BTE.qpoints, BTE.omega, output/}
    └── Graphene/{CONTROL, BTE.qpoints, BTE.omega, output/}
```

---

## References

* ShengBTE — <https://bitbucket.org/sousaw/shengbte>
* seekpath — <https://github.com/giovannipizzi/seekpath> · 
* spglib — <https://spglib.readthedocs.io/> · 
* ASE — <https://wiki.fysik.dtu.dk/ase/>
