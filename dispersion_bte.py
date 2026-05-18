#!/usr/bin/env python3
import argparse
import os
import re
import sys
import numpy as np
import matplotlib.pyplot as plt
import seekpath
import spglib
from ase.data import atomic_numbers


# ShengBTE writes omega in rad/ps; nu[THz] = omega / (2*pi)
RAD_PS_TO_THZ = 1.0 / (2.0 * np.pi)
THZ_TO_CM1 = 33.356410

# Periodic table for element-symbol
ELEMENT_Z = {symbol: z for symbol, z in atomic_numbers.items() if z > 0}

def fnum(s):
    """Float that understands Fortran 'd'/'D' exponents."""
    return float(s.replace('d', 'e').replace('D', 'e'))

def floats_after(text, regex, n):
    m = re.search(regex, text, re.IGNORECASE)
    if not m:
        return None
    rest = text[m.end():].splitlines()[0]
    nums = re.findall(r'[-+]?\d*\.?\d+(?:[eEdD][-+]?\d+)?', rest)
    if len(nums) < n:
        return None
    return [fnum(x) for x in nums[:n]]



## Read crystal structure + mesh from the CONTROL file

def parse_control(filename='CONTROL'):
    with open(filename) as f:
        raw = f.read()
    text = re.sub(r'!.*$', '', raw, flags=re.MULTILINE)  # strip comments

    # lfactor: lattice unit in nm (default 1.0)
    m = re.search(r'\blfactor\s*=\s*([-+]?[\d.]+(?:[eEdD][-+]?\d+)?)',
                  text, re.IGNORECASE)
    lfactor_nm = fnum(m.group(1)) if m else 1.0

    # lattvec(:,i) - the i-th lattice vector as a column
    lattvec = np.zeros((3, 3))
    for i in range(1, 4):
        v = floats_after(text, rf'lattvec\s*\(\s*:\s*,\s*{i}\s*\)\s*=', 3)
        if v is None:
            raise ValueError(f"could not parse lattvec(:,{i}) in {filename}")
        lattvec[:, i - 1] = v
    # rows = lattice vectors, in Angstrom (lfactor is in nm)
    lattice = (lattvec * lfactor_nm * 10.0).T

    # elements
    m = re.search(r'\belements\s*=\s*((?:"[^"]+"\s*,?\s*)+)', text)
    if not m:
        raise ValueError(f"could not parse elements in {filename}")
    elements = re.findall(r'"([^"]+)"', m.group(1))

    # types (one integer per atom)
    m = re.search(r'\btypes\s*=\s*([\d ,\t]+)', text)
    if not m:
        raise ValueError(f"could not parse types in {filename}")
    types = [int(t) for t in re.findall(r'\d+', m.group(1))]

    # positions(:,i) - one column per atom
    positions = np.zeros((len(types), 3))
    for i in range(1, len(types) + 1):
        v = floats_after(text, rf'positions\s*\(\s*:\s*,\s*{i}\s*\)\s*=', 3)
        if v is None:
            raise ValueError(f"could not parse positions(:,{i}) in {filename}")
        positions[i - 1] = v

    # ngrid
    m = re.search(r'\bngrid\s*\(\s*:\s*\)\s*=\s*([\d ,\t]+)', text)
    if not m:
        raise ValueError(f"could not parse ngrid in {filename}")
    mesh = tuple(int(t) for t in re.findall(r'\d+', m.group(1))[:3])

    symbols = [elements[t - 1] for t in types]
    numbers = [ELEMENT_Z[s] for s in symbols]

    return {
        'lattice':   lattice,
        'positions': positions,
        'numbers':   numbers,
        'symbols':   symbols,
        'mesh':      mesh,
    }

def display_label(label):
    return 'G' if label == 'GAMMA' else label

# Accept G / g / Gamma / GAMMA as Gamma
def canonical_label(label):
    return 'GAMMA' if label.upper() in ('G', 'GAMMA') else label



def tex_label(label):
    if label == 'G':
        return r'$\Gamma$'
    if '_' in label:
        base, sub = label.split('_', 1)
        return rf'${base}_{{{sub}}}$'
    return label


# Get high-symmetry points expressed in the reciprocal basis
def get_high_sym_points(cell):
    sk = seekpath.get_path(cell, with_time_reversal=True)
    user_lattice = np.asarray(cell[0])
    std_lattice = np.asarray(sk['primitive_lattice'])

    # q_user = q_std @ B_std @ inv(B_user) with B = 2*pi * inv(A)^T
    B_user = 2.0 * np.pi * np.linalg.inv(user_lattice).T
    B_std = 2.0 * np.pi * np.linalg.inv(std_lattice).T
    M = B_std @ np.linalg.inv(B_user)

    coords = {label: np.asarray(c) @ M
              for label, c in sk['point_coords'].items()}

    suggested = []
    for a, b in sk['path']:
        if not suggested or suggested[-1] != display_label(a):
            suggested.append(display_label(a))
        suggested.append(display_label(b))

    return coords, suggested, sk

##

# Return a function lookup(ix, iy, iz) -> row index in BTE.omega
def build_lookup(cell, mesh, qpoints_sb):
    ir_map, ir_addr = spglib.get_ir_reciprocal_mesh(mesh, cell)
    nx, ny, nz = mesh

    # (ix, iy, iz) in [0, N) -> spglib's representative index
    grid_to_rep = {}
    for flat, addr in enumerate(ir_addr):
        key = (int(addr[0]) % nx, int(addr[1]) % ny, int(addr[2]) % nz)
        grid_to_rep[key] = int(ir_map[flat])

    # row in BTE.qpoints / BTE.omega
    rep_to_row = {}
    for row, q in enumerate(qpoints_sb):
        ix = int(round(q[0] * nx)) % nx
        iy = int(round(q[1] * ny)) % ny
        iz = int(round(q[2] * nz)) % nz
        rep = grid_to_rep[(ix, iy, iz)]
        rep_to_row.setdefault(rep, row)

    if len(rep_to_row) != len(set(ir_map.tolist())):
        print("  WARNING: spglib and ShengBTE disagree on the number of "
              "irreducible q-points -- check that CONTROL matches BTE.qpoints.")

    def lookup(ix, iy, iz):
        key = (int(ix) % nx, int(iy) % ny, int(iz) % nz)
        return rep_to_row[grid_to_rep[key]]

    return lookup


def build_dispersion(labels, coords, lookup, omega, mesh, lattice, npts_per_seg):
    nx, ny, nz = mesh
    B = 2.0 * np.pi * np.linalg.inv(lattice).T

    def trilinear(q):
        fx = q[0] * nx
        fy = q[1] * ny
        fz = q[2] * nz
        ix0 = int(np.floor(fx))
        iy0 = int(np.floor(fy))
        iz0 = int(np.floor(fz))
        tx = fx - ix0
        ty = fy - iy0
        tz = fz - iz0
        
        c000 = omega[lookup(ix0,     iy0,     iz0    )]
        c100 = omega[lookup(ix0 + 1, iy0,     iz0    )]
        c010 = omega[lookup(ix0,     iy0 + 1, iz0    )]
        c110 = omega[lookup(ix0 + 1, iy0 + 1, iz0    )]
        c001 = omega[lookup(ix0,     iy0,     iz0 + 1)]
        c101 = omega[lookup(ix0 + 1, iy0,     iz0 + 1)]
        c011 = omega[lookup(ix0,     iy0 + 1, iz0 + 1)]
        c111 = omega[lookup(ix0 + 1, iy0 + 1, iz0 + 1)]
        c00 = c000 * (1.0 - tx) + c100 * tx
        c01 = c001 * (1.0 - tx) + c101 * tx
        c10 = c010 * (1.0 - tx) + c110 * tx
        c11 = c011 * (1.0 - tx) + c111 * tx
        c0 = c00 * (1.0 - ty) + c10 * ty
        c1 = c01 * (1.0 - ty) + c11 * ty
        return c0 * (1.0 - tz) + c1 * tz

    dist, freq, qpts = [], [], []
    ticks = [0.0]
    tick_labels = [display_label(labels[0])]
    cum = 0.0

    for i in range(len(labels) - 1):
        p0 = coords[labels[i]]
        p1 = coords[labels[i + 1]]
        seg_len = float(np.linalg.norm((p1 - p0) @ B))

        ts = np.linspace(0.0, 1.0, npts_per_seg)
        for j, t in enumerate(ts):
            if i > 0 and j == 0:
                continue  # don't duplicate the junction with the previous segment..
            q = p0 + t * (p1 - p0)
            qpts.append(q)
            freq.append(trilinear(q))
            dist.append(cum + t * seg_len)

        cum += seg_len
        ticks.append(cum)
        tick_labels.append(display_label(labels[i + 1]))

    return (np.asarray(dist), np.asarray(freq), np.asarray(qpts),
            ticks, tick_labels)


def write_data_file(filename, dist, freq, qpts, ticks, tick_labels,
                    unit_tag, path_str):
    n_bands = freq.shape[1]
    with open(filename, 'w') as f:
        f.write("# Phonon dispersion from ShengBTE\n")
        f.write(f"# Path: {path_str}\n")
        f.write("# Tick positions:\n")
        for lbl, x in zip(tick_labels, ticks):
            f.write(f"#   {lbl:>5s}  d = {x:.6f}\n")
        cols = (['distance', 'qx', 'qy', 'qz'] +
                [f'omega_{b + 1}({unit_tag})' for b in range(n_bands)])
        f.write("#" + "  ".join(f"{c:>15s}" for c in cols) + "\n")
        for i in range(len(dist)):
            row = [f"{dist[i]:16.6f}",
                   f"{qpts[i, 0]:16.6f}",
                   f"{qpts[i, 1]:16.6f}",
                   f"{qpts[i, 2]:16.6f}"]
            row += [f"{freq[i, b]:16.6f}" for b in range(n_bands)]
            f.write(" ".join(row) + "\n")


def make_plot(dist, freq, ticks, tick_labels, ylabel, out_png):
    plt.rcParams.update({
        'font.size': 13,
        'axes.linewidth': 1.4,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.major.size': 5,
        'ytick.major.size': 5,
    })

    fig, ax = plt.subplots(figsize=(8.5, 6.0))
    for b in range(freq.shape[1]):
        ax.plot(dist, freq[:, b], color='royalblue', lw=1.6)

    for x in ticks:
        ax.axvline(x, color='black', lw=0.8)
    ax.axhline(0, color='gray', lw=0.6, ls='--')

    ax.set_xticks(ticks)
    ax.set_xticklabels([tex_label(l) for l in tick_labels])
    ax.set_xlim(ticks[0], ticks[-1])
    ax.set_xlabel('Wave Vector')
    ax.set_ylabel(ylabel)
    ax.tick_params(top=True, right=True)

    fig.tight_layout()
    fig.savefig(out_png, dpi=300)
    return fig


def path_tag(labels):
    """Turn ['GAMMA', 'X', 'K', 'GAMMA', 'L'] into 'G_X_K_G_L'."""
    return "_".join(display_label(l) for l in labels)


def main():
    parser = argparse.ArgumentParser(
        description="Plot phonon dispersion from ShengBTE output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--control', default='CONTROL')
    parser.add_argument('--qpoints', default='BTE.qpoints')
    parser.add_argument('--omega',   default='BTE.omega')
    parser.add_argument('--unit',    choices=['THz', 'cm-1'], default='THz',
                        help='frequency unit')
    parser.add_argument('--path',    default=None,
                        help='space-separated path labels; skips the prompt')
    parser.add_argument('--npts',    type=int, default=200,
                        help='sample points per segment')
    parser.add_argument('--output', default=None,
                        help='output directory')
    parser.add_argument('--no-show', action='store_true',
                        help='save the figure without opening a window')
    args = parser.parse_args()

    # Crystal ##
    print(f"== Crystal structure from {args.control} ==")
    cs = parse_control(args.control)
    cell = (cs['lattice'], cs['positions'], cs['numbers'])
    print(f"  Atoms ({len(cs['symbols'])}): {' '.join(cs['symbols'])}")
    print("  Lattice (Angstrom):")
    for v in cs['lattice']:
        print(f"     {v[0]:10.5f}  {v[1]:10.5f}  {v[2]:10.5f}")
    print(f"  Mesh: {cs['mesh'][0]} x {cs['mesh'][1]} x {cs['mesh'][2]}")

    # High-symmetry analysis ##
    print("\n== High-symmetry analysis (seekpath) ==")
    coords, suggested, sk = get_high_sym_points(cell)
    print(f"  Bravais lattice : {sk['bravais_lattice_extended']}")
    print(f"  Space group     : {sk['spacegroup_international']} "
          f"(#{sk['spacegroup_number']})")
    print("\n  High-symmetry points (in primitive reciprocal basis):")
    for lbl in sorted(coords.keys()):
        c = coords[lbl]
        print(f"     {display_label(lbl):>5s}  =  "
              f"({c[0]:8.4f}, {c[1]:8.4f}, {c[2]:8.4f})")
    print(f"\n  Seekpath's suggested path:  {' '.join(suggested)}")

    # path ##
    if args.path:
        labels_in = args.path.split()
        print(f"\n  Using path from --path: {' '.join(labels_in)}")
    else:
        labels_in = input(
            "\nType the path you want, separated by spaces "
            "(e.g.  G X K G L):\n> ").split()
    if len(labels_in) < 2:
        sys.exit("ERROR: need at least 2 labels for a path.")

    labels = [canonical_label(l) for l in labels_in]
    bad = [display_label(l) for l in labels if l not in coords]
    if bad:
        sys.exit(f"ERROR: label(s) not in the high-symmetry list: {bad}")

    # ShengBTE data ##
    print("\n== Loading ShengBTE output ==")
    qpoints_sb = np.loadtxt(args.qpoints, usecols=(3, 4, 5))
    omega = np.loadtxt(args.omega)
    print(f"  {qpoints_sb.shape[0]} IBZ q-points, {omega.shape[1]} bands")

    print("  Building full BZ -> IBZ symmetry map ...")
    lookup = build_lookup(cell, cs['mesh'], qpoints_sb)

    print("  Sampling dispersion along the path ...")
    dist, freq_rad, qpts, ticks, tick_labels = build_dispersion(
        labels, coords, lookup, omega, cs['mesh'], cs['lattice'], args.npts)

    # Units ##
    if args.unit == 'THz':
        freq = freq_rad * RAD_PS_TO_THZ
        ylabel = 'Frequency (THz)'
        unit_tag = 'THz'
    else:
        freq = freq_rad * RAD_PS_TO_THZ * THZ_TO_CM1
        ylabel = r'Frequency (cm$^{-1}$)'
        unit_tag = 'cm-1'

    # output filenames default to the path tag ##
    tag = path_tag(labels)
    output_dir = args.output if args.output else os.path.dirname(os.path.abspath(__file__))
    out_png = os.path.join(output_dir, f"phonon_dispersion_{tag}.png")
    out_dat = os.path.join(output_dir, f"dispersion_{tag}.out")

    path_str = " -> ".join(display_label(l) for l in labels)
    fig = make_plot(dist, freq, ticks, tick_labels, ylabel, out_png)
    write_data_file(out_dat, dist, freq, qpts, ticks, tick_labels,
                    unit_tag, path_str)

    print(f"\n  Saved plot : {out_png}")
    print(f"  Saved data : {out_dat}")
    print(f"  Path       : {path_str}")
    print(f"  Unit       : {unit_tag}")

    if not args.no_show:
        plt.show()
    else:
        plt.close(fig)

    print("\nDONE.")

if __name__ == '__main__':
    main()