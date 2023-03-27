import numpy, os
import matplotlib.pyplot as plt
from labutil.plugins.pwscf import run_qe_pwscf, PWscf_inparam, parse_qe_pwscf_output
from labutil.objects import Struc, Dir, ase2struc, Kpoints, Constraint, PseudoPotential
from ase.io import write
from ase import Atoms
from ase.build import bulk


def make_struc(alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    au_cell = bulk("Au", "fcc", a=alat)
    au_cell.set_atomic_numbers([79])
    print(au_cell, au_cell.get_atomic_numbers())
#     lattice = alat * numpy.identity(3)
#     symbols = ["Pb", "Ti", "O", "O", "O"]
#     sc_pos = [
#         [0, 0, 0],
#         [0.5, 0.5, 0.5 + displacement],
#         [0, 0.5, 0.5],
#         [0.5, 0, 0.5],
#         [0.5, 0.5, 0],
#     ]
#     perov = Atoms(symbols=symbols, scaled_positions=sc_pos, cell=lattice)
    # check how your cell looks like
    # write('s.cif', perov)
    structure = Struc(ase2struc(au_cell))
    print(structure.species)
    return structure


def compute_energy(alat, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    pseudopots = {
        "Au": PseudoPotential(
            ptype="uspp", element="Au", functional="LDA", name="Au.pz-d-rrkjus.UPF"
        ),
    }
    struc = make_struc(alat=alat)
    # fix the Pb and Ti atoms in place during relaxation
#     constraint = Constraint(atoms={"0": [0, 0, 0], "1": [0, 0, 0]})
    kpts = Kpoints(gridsize=[nk, nk, nk], option="automatic", offset=False)
    dirname = "Au_a_{}_ecut_{}_nk_{}".format(alat, ecut, nk)
    runpath = Dir(path=os.path.join(os.environ["WORKDIR"], "Lab4/Problem3", dirname))
    input_params = PWscf_inparam(
        {
            "CONTROL": {
                "calculation": "vc-relax",
                "pseudo_dir": os.environ["QE_POTENTIALS"],
                "outdir": runpath.path,
                "tstress": True,
                "tprnfor": True,
                "disk_io": "none",
            },
            "SYSTEM": {
                "ecutwfc": ecut,
                "ecutrho": ecut * 8,
                "occupations": "smearing",
                "smearing": "mp",
                "degauss": 0.02,
            },
            "ELECTRONS": {
                "diagonalization": "david",
                "mixing_beta": 0.7,
                "conv_thr": 1e-7,
            },
            "IONS": {"ion_dynamics": "bfgs"},
            "CELL": {"cell_dynamics": "bfgs"},
        }
    )

    output_file = run_qe_pwscf(
        runpath=runpath,
        struc=struc,
        pseudopots=pseudopots,
        params=input_params,
        kpoints=kpts,
        ncpu=2,
    )
    output = parse_qe_pwscf_output(outfile=output_file)
    return output


def lattice_scan():
#     nk = 4
    ecut = 40
    alat = 3
    nk_list = [3,5,7,9,11,13,15,17]
    energy_list = []
    for nk in nk_list:
        output = compute_energy(alat=alat, ecut=ecut, nk=nk)
        energy_list.append(output["energy"])
        print(output)
            
    plt.plot(nk_list, energy_list,'.-')
    plt.xlabel('K-points')
    plt.ylabel('Energy of Au (Ry)')
    plt.title('Convergence for cutoff')
    plt.show()
 


if __name__ == "__main__":
    # put here the function that you actually want to run
    lattice_scan()
