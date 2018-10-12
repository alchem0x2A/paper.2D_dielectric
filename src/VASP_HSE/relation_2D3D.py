import ase.db
from ase.io.trajectory import Trajectory
from gpaw import GPAW
import warnings
import numpy
import matplotlib.pyplot as plt
from ase.data import covalent_radii
from scipy.stats import linregress
from scipy.optimize import curve_fit
import os, os.path
from scipy.constants import pi, epsilon_0
import scipy
import csv

"""
Extract the alpha from the HSE xlsx files
"""

db_file = "../../data/2D-bulk/bulk.db"


db = ase.db.connect(db_file)

materials = []
eps_x = []
eps_z = []
alpha_x = []
alpha_z = []
eps_x_3D = []                   
eps_z_3D = []
Eg_HSE = []
thick = []

def get_thick(atom_row):
    # Get thickness from 2D material
    pos = atom_row.positions[:, -1]
    diff = covalent_radii[atom_row.numbers]
    zmax = numpy.max(pos + diff) - numpy.min(pos - diff)
    return zmax

def get_bulk(name, proto, id=None):
    # Get bulk properties
    if id is None:
        res = list(db.select(formula=name, prototype=proto))
        if len(res) == 0:
            return None
        r = res[0]
    else:
        r = db.get(id)    
    try:
        L = r.bulk_L
        eps_para = (r.bulk_eps_x + r.bulk_eps_y) / 2
        eps_perp = r.bulk_eps_z
        if eps_para < 0 or eps_perp < 0:
            return None
    except Exception:
        return None
    return L, eps_para, eps_perp
    
    # dir_root = os.path.join("../../data/2D-bulk/{}-{}".format(name, proto))
    # gpw_file = os.path.join(dir_root, "gs.gpw")
    # traj_file = os.path.join(dir_root, "{}-{}_relax.traj".format(name, proto))
    # eps_file = os.path.join(dir_root, "eps_df.npz")
    # if (not os.path.exists(gpw_file)) or (not os.path.exists(eps_file)):
    #     return None
    # else:
    #     try:
    #         d = Trajectory(traj_file)[-1].cell[-1][-1]
    #     except Exception:
    #         return None
    #         # calc = GPAW(gpw_file)
    #         # d = calc.get_atoms().cell[-1][-1]
    #     f = numpy.load(eps_file)
    #     eps_xx = f["eps_x"][0].real
    #     eps_zz = f["eps_z"][0].real
    #     return d, eps_xx, eps_zz

reader = csv.reader(open("../../data/HSE-data/2D_HSE.csv", encoding="utf8"))
next(reader)                    # skip line1
for row in reader:
    if row[4] != "":
        name, proto = row[: 2]
        print(name, proto)
        L, E, ex, ey, ez, *_ = map(float, row[2:])
        if ez < ex:
            e_xy = numpy.sqrt(ex * ey)
            ax = (e_xy - 1) / (4 * pi) * L
            az = (1 - 1/ez) * L / (4 * pi)
            bulk_res = get_bulk(name, proto)
            if bulk_res is not None:
                materials.append("-".join((name, proto)))
                eps_x.append(e_xy); eps_z.append(ez)
                alpha_x.append(ax); alpha_z.append(az)
                L_3D, ex_3D, ez_3D = bulk_res
                print(L_3D, ex_3D, ez_3D)
                ex_simu = 1 + 4 * pi * ax / L_3D
                ez_simu = 1 / (1 - 4 * pi * az / L_3D)
                eps_x_3D.append((ex_3D, ex_simu))
                eps_z_3D.append((ez_3D, ez_simu))
                Eg_HSE.append(E)
            # mol = list(db.select(formula=name, prototype=proto))[0]
            # thick.append(get_thick(mol))

alpha_x = numpy.array(alpha_x)
alpha_z = numpy.array(alpha_z)
eps_x_3D = numpy.array(eps_x_3D)
eps_z_3D = numpy.array(eps_z_3D)
Eg_HSE = numpy.array(Eg_HSE)
# thick = numpy.array(thick)

eps_x_gpaw = []
eps_z_gpaw = []
for db_id in range(1, db.count() + 1):  # db index starts with 1
    mol = db.get(db_id)
    # if any(hasattr(mol, key) is False for key in ["alphax", "alphay", "alphaz",
                                         # "bulk_L", "bulk_eps_x", "bulk_eps_y",
                                         # "bulk_eps_z"]):
        # continue
    # if mol.bulk_calculated is False:
        # continue
    try:
        ax = (mol.alphax +  mol.alphay) / 2
        az = mol.alphaz
        L, ex, ez = get_bulk(None, None, db_id)
        ex_simu = 1 + 4 * pi * ax / L
        ez_simu = 1 / (1 - 4 * pi * az / L)
        eps_x_gpaw.append((ex, ex_simu))
        eps_z_gpaw.append((ez, ez_simu))
    except Exception:
        continue
eps_x_gpaw = numpy.array(eps_x_gpaw)
eps_z_gpaw = numpy.array(eps_z_gpaw)
print(eps_x_gpaw.shape, eps_z_gpaw.shape)

# cond = numpy.where(Eg_HSE > 0.6)
# Eg_HSE = Eg_HSE[cond]
# alpha_x = alpha_x[cond]
# alpha_z = alpha_z[cond]
# thick = thick[cond]

img_path = "../../tmp_img/"
plt.style.use("science")


def fit_func(x, a,b):
    return a+b / x



# x-direction
# MAE=numpy.mean(numpy.abs(eps_x_3D[:, 1][eps_x_3D[:, 0] < 30] - eps_x_3D[:, 0][eps_x_3D[:, 0]<30]))
# print(MAE)
plt.figure(figsize=(3, 2.8))
plt.plot(eps_x_gpaw[:, 1], eps_x_gpaw[:, 0], "s",
         alpha=0.2, color="grey")
plt.scatter(eps_x_3D[:, 1], eps_x_3D[:, 0],
            c=Eg_HSE,
            alpha=0.5,
            cmap="jet")
plt.plot(numpy.linspace(0, 30), numpy.linspace(0, 30), "--")
plt.xlim(0, 30)
plt.ylim(0, 30)
cb = plt.colorbar()
cb.ax.set_title("$E_{\mathrm{g}}$ (eV)")
plt.ylabel("eps bulk xx from DFT")
plt.xlabel("eps bulk xx from alpha")
plt.tight_layout()
plt.savefig(os.path.join(img_path, "eps_3D_xx.svg"))

# z-direction
plt.figure(figsize=(3, 2.8))
plt.plot(eps_z_gpaw[:, 1], eps_z_gpaw[:, 0],
         "s",
         alpha=0.2, color="brown")
plt.scatter(eps_z_3D[:, 1],
            eps_z_3D[:, 0], alpha=0.5,
            c=Eg_HSE,
            cmap="jet"
)
plt.plot(numpy.linspace(0, 10), numpy.linspace(0, 10), "--")
plt.ylim(0, 10)
plt.xlim(0, 10)
cb = plt.colorbar()
cb.ax.set_title("$E_{\mathrm{g}}$ (eV)")
plt.ylabel("eps bulk zz from DFT")
plt.xlabel("eps bulk zz from alpha")
plt.tight_layout()
plt.savefig(os.path.join(img_path, "eps_3D_zz.svg"))
