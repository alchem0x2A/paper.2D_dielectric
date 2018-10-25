import ase.db
import warnings
import numpy
import matplotlib.pyplot as plt
from ase.data import covalent_radii
from scipy.stats import linregress
from scipy.optimize import curve_fit
import os, os.path
from scipy.constants import pi, epsilon_0
from gpaw_data import get_data
import scipy
import csv

"""
Extract the alpha from the HSE xlsx files
"""

db_file = "../../data/gpaw_data/c2db.db"
bulk_file = "../../data/2D-bulk/bulk.db"
if not os.path.exists(db_file):
    raise FileExistsError(("Please download the c2db data into ../../data/gpaw_data/ folder,"
                   "from https://cmr.fysik.dtu.dk/_downloads/c2db.db"))

def get_bulk(name, proto, id=None, method="gpaw"):
    # Get bulk properties
    if id is None:
        res = list(db.select(formula=name, prototype=proto))
        if len(res) == 0:
            return None
        r = res[0]
    else:
        r = db.get(id)    
    try:
        if method.lower() == "gpaw":
            L = r.bulk_L
            eps_para = (r.bulk_eps_x + r.bulk_eps_y) / 2
            eps_perp = r.bulk_eps_z
            e = r.gap_hse
        # VASP version below:
        elif method.lower() == "vasp":
            L = r.bulk_L_vasp
            eps_para = (r.bulk_eps_x_vasp + r.bulk_eps_y_vasp) / 2
            eps_perp = r.bulk_eps_z_vasp
            if r.bulk_gap_vasp < 0:
                r = r.gap_hse
            else:
                r = r.bulk_gap_vasep
        else:
            return None
        if eps_para < 0 or eps_perp < 0:
            return None
    except Exception:
        return None
    return L, eps_para, eps_perp, e

db = ase.db.connect(db_file)
bulk_db = ase.db.connect(bulk_file)

materials = []
eps_x = []
eps_z = []
alpha_x = []
alpha_z = []
Eg_HSE = []
thick = []

def get_thick(atom_row):
    pos = atom_row.positions[:, -1]
    diff = covalent_radii[atom_row.numbers]
    zmax = numpy.max(pos + diff) - numpy.min(pos - diff)
    return zmax

reader = csv.reader(open("../../data/HSE-data/2D_HSE.csv", encoding="utf8"))
next(reader)                    # skip line1
for row in reader:
    if row[4] != "":
        name, proto = row[: 2]
        print(name, proto)
        L, E, ex, ey, ez, *_ = map(float, row[2:])
        if ez < ex:
            eps_z.append(ez)
            materials.append("-".join((name, proto)))
            e_xy = numpy.sqrt(ex * ey)
            ax = (e_xy - 1) / (4 * pi) * L
            az = (1 - 1/ez) * L / (4 * pi)
            ax = max(1 / 1.2, ax)
            eps_x.append(e_xy); eps_z.append(ez)
            alpha_x.append(ax); alpha_z.append(az)
            Eg_HSE.append(E)
            mol = list(db.select(formula=name, prototype=proto))[0]
            thick.append(get_thick(mol))

print(len(alpha_x))
alpha_x = numpy.array(alpha_x)
alpha_z = numpy.array(alpha_z)
Eg_HSE = numpy.array(Eg_HSE)
thick = numpy.array(thick)

cond = numpy.where(Eg_HSE > 0.6)
Eg_HSE = Eg_HSE[cond]
alpha_x = alpha_x[cond]
alpha_z = alpha_z[cond]
thick = thick[cond]

colors = {"MoS2": "#AA0000", "CdI2": "#2C5AA0",
          "GaS": "#FFCC00", "BN": "#A05A2C",
          "P": "#447821", "CH": "#FF6600"}

cs = [colors[mat.split("-")[-1]] for mat in materials]

img_path = "../../tmp_img/"
plt.style.use("science")


gp_data = get_data()

import relation_2D3D as bulk

# x-direction
fig = plt.figure(figsize=(3, 3))
ax = fig.add_subplot(111)

cnt_r = numpy.array([3.94, 4.33, 5.12, 5.52, 6.30, 6.70, 7.49, 2.09])
cnt_x = numpy.array([174.7, 171.6, 292.4, 268.3, 445.5, 401.4, 651.1, 49.1])
cnt_z = numpy.array([10.3, 12.1, 15.8, 17.9, 22.4, 24.9, 30.2, 4.2])
# def get_cnt(Eg, R):
#     A = 52.0; B = 18.2
#     C = 9.5; D = 2.59; k = 5.0
#     ax = A + (R * B) / Eg** 2
#     az = (C + D * R ** 2) / k
#     return az / ax
A = 52.0; B = 18.2
# cnt_eg = numpy.linspace(0.1, 2, 8)
# cnt_r = numpy.linspace(2.0, 8.0, 8)
# cnt_data = numpy.array([(eg, get_cnt(eg, r)) for eg in cnt_eg for r in cnt_r])

cnt_eg = numpy.sqrt(cnt_r * B / (cnt_x - A))
# gpaw
ax.scatter(gp_data[2], gp_data[1] / gp_data[0], marker="^", alpha=0.1,
           c="#7289af",
           label="C2DB@2D")
# hse
ax.scatter(Eg_HSE, alpha_z / alpha_x, marker="^",
            edgecolors=None, alpha=0.1,
           label="VASP@2D")

ax.scatter(bulk.Eg_HSE, numpy.min([bulk.eps_z_3D[:, 0] / bulk.eps_x_3D[:, 0],
                                   bulk.eps_x_3D[:, 0] / bulk.eps_z_3D[:, 0]],
                                  axis=0),
           marker="s", alpha=0.1,
           label="VASP@3D")
ax.scatter(bulk.Eg_gpaw, numpy.min([bulk.eps_z_gpaw[:, 0] / bulk.eps_x_gpaw[:, 0],
                                    bulk.eps_x_gpaw[:, 0] / bulk.eps_z_gpaw[:, 0]],
                                   axis=0),
           marker="s",
           alpha=0.1,
           label="C2DB@3D")

# ax.scatter(cnt_data[:, 0], cnt_data[:, 1], alpha=0.5,
           # marker="*", label="CNT-1D")
ax.scatter(cnt_eg, cnt_z / cnt_x, alpha=0.2,
           marker="*", label="CNT-1D")


xx = yy = numpy.linspace(0, 8, 100)
ax.plot(xx, numpy.ones_like(xx), "--")

# ax.set_title("$y={0:.4f}x+{1:.4f},\ R^2={2:.4f}$".format(res.slope, res.intercept, res.rvalue))
ax.set_xlabel("$E_{\mathrm{g}}$")
ax.set_ylabel("Dielectric Anisotropy")
# ax.set_ylabel("$\\alpha_{zz}/(4\\pi \\varepsilon_0)$ ($\\AA$)")
ax.set_xlim(0, 8)
ax.set_ylim(0, 1.05)
ax.legend()
# ax.set_xticks([1,3,5,7])
fig.tight_layout()
fig.savefig(os.path.join(img_path, "alpha_alpha.svg"))
