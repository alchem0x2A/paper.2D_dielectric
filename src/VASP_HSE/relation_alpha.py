import ase.db
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

db_file = "../../data/gpaw_data/c2db.db"
if not os.path.exists(db_file):
    raise FileExistsError(("Please download the c2db data into ../../data/gpaw_data/ folder,"
                   "from https://cmr.fysik.dtu.dk/_downloads/c2db.db"))


db = ase.db.connect(db_file)

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

plt.figure(figsize=(7, 3.5))
plt.subplot(121)
plt.scatter(Eg_HSE, alpha_x, marker="o", alpha=0.5, c=cs)
def fit_func(x, a,b):
    return b / x

fit_param, cov = curve_fit(fit_func,
                      xdata=Eg_HSE, ydata=alpha_x,
                      p0=(0, 10))
print(fit_param)
xx = numpy.linspace(0.5, 8)
yy = fit_func(xx, *fit_param)
plt.plot(xx, yy, "--")
plt.xlabel("$E_{\\rm{g}}$ (eV)")
plt.ylabel("$\\alpha^{\\parallel}/(4 \\pi \\varepsilon_0)$ ($\\AA$)")

plt.subplot(122)
plt.scatter(Eg_HSE, alpha_z, marker="o", alpha=0.5, c=cs)
plt.xlabel("$E_{\\rm{g}}$ (eV)")
plt.ylabel("$\\alpha^{\\perp} / (4 \pi \\varepsilon_0)$ ($\\AA$)")

plt.tight_layout()
plt.savefig(os.path.join(img_path, "alpha_Eg_HSE.svg"))

plt.figure(figsize=(3.5, 3.5))
def fit_func(x, a,b):
    return b / x

fit_param, cov = curve_fit(fit_func,
                      xdata=Eg_HSE, ydata=alpha_x,
                      p0=(0, 10))
print(fit_param)
xx = numpy.linspace(0.5, 8)
yy = fit_func(xx, *fit_param)
plt.plot(Eg_HSE, alpha_x, "o", alpha=0.5)
plt.plot(xx, yy, "--")
plt.savefig(os.path.join(img_path, "alpha_fit_naive.pdf"))



# x-direction
plt.figure(figsize=(3.5, 3.5))
plt.scatter(Eg_HSE, 1 / (alpha_x), marker="o",
            edgecolors=None, alpha=0.5,
            c=cs)
res = linregress(x=Eg_HSE, y=1 / (alpha_x))
print(res)
xx = numpy.linspace(min(Eg_HSE), max(Eg_HSE))
yy = res.slope * xx + res.intercept
# for mat, x, y in zip(materials, Eg_HSE ** -p, alpha_x):
    # plt.text(x=x, y=y * 4 * pi, s=mat)
plt.plot(xx, yy, "--")
plt.title("$y={0:.4f}x+{1:.4f},\ R^2={2:.4f}$".format(res.slope, res.intercept, res.rvalue))
plt.xlabel("$E_{\\rm{g}}$ (eV)")
plt.ylabel("$1/\\alpha_{xx} \\varepsilon_0$ ($1/\\AA$)")
plt.savefig(os.path.join(img_path, "alpha_xx_1_Eg_HSE.svg"))

# z-direction
plt.figure(figsize=(3.5, 3.5))
plt.scatter(thick, alpha_z, marker="o",
            edgecolors=None,
            alpha=0.5, c=cs)
res = linregress(x=thick, y=alpha_z)
print(res)
xx = numpy.linspace(min(thick), max(thick))
yy = res.slope * xx + res.intercept
plt.plot(xx, yy, "--")
# plt.colorbar()
plt.title("$y={0:.2f}x+{1:.2f},\ R^2={2:.2f}$".format(res.slope, res.intercept, res.rvalue))
plt.xlabel("Thickness ($\\AA$)")
plt.ylabel("$\\alpha_{zz} / \\varepsilon_0$ ($\\AA$)")
plt.savefig(os.path.join(img_path, "alpha_zz_thick_HSE.svg"))

