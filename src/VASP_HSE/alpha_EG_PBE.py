import matplotlib.pyplot as plt
from scipy.constants import pi
from scipy.stats import linregress
import json
import numpy

materials = json.load(open("../../data/HSE-data/all_data.json", "r"))
data = []
for m in materials:
    E_2D = m["gap_2D"]
    E_2D_direct = m["gap_2D_direct"]
    E_2D_min = m["gap_2D_min"]
    alpha_2D = m["alphax_2D"]
    # E_3D = m["gap_3D"]
    # L = m["L_3D"]
    # eps_x = m["epsx_3D"]
    n_2D = m["n_2D"]
    if None not in (E_2D, alpha_2D, n_2D):
        # alpha_3D = (eps_x - 1) * L / (4 * pi)
        data.append((E_2D, E_2D_direct, E_2D_min, alpha_2D, n_2D))

data = numpy.array(data)
# plt.plot(data[:, 0], data[:, 2], "o")  # 2D
# plt.plot(data[:, 1] ** -0.5, data[:, 4], "o")  # 2D
xx = numpy.linspace(0.5, 8)

plt.figure(figsize=(3.5, 3.5))
plt.style.use("science")
l, *_ =  plt.plot(data[:, 0], 1 / data[:, -2], "o",
                  alpha=0.5, label="HSE")
k, b, r, *_ = linregress(data[:, 0], 1/data[:, -2])
print(k, b, r)
plt.plot(xx, k * xx + b, "--", color=l.get_c())

l, *_ = plt.plot(data[:, 1], 1/data[:, -2], "o",
                 alpha=0.5, label="PBE-direct")
k, b, r, *_ = linregress(data[:, 1], 1/data[:, -2])
print(k, b, r)
plt.plot(xx, k * xx + b, "--", color=l.get_c())

l, *_ = plt.plot(data[:, 2], 1/data[:, -2], "o",
                 alpha=0.5, label="PBE-min")
k, b, r, *_ = linregress(data[:, 2], 1/data[:, -2])
print(k, b, r)
plt.plot(xx, k * xx + b, "--", color=l.get_c())
# plt.xlim(0, 2)

plt.xlabel("$E_{\\mathrm{g}}$ (eV)")
plt.ylabel("$(4 \\pi \\varepsilon_0) / \\alpha^{\parallel}$ ($\\AA^{-1}$)")
plt.legend()
plt.savefig("../../tmp_img/compare_alpha_Eg.svg")
    
