import matplotlib
import matplotlib.pyplot as plt
import numpy as np

fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8 = 506, 507, 508, 509, 510, 511, 512, 513
list_dist = [30, 70, 110, 150]
list_CIV_ew, list_CIV_std, list_MgII_ew, list_MgII_std = [], [], [], []

var_mode = "div"
file1 = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\Spectra-med\\"
image_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Bachelorarbeit\\Images\\Stacking\\"

list_fac, count_fac = [fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8], 1


bin_arr = np.linspace(1, 10, 10000, dtype=int)
bin_lst = bin_arr.tolist()
direc = {}
for ele in bin_lst:
    direc["x%s" % ele] = []
    direc["v%s" % ele] = []
    direc["z%s" % ele] = []
    direc["std%s" % ele] = []

for ele in list_fac:
    if ele != 0:
        x_fac, v_fac, z_fac, su_fac = "x" + str(count_fac), "v" + str(count_fac), "z" + str(count_fac), "std" + str(
            count_fac)
        for l in open(file1 + "Spec-med, bin, (" + str(ele) + ").txt", "r"):
            ls = [i for i in l.split()]
            direc["x" + str(count_fac)].append(float(ls[0]))
        for l in open(file1 + "Spec-med, median, (" + str(ele) + ").txt", "r"):
            ls = [i for i in l.split()]
            direc["v" + str(count_fac)].append(float(ls[0]))
        for l in open(file1 + "Spec-med, count, (" + str(ele) + ").txt", "r"):
            ls = [i for i in l.split()]
            direc["z" + str(count_fac)].append(float(ls[0]))
        for l in open(file1 + "Spec-med, std, (" + str(ele) + ").txt", "r"):
            ls = [i for i in l.split()]
            direc["std" + str(count_fac)].append(float(ls[0]))
        count_fac = count_fac + 1

# CIV
for ele6 in range(4):
    ele_06 = ele6 + 1
    ew = 0
    list_cont, list_flux, list_std, list_std_cont, list_std_Ilam, list_std_Icont = [], [], [], [], [], []
    for ele in direc["x" + str(ele_06)]:
        if 1538 <= ele + 1 <= 1545 or 1553 <= ele + 1 <= 1558:
            list_cont.append(direc["v" + str(ele_06)][int(ele + 1)])
        if 1538 <= ele + 1 <= 1545 or 1553 <= ele + 1 <= 1558:
            list_std_cont.append(direc["std" + str(ele_06)][int(ele + 1)])
        if 1546 <= ele + 1 <= 1552:
            list_flux.append(direc["v" + str(ele_06)][int(ele + 1)])
        if 1546 <= ele + 1 <= 1552:
            list_std.append(direc["std" + str(ele_06)][int(ele + 1)])
    cont = np.median(list_cont)
    std_cont = np.median(list_std_cont)

    count_std = 0
    for ele8 in list_std:
        list_std_Ilam.append((1 - 1/cont) * ele8)
        list_std_Icont.append((1 - list_flux[count_std]/(cont ** 2)) * std_cont)
        count_std = count_std + 1

    list_CIV_std.append(abs(sum(list_std_Ilam)) + abs(sum(list_std_Icont)))

    print(cont)
    for ele in list_flux:
        ew = ew + cont - ele
    print("CIV ew: ", ew / cont)
    list_CIV_ew.append(ew / cont)

# MgII
for ele6 in range(4):
    ele_06 = ele6 + 5
    ew = 0
    list_cont, list_flux, list_std, list_std_cont, list_std_Ilam, list_std_Icont = [], [], [], [], [], []
    for ele in direc["x" + str(ele_06)]:
        if 2784 <= ele + 1 <= 2792 or 2806 <= ele + 1 <= 2815:
            list_cont.append(direc["v" + str(ele_06)][int(ele + 1)])
        if 2784 <= ele + 1 <= 2792 or 2806 <= ele + 1 <= 2815:
            list_std_cont.append(direc["std" + str(ele_06)][int(ele + 1)])
        if 2793 <= ele + 1 <= 2806:
            list_flux.append(direc["v" + str(ele_06)][int(ele + 1)])
        if 2793 <= ele + 1 <= 2806:
            list_std.append(direc["std" + str(ele_06)][int(ele + 1)])
    cont = np.median(list_cont)
    std_cont = np.median(list_std_cont)

    count_std = 0
    for ele8 in list_std:
        list_std_Ilam.append((1 - 1 / cont) * ele8)
        list_std_Icont.append((1 - list_flux[count_std] / (cont ** 2)) * std_cont)
        count_std = count_std + 1

    list_MgII_std.append(abs(sum(list_std_Ilam)) + abs(sum(list_std_Icont)))

    print(cont)
    for ele in list_flux:
        ew = ew + cont - ele
    print("MgII ew: ", ew / cont)
    list_MgII_ew.append(ew / cont)

#
# fig1, ax = plt.subplots()
# plt.xlim(20, 300), plt.ylim(0.01, 3)
# plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("Äquivalentbreite EW [Å]")
# plt.errorbar(list_dist, list_CIV_ew, yerr=list_CIV_std, fmt="none", ecolor="xkcd:faded red", elinewidth=0.83, capsize=4, zorder=1)
# ax.scatter(list_dist, list_CIV_ew, color="xkcd:faded blue", s=30, label="CIV 1549"), ax.plot(list_dist, list_CIV_ew, linewidth=0.94, color="xkcd:faded blue")
# ax.set_xscale("log"), ax.set_yscale("log")
# ax.set_xticks([20, 50, 100, 200]), ax.set_yticks([0.01, 0.1, 1, 10])
# ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter()), ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
#
# plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("Äquivalentbreite EW [Å]")
# plt.errorbar(list_dist, list_MgII_ew, yerr=list_MgII_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
# ax.scatter(list_dist, list_MgII_ew, color="xkcd:faded green", s=30, label="MgII 2796"), ax.plot(list_dist, list_MgII_ew, linewidth=0.94, color="xkcd:faded green")
# plt.legend()
# plt.show()


fig1, ax = plt.subplots()
plt.grid(linestyle="-.", linewidth=0.35, zorder=1)
plt.title("EW CIV-Absorptionslinie")
plt.xlim(20, 300), plt.ylim(0.01, 3)
plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("Äquivalentbreite EW [Å]")
plt.errorbar(list_dist, list_CIV_ew, yerr=list_CIV_std, fmt="none", ecolor="xkcd:faded red", elinewidth=0.83, capsize=4, zorder=3)
ax.scatter(list_dist, list_CIV_ew, s=40, label="CIV 1549", zorder=4), ax.plot(list_dist, list_CIV_ew, linewidth=0.94, color="xkcd:faded blue", zorder=2)
ax.set_xscale("log"), ax.set_yscale("log")
ax.set_xticks([20, 50, 100, 200]), ax.set_yticks([0.01, 0.1, 1, 10])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter()), ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.savefig(image_path + "EW_CIV", dpi=200), plt.close()


fig2, ax = plt.subplots()
plt.grid(linestyle="-.", linewidth=0.35, zorder=1)
plt.title("EW MgII-Absorptionslinie")
plt.xlim(20, 300), plt.ylim(0.01, 10)
plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("Äquivalentbreite EW [Å]")
plt.errorbar(list_dist, list_MgII_ew, yerr=list_MgII_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
ax.scatter(list_dist, list_MgII_ew, s=40, label="MgII 2796"), ax.plot(list_dist, list_MgII_ew, linewidth=0.94, color="xkcd:faded blue")
ax.set_xscale("log"), ax.set_yscale("log")
ax.set_xticks([20, 50, 100, 200]), ax.set_yticks([0.01, 0.1, 1, 10])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter()), ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.savefig(image_path + "EW_MgII", dpi=200), plt.close()