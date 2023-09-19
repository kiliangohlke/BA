import matplotlib.pyplot as plt
import numpy as np


def abline_func(double, line_name, line_wl, line_wl2, end_dashed, text_loc):
    axs[0, 0].vlines(x=line_wl, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
    axs[0, 1].vlines(x=line_wl, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
    axs[1, 0].vlines(x=line_wl, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
    axs[1, 1].vlines(x=line_wl, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
    axs[2, 0].vlines(x=line_wl, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
    axs[2, 1].vlines(x=line_wl, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
    axs[3, 0].vlines(x=line_wl, ymin=end_dashed, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
    axs[3, 1].vlines(x=line_wl, ymin=end_dashed, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)

    if fac4 == 0:
        if double == 0:
            axs[2, 0].text(line_wl, text_loc, line_name, ha="center", va="bottom", color="dimgrey", fontsize=9)
            axs[2, 1].text(line_wl, text_loc, line_name, ha="center", va="bottom", color="dimgrey", fontsize=9)
        if double == 1:
            axs[0, 0].vlines(x=line_wl2, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[0, 1].vlines(x=line_wl2, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[1, 0].vlines(x=line_wl2, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[1, 1].vlines(x=line_wl2, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[2, 0].vlines(x=line_wl2, ymin=end_dashed, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[2, 1].vlines(x=line_wl2, ymin=end_dashed, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[2, 0].text((line_wl + line_wl2) / 2, text_loc, line_name, ha="center", va="bottom", color="dimgrey",
                           fontsize=9)
            axs[2, 1].text((line_wl + line_wl2) / 2, text_loc, line_name, ha="center", va="bottom", color="dimgrey",
                           fontsize=9)
    else:
        if double == 0:
            axs[3, 0].text(line_wl, text_loc, line_name, ha="center", va="bottom", color="dimgrey", fontsize=9)
            axs[3, 1].text(line_wl, text_loc, line_name, ha="center", va="bottom", color="dimgrey", fontsize=9)
        if double == 1:
            axs[0, 0].vlines(x=line_wl2, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[0, 1].vlines(x=line_wl2, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[1, 0].vlines(x=line_wl2, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[1, 1].vlines(x=line_wl2, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[2, 0].vlines(x=line_wl2, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[2, 1].vlines(x=line_wl2, ymin=0.5, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[3, 0].vlines(x=line_wl2, ymin=end_dashed, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[3, 1].vlines(x=line_wl2, ymin=end_dashed, ymax=2, color="dimgrey", linestyle="dotted", linewidth=0.6)
            axs[3, 0].text((line_wl + line_wl2) / 2, text_loc, line_name, ha="center", va="bottom", color="dimgrey",
                           fontsize=9)
            axs[3, 1].text((line_wl + line_wl2) / 2, text_loc, line_name, ha="center", va="bottom", color="dimgrey",
                           fontsize=9)

# 500, 501, 502, 503, div, 0-100, 50-150, 100-200, 150-250, magi and magy < 26.5, mask y, ang 1, z 0-10, rej-out 2
# 506, 507, 508, 509, div, 0-60, 40-100, 80-140, 120-180, magi and magy < 27, mask y, ang 1, z 0-10, rej-out 2
# 510, 511, 512, 513, div, 0-60, 40-100, 80-140, 120-180, magi and magy < 26, mask y, ang 1, z 0-10, rej-out 2

fac1 = 506
# fac2, fac3, fac4 = 0, 0, 0
# fac5, fac6, fac7, fac8 = 0, 0, 0, 0
# fac2, fac3, fac4 = 500, 0, 0
# fac2, fac3, fac4 = 501, 502, 503
# fac5, fac6, fac7, fac8 = 500, 501, 502, 503
fac2, fac3, fac4 = 507, 508, 509
fac5, fac6, fac7, fac8 = 510, 511, 512, 513


var_mode = "div"
file1 = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\Spectra-med\\"
image_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Bachelorarbeit\\Images\\BV\\BV-python\\"
bin_arr = np.linspace(1, 10, 10000, dtype=int)
bin_lst = bin_arr.tolist()
direc = {}
for ele in bin_lst:
    direc["x%s" % ele] = []
    direc["v%s" % ele] = []
    direc["z%s" % ele] = []
    direc["su%s" % ele] = []

list_fac, count_fac = [fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8], 1

for ele in list_fac:
    if ele != 0:
        x_fac, v_fac, z_fac, su_fac = "x" + str(count_fac), "v" + str(count_fac), "z" + str(count_fac), "su" + str(
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
        for l in open(file1 + "Spec-med, sum, (" + str(ele) + ").txt", "r"):
            ls = [i for i in l.split()]
            direc["su" + str(count_fac)].append(float(ls[0]))
        count_fac = count_fac + 1

# CIV
ew, ew2 = 0, 0
list_cont, list_flux, list_flux2 = [], [], []
for ele in direc["x1"]:
    if 1538 <= ele + 1 <= 1545 or 1553 <= ele + 1 <= 1558:
        list_cont.append(direc["v1"][int(ele + 1)])
    if 1546 <= ele + 1 <= 1552:
        list_flux.append(direc["v1"][int(ele + 1)])
cont = np.mean(list_cont)
for ele in list_flux:
    ew = ew + cont - ele
print("CIV ew: ", ew / cont)

# MgII
ew, ew2 = 0, 0
list_cont, list_flux, list_flux2 = [], [], []
for ele in direc["x1"]:
    if 2784 <= ele + 1 <= 2792 or 2806 <= ele + 1 <= 2815:
        list_cont.append(direc["v1"][int(ele + 1)])
    if 2799 <= ele + 1 <= 2806:
        list_flux.append(direc["v1"][int(ele + 1)])
    if 2793 <= ele + 1 <= 2799:
        list_flux2.append(direc["v1"][int(ele + 1)])
cont = np.mean(list_cont)
for ele in list_flux:
    ew = ew + cont - ele
for ele in list_flux2:
    ew2 = ew2 + cont - ele
print("MgII ew: ", ew / cont)
print("MgII ew: ", ew2 / cont)

count_fac = 1
for ele_fac in list_fac:
    ele_0 = 0
    for ele in range(len(direc["z" + str(count_fac)])):
        ele_1 = ele - ele_0
        if direc["z" + str(count_fac)][ele_1] == 0:
            direc["x" + str(count_fac)].pop(ele_1)
            direc["v" + str(count_fac)].pop(ele_1)
            direc["z" + str(count_fac)].pop(ele_1)
            direc["su" + str(count_fac)].pop(ele_1)
            ele_0 = ele_0 + 1
    count_fac = count_fac + 1

plt.bar(direc["x1"], direc["z1"], width=0.5)
plt.title("Anzahl Werte pro lambda"), plt.xlabel("lambda"), plt.ylabel("Count pro lambda")
plt.show()

if var_mode == "sub":
    plt.xlabel("Ruhewellenlänge Vordergrundgalaxie [Å]"), plt.ylabel("Relativer Fluss"), plt.xlim(1000, 2000)
    plt.step(direc["x1"], direc["su1"], zorder=3, linewidth=0.6, label="median")

    if fac2 != 0:
        plt.step(direc["x2"], direc["su2"], color="indianred", zorder=2, linewidth=0.6, label="median")
    if fac3 != 0:
        plt.step(direc["x3"], direc["su3"], color="gold", zorder=1, linewidth=0.6, label="median")

    plt.axvline(x=1548.2, color="xkcd:jade green", linewidth=0.65, linestyle="dashed"), plt.text(1549.5, 0.85, "CIV", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.axvline(x=1550.8, color="xkcd:jade green", linewidth=0.65, linestyle="dashed")
    plt.axvline(x=2796.4, color="xkcd:jade green", linewidth=0.65, linestyle="dashed")
    plt.axvline(x=2803.5, color="xkcd:jade green", linewidth=0.65, linestyle="dashed"), plt.text(2799.95, 0.85, "MgII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.axvline(x=2600.2, color="xkcd:jade green", linewidth=0.65, linestyle="dashed"), plt.text(2590.0, 0.85, "FeII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.axvline(x=2586.7, color="xkcd:jade green", linewidth=0.65, linestyle="dashed"), plt.axvline(x=2344.2, color="xkcd:jade green", linewidth=0.65, linestyle="dashed"), plt.axvline(x=2374.5, color="xkcd:jade green", linewidth=0.65, linestyle="dashed"), plt.axvline(x=2382.8, color="xkcd:jade green", linewidth=0.65, linestyle="dashed"),
    plt.show()

if var_mode == "div":
    plt.xlabel("Ruhewellenlänge Vordergrundgalaxie [Å]"), plt.ylabel("Relativer Fluss"), plt.xlim(1000, 2000), plt.ylim(
        0.5, 1.2)
    plt.step(direc["x1"], direc["v1"], zorder=3, linewidth=0.6, label="median")

    if fac2 != 0:
        plt.step(direc["x2"], direc["v2"], color="indianred", zorder=2, linewidth=0.6, label="median")
    if fac3 != 0:
        plt.step(direc["x3"], direc["v3"], color="gold", zorder=1, linewidth=0.6, label="median")

    plt.vlines(x=1548.2, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1549.5, 0.85, "CIV", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1550.8, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65)
    plt.vlines(x=1333.1, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1333.1, 0.85, "CII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1656.9, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1656.9, 0.85, "CI", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1741.5, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1741.5, 0.85, "NiII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=2796.4, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65),
    plt.vlines(x=2803.5, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(2799.95, 0.85, "MgII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=2600.2, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(2600.2, 0.85, "FeII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=2586.7, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(2586.7, 0.85, "FeII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1854.7, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1854.7, 0.85, "AlII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1670.8, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1670.8, 0.85, "AlII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1808.0, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1808.0, 0.85, "SiII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1417.0, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1417.0, 0.85, "SiIII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1640.0, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1640.0, 0.85, "HeII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1808.0, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1808.0, 0.85, "SiII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1207.0, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1207.0, 0.85, "SiII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1501.8, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1501.8, 0.85, "SV", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1393.8, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1393.8, 0.85, "SiV", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1239.0, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1239.0, 0.85, "NV", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.vlines(x=1303.2, ymin=-0.5, ymax=1.1, color="xkcd:jade green", linewidth=0.65), plt.text(1303.2, 0.85, "SiII", ha="center", va="bottom", color="xkcd:jade green", fontsize=8)
    plt.show()

if var_mode == "div":
    if fac2 != 0:
        fig, axs = plt.subplots(4, 2, sharex="col", sharey="col", figsize=(11, 5), dpi=200)
        fig.subplots_adjust(hspace=0)
        cm = 1 / 2.54
        fig.text(0.5, 0.04, "Ruhewellenlänge Vordergrundgalaxie [Å]", ha="center", va="center")
        fig.text(0.06, 0.5, "Relativer Fluss", ha="center", va="center", rotation="vertical")

        abline_func(1, "CIV", 1547.7, 1550.3, 0.875, 0.79)
        abline_func(1, "MgII", 2795.5, 2802.7, 0.945, 0.90)
        # abline_func(1, "CIV", 1548.2, 1550.8, 0.875, 0.815)
        # abline_func(1, "MgII", 2796.4, 2803.5, 0.955, 0.93)
        # abline_func(0, "FeII", 2600.2, 0, 0.98, 0.97)
        # abline_func(0, "FeII", 2586.7, 0, 0.98, 0.97)
        # abline_func(0, "FeII", 2382.8, 0, 0.975, 0.96)
        # abline_func(0, "FeII", 2374.5, 0, 0.975, 0.96)
        # abline_func(0, "FeII", 2344.2, 0, 0.975, 0.96)
        # abline_func(0, "NiII", 1741.5, 0, 0.941, 0.927)
        # abline_func(0, "NiII", 1751.9, 0, 0.941, 0.927)
        # abline_func(1, "SiIV", 1393.8, 1402.8, 0.905, 0.875)
        # abline_func(0, "SiIII", 1417, 0, 0.905, 0.875)
        # abline_func(0, "HeII", 1640, 0, 0.941, 0.927)
        # abline_func(0, "CI", 1656.9, 0, 0.941, 0.927)
        # abline_func(0, "AlII", 1670.7, 0, 0.941, 0.927)

        plot1, = axs[0, 0].step(direc["x1"], direc["v1"], zorder=3, linewidth=0.6, label="median",
                                color="xkcd:faded blue")
        axs[0, 0].set_ylim(0.83, 1.05)
        axs[0, 0].set_xlim(1300, 2800)
        axs[0, 0].spines.bottom.set_visible(False)
        # axs[0, 0].set_yticks(np.arange(0.82, 1.06, 0.02))
        axs[0, 0].set_yticks(np.arange(0.74, 1.10, 0.06))

        plot2, = axs[0, 1].step(direc["x5"], direc["v5"], zorder=3, linewidth=0.6, label="median",
                                color="xkcd:faded blue")
        axs[0, 1].set_ylim(0.83, 1.05)
        axs[0, 1].set_xlim(1300, 2800)
        axs[0, 1].spines.bottom.set_visible(False)
        # axs[0, 1].set_yticks(np.arange(0.82, 1.06, 0.02))
        axs[0, 1].set_yticks(np.arange(0.82, 1.03, 0.03))

        plot3, = axs[1, 0].step(direc["x2"], direc["v2"], zorder=3, linewidth=0.6, label="median",
                                color="xkcd:faded red")
        axs[1, 0].set_ylim(0.83, 1.05)
        axs[1, 0].set_xlim(1300, 2800)
        axs[1, 0].spines.bottom.set_visible(False)
        axs[1, 0].spines.top.set_visible(False)
        # axs[1, 0].set_yticks(np.arange(0.82, 1.06, 0.02))
        axs[1, 0].set_yticks(np.arange(0.74, 1.10, 0.06))

        plot4, = axs[1, 1].step(direc["x6"], direc["v6"], zorder=3, linewidth=0.6, label="median",
                                color="xkcd:faded red")
        axs[1, 1].set_ylim(0.83, 1.05)
        axs[1, 1].set_xlim(1300, 2800)
        axs[1, 1].spines.bottom.set_visible(False)
        axs[1, 1].spines.top.set_visible(False)
        # axs[1, 1].set_yticks(np.arange(0.82, 1.06, 0.02))
        axs[1, 1].set_yticks(np.arange(0.82, 1.03, 0.03))

        plot5, = axs[2, 0].step(direc["x3"], direc["v3"], zorder=3, linewidth=0.6, label="median",
                                color="xkcd:faded green")
        axs[2, 0].set_ylim(0.83, 1.05)
        axs[2, 0].set_xlim(1300, 2800)
        axs[2, 0].spines.top.set_visible(False)
        axs[2, 0].spines.bottom.set_visible(False)
        # axs[2, 0].set_yticks(np.arange(0.82, 1.06, 0.02))
        axs[2, 0].set_yticks(np.arange(0.74, 1.10, 0.06))

        plot6, = axs[2, 1].step(direc["x7"], direc["v7"], zorder=3, linewidth=0.6, label="median",
                                color="xkcd:faded green")
        axs[2, 1].set_ylim(0.83, 1.05)
        axs[2, 1].set_xlim(1300, 2800)
        axs[2, 1].spines.top.set_visible(False)
        axs[2, 1].spines.bottom.set_visible(False)
        # axs[2, 1].set_yticks(np.arange(0.82, 1.06, 0.02))
        axs[2, 1].set_yticks(np.arange(0.82, 1.03, 0.03))

        plot7, = axs[3, 0].step(direc["x4"], direc["v4"], zorder=3, linewidth=0.6, label="median", color="xkcd:stone")
        axs[3, 0].set_ylim(0.83, 1.05)
        axs[3, 0].set_xlim(1300, 2800)
        axs[3, 0].spines.top.set_visible(False)
        # axs[3, 0].set_yticks(np.arange(0.82, 1.06, 0.02))
        axs[3, 0].set_yticks(np.arange(0.74, 1.10, 0.06))

        plot8, = axs[3, 1].step(direc["x8"], direc["v8"], zorder=3, linewidth=0.6, label="median", color="xkcd:stone")
        axs[3, 1].set_ylim(0.83, 1.05)
        axs[3, 1].set_xlim(1300, 2800)
        axs[3, 1].spines.top.set_visible(False)
        # axs[3, 1].set_yticks(np.arange(0.82, 1.06, 0.02))
        axs[3, 1].set_yticks(np.arange(0.82, 1.03, 0.03))

        axs[0, 0].legend((plot1, plot3, plot5, plot7), ("0<kpc<60", "40<kpc<100", "80<kpc<140", "120<kpc<180"), bbox_to_anchor=(1.1, 1.4), loc="upper center", frameon=False, ncol=4, markerscale=2)
        plt.plot(figsize=(15 * cm, 5 * cm))

        plt.show()

if var_mode == "sub":
    if fac2 != 0:
        fig, axs = plt.subplots(4, 2, sharex="col", sharey="col", figsize=(11, 5), dpi=200)
        fig.subplots_adjust(hspace=0)
        cm = 1 / 2.54
        fig.text(0.5, 0.04, "Ruhewellenlänge Vordergrundgalaxie [Å]", ha="center", va="center")
        fig.text(0.06, 0.5, "Relativer Fluss", ha="center", va="center", rotation="vertical")

        # abline_func(1, "CIV", 1548.2, 1550.8, 0.91, 0.88)
        # abline_func(1, "MgII", 2796.4, 2803.5, 0.941, 0.927)
        # abline_func(0, "FeII", 2600.2, 0, 0.96, 0.95)
        # abline_func(0, "FeII", 2586.7, 0, 0.96, 0.95)
        # abline_func(0, "FeII", 2382.8, 0, 0.975, 0.96)
        # abline_func(0, "FeII", 2374.5, 0, 0.975, 0.96)
        # abline_func(0, "FeII", 2344.2, 0, 0.975, 0.96)
        # abline_func(0, "NiII", 1741.5, 0, 0.941, 0.927)
        # abline_func(0, "NiII", 1751.9, 0, 0.941, 0.927)
        # abline_func(1, "SiIV", 1393.8, 1402.8, 0.905, 0.875)
        # abline_func(0, "SiIII", 1417, 0, 0.905, 0.875)
        # abline_func(0, "HeII", 1640, 0, 0.941, 0.927)
        # abline_func(0, "CI", 1656.9, 0, 0.941, 0.927)
        # abline_func(0, "AlII", 1670.7, 0, 0.941, 0.927)

        plot1, = axs[0, 0].step(direc["x1"], direc["su1"], zorder=3, linewidth=0.6, label="median", color="xkcd:faded blue")
        axs[0, 0].set_ylim(-300, 300)
        axs[0, 0].set_xlim(1300, 2800)
        axs[0, 0].spines.bottom.set_visible(False)
        axs[0, 0].set_yticks(np.arange(-300, 300, 50))

        plot2, = axs[0, 1].step(direc["x1"], direc["su1"], zorder=3, linewidth=0.6, label="median", color="xkcd:faded blue")
        axs[0, 1].set_ylim(-300, 300)
        axs[0, 1].set_xlim(1300, 2800)
        axs[0, 1].spines.bottom.set_visible(False)
        axs[0, 1].set_yticks(np.arange(-300, 300, 50))

        plot3, = axs[1, 0].step(direc["x2"], direc["su2"], zorder=3, linewidth=0.6, label="median",
                                color="xkcd:faded red")
        axs[1, 0].set_ylim(-300, 300)
        axs[1, 0].set_xlim(1300, 2800)
        axs[1, 0].spines.bottom.set_visible(False)
        axs[1, 0].spines.top.set_visible(False)
        axs[1, 0].set_yticks(np.arange(-300, 300, 50))

        plot4, = axs[1, 1].step(direc["x2"], direc["su2"], zorder=3, linewidth=0.6, label="median",
                                color="xkcd:faded red")
        axs[1, 1].set_ylim(-300, 300)
        axs[1, 1].set_xlim(1300, 2800)
        axs[1, 1].spines.bottom.set_visible(False)
        axs[1, 1].spines.top.set_visible(False)
        axs[1, 1].set_yticks(np.arange(-300, 300, 50))

        plot5, = axs[2, 0].step(direc["x3"], direc["su3"], zorder=3, linewidth=0.6, label="median",
                                color="xkcd:faded green")
        axs[2, 0].set_ylim(-300, 300)
        axs[2, 0].set_xlim(1300, 2800)
        axs[2, 0].spines.top.set_visible(False)
        axs[2, 0].spines.bottom.set_visible(False)
        axs[2, 0].set_yticks(np.arange(-300, 300, 50))

        plot6, = axs[2, 1].step(direc["x3"], direc["su3"], zorder=3, linewidth=0.6, label="median",
                                color="xkcd:faded green")
        axs[2, 1].set_ylim(-300, 300)
        axs[2, 1].set_xlim(1300, 2800)
        axs[2, 1].spines.top.set_visible(False)
        axs[2, 1].spines.bottom.set_visible(False)
        axs[2, 1].set_yticks(np.arange(-300, 300, 50))

        plot7, = axs[3, 0].step(direc["x4"], direc["su4"], zorder=3, linewidth=0.6, label="median", color="xkcd:stone")
        axs[3, 0].set_ylim(-300, 300)
        axs[3, 0].set_xlim(1300, 2800)
        axs[3, 0].spines.top.set_visible(False)
        axs[3, 0].set_yticks(np.arange(-300, 300, 50))

        plot8, = axs[3, 1].step(direc["x4"], direc["su4"], zorder=3, linewidth=0.6, label="median", color="xkcd:stone")
        axs[3, 1].set_ylim(-300, 300)
        axs[3, 1].set_xlim(1300, 2800)
        axs[3, 1].spines.top.set_visible(False)
        axs[3, 1].set_yticks(np.arange(-300, 300, 50))

        axs[0, 0].legend((plot1, plot3, plot5, plot7), ("0<kpc<100", "50<kpc<150", "100<kpc<200", "150<kpc<250"), bbox_to_anchor=(1.1, 1.4), loc="upper center", frameon=False, ncol=4, markerscale=2)
        plt.plot(figsize=(15 * cm, 5 * cm))

        plt.show()