import numpy as np
import math
from astropy.io import fits
import matplotlib.pyplot as plt
from mpdaf.obj import Cube

# var_dist_kpc_min = 0 var_dist_kpc = 250 var_z_max = 12 var_snr = 15 var_err = 0.5 var_mag = 25
# remove outlier, var_err
# try var_mag = 28

var_dist_kpc_min = 0
var_dist_kpc = 300
var_z_max = 0.91
var_snr = 5
var_err = 100
var_mag = 100


image_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\Images\\w7\\"
fit_main_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\UDF_HUDF_Catalogs_FULL\\dr2_main_09.fits"
cube_UDF_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\DATACUBE_UDF-10.fits"
cube_MXDF_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\DATACUBE_MXDF.fits"
cube_MOSAIC_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\DATACUBE_UDF-MOSAIC.fits"

# function kpc/arcsec at redshift z, conv_func(z)
a, b, c, d, e, f, g, h, o, j, k, l, m, n = [-24.862278979480052, 21.267147446426016, 0.005306164849709134, 19.600027472637663, -11.409572942008737, 4.751202039434343, -1.3170007020351893, 0.2021612052937083, -0.007249014401993737, 0.0016127139142351095, -0.00017642802527037128, 9.986804760095554e-06, -2.1958087984210742e-07, -1.148703972276188e-09]

def conv_func(x):
    return n * x**14 + m * x**13 + l * x**12 + k * x**11 + j * x**10 + o * x**9 + h * x**7 + g * x**6 + f * x**5 + e * x**4 + d * x**3 + a * x**2 + b * x + c

# def cubes
cube_UDF, cube_MXDF, cube_MOSAIC = Cube(cube_UDF_path, ext=1), Cube(cube_MXDF_path, ext=1), Cube(cube_MOSAIC_path, ext=1)
# open region and fits file
fits_info = fits.open(fit_main_path)

# right ascension
ra = fits_info[1].data['RA']
# declination
dec = fits_info[1].data['DEC']
# red shift
z = fits_info[1].data['Z']
# max flux
flux = fits_info[1].data["FLUX_MAX"]
# H-alpha emission flux line
h_a = fits_info[1].data["HALPHA_EMI_FLUX"]
# signal-to-noise Ha
h_a_snr = fits_info[1].data["HALPHA_EMI_SNR"]
# H-alpha absorption flux
ha_ab = fits_info[1].data["HALPHA_ABS_FLUX"]
# H-beta emission flux line
h_b = fits_info[1].data["HBETA_EMI_FLUX"]
# signal-to-noise Hb
h_b_snr = fits_info[1].data["HBETA_EMI_SNR"]
# H-beta absorption flux
hb_ab = fits_info[1].data["HBETA_ABS_FLUX"]
# H_gamma emission flux line
h_g = fits_info[1].data["HGAMMA_EMI_FLUX"]
# signal-to-noise Ha
h_g_snr = fits_info[1].data["HGAMMA_EMI_SNR"]
# H-Beta absorption flux
hg_ab = fits_info[1].data["HGAMMA_ABS_FLUX"]
# mag in F606W
mag606 = fits_info[1].data["MAG_F606W"]
# mag error in F606W
mag606_err = fits_info[1].data["MAGERR_F606W"]
# mag in F775W
mag775 = fits_info[1].data["MAG_F775W"]
# mag error in F775W
mag775_err = fits_info[1].data["MAGERR_F775W"]
# mag in F850LP
mag850 = fits_info[1].data["MAG_F850LP"]
# mag error in F775W
mag850_err = fits_info[1].data["MAGERR_F850LP"]

count_pairs, count_hahb, count_hghb = 0, 0, 0

list_HaHb_value, list_HaHb_dist_kpc, list_HgHb_value, list_HgHb_dist_kpc = [], [], [], []
list_fg_HaHb_value, list_fg_HaHb_dist_kpc, list_fg_HgHb_value, list_fg_HgHb_dist_kpc = [], [], [], []
list_F6F7_value, list_F6F7_dist_kpc, list_F6F8_value, list_F6F8_dist_kpc = [], [], [], []
list_add = [14, 202, 207, 417]

list_F6F7_value2, list_F6F7_value3, list_HgHb_value2, list_HaHb_value2 = [], [], [], []

# goes through list table, saves right ascension and declination for every object in var r, d
for i in range(len(ra)):
    ri, di, zi, flux_i = ra[i], dec[i], float(z[i]), flux[i]
    for y in range(len(ra)):
        # goes through list table, saves right ascension and declination, flux, emission lines for
        ry, dy, zy, mag606_y = ra[y], dec[y], float(z[y]), mag606[y]

        if i != y:
            if zi > zy and mag606_y < var_mag:
                # if 0.01 < zi < var_z_max and zy > 0.01:


                if 0.01 < zi < var_z_max and zy > 0.01:
                    dist = ((ry - ri) ** 2 + (dy - di) ** 2) ** (1 / 2)
                    dist_kpc = dist * 3600 * conv_func(zy)

                    if var_dist_kpc_min < dist_kpc <= var_dist_kpc:
                        print("")
                        print("i: " + str(i+1), "y: " + str(y+1), "dist-kpc", dist_kpc, dist)
                        count_pairs = count_pairs + 1

                        h_a_i, h_b_i, h_g_i = h_a[i], h_b[i], h_g[i],
                        ha_snr_i, hb_snr_i, hg_snr_i = h_a_snr[i], h_b_snr[i], h_g_snr[i]
                        ha_ab_i, hb_ab_i, hg_ab_i = ha_ab[i], hb_ab[i], hg_ab[i]

                        h_a_y, h_b_y, h_g_y = h_a[y], h_b[y], h_g[y],
                        ha_snr_y, hb_snr_y, hg_snr_y = h_a_snr[y], h_b_snr[y], h_g_snr[y]
                        ha_ab_y, hb_ab_y, hg_ab_y = ha_ab[y], hb_ab[y], hg_ab[y]

                        # Blamer Decrement
                        if var_snr < ha_snr_i < 1e+15 and var_snr < hb_snr_i < 1e+15:
                            if (h_b_i + hb_ab_i > 0 and h_a_i + ha_ab_i > 0) or i in list_add:
                                count_hahb = count_hahb + 1
                                # list_HaHb_value.append(h_a_i/(h_b_i * 2.85))
                                list_HaHb_value.append(-2.5 * math.log(h_a_i / (h_b_i * 2.85)))
                                list_HaHb_dist_kpc.append(dist_kpc)
                        if var_snr < ha_snr_y < 1e+15 and var_snr < hb_snr_y < 1e+15:
                            if (h_b_y + hb_ab_y > 0 and h_a_y + ha_ab_y > 0) or i in list_add:
                                list_fg_HaHb_value.append(-2.5 * math.log(h_a_y / (h_b_y * 2.85)))
                                list_fg_HaHb_dist_kpc.append(dist_kpc)

                        if var_snr < hg_snr_i < 1e+15 and var_snr < hb_snr_i < 1e+15:
                            if (h_b_i + hb_ab_i > 0 and h_g_i + hg_ab_i > 0) or i in list_add:
                                count_hghb = count_hghb + 1
                                # list_HgHb_value.append(h_g_i/(h_b_i * 0.469))
                                list_HgHb_value.append(-2.5 * math.log(h_g_i / (h_b_i * 0.47)))
                                list_HgHb_dist_kpc.append(dist_kpc)
                        if var_snr < hg_snr_y < 1e+15 and var_snr < hb_snr_y < 1e+15:
                            if (h_b_y + hb_ab_y > 0 and h_g_y + hg_ab_y > 0) or i in list_add:
                                list_fg_HgHb_value.append(-2.5 * math.log(h_g_y / (h_b_y * 0.47)))
                                list_fg_HgHb_dist_kpc.append(dist_kpc)

                        # if 0.001 < mag606_i < var_mag and 0.001 < mag775_i < var_mag and mag606_err_i < var_err and mag775_err_i < var_err:
                        #     list_F6F7_value.append(mag606_i - mag775_i)
                        #     list_F6F7_dist_kpc.append(dist_kpc)
                        #
                        # if 0.001 < mag606_i < var_mag and 0.001 < mag850_i < var_mag and mag606_err_i < var_err and mag850_err_i < var_err:
                        #     list_F6F8_value.append(mag606_i - mag850_i)
                        #     list_F6F8_dist_kpc.append(dist_kpc)
                        #
                        # if 0.001 < mag606_i < var_mag and 0.001 < mag775_i < var_mag and mag606_err_i < var_err and mag775_err_i < var_err:
                        #     if var_snr < ha_snr_i < 1e+15 and var_snr < hb_snr_i < 1e+15:
                        #         if (h_b_i + hb_ab_i > 0 and h_a_i + ha_ab_i > 0) or i in list_add:
                        #             list_F6F7_value2.append(mag606_i - mag775_i)
                        #             list_HaHb_value2.append(-2.5 * math.log(h_a_i/(h_b_i * 2.85)))
                        #
                        # if 0.001 < mag606_i < var_mag and 0.001 < mag775_i < var_mag and mag606_err_i < var_err and mag775_err_i < var_err:
                        #     if var_snr < hg_snr_i < 1e+15 and var_snr < hb_snr_i < 1e+15:
                        #         if (h_b_i + hb_ab_i > 0 and h_g_i + ha_ab_i > 0) or i in list_add:
                        #             list_F6F7_value3.append(mag606_i - mag775_i)
                        #             list_HgHb_value2.append(h_g_i / (h_b_i * 0.469))

print(count_pairs, count_hahb, count_hghb)

plt.scatter(list_HaHb_dist_kpc, list_HaHb_value, s=2.5)
plt.title("Hα/Hβ-Dist"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("Hα/Hβ/0.47")
plt.savefig(image_path + "HaHb_scatter", dpi=200), plt.close()

plt.scatter(list_HgHb_dist_kpc, list_HgHb_value, s=1.5)
plt.title("Hγ/Hβ-Dist"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("Hγ/Hβ/0.47")
plt.savefig(image_path + "HgHb_scatter", dpi=200), plt.close()

# Binning Ha/Hb
bin_arr = np.linspace(10, var_dist_kpc, int(var_dist_kpc / 10), dtype=int)
bin_lst = bin_arr.tolist()
direc = {}
for ele in bin_lst:
    direc["kpc_%s" % ele] = []

for ele in range(len(list_HaHb_value)):
    direc["kpc_" + str(int(np.ceil(list_HaHb_dist_kpc[ele] / 10)) * 10)].append(list_HaHb_value[ele])

list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
for ele in range(int(var_dist_kpc / 10)):
    list_bin_median.append(np.median(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_mean.append(np.mean(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_std.append(np.std(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_steps.append((ele + 1) * 10)
    list_bin_count.append(len(direc["kpc_" + str((ele + 1) * 10)]))

plt.bar(list_bin_steps, list_bin_count, width=4)
plt.title("Anzahl Galaxien pro Bin, E(Hα-Hβ)"), plt.xlabel("Bins"), plt.ylabel("Anzahl Galaxien pro Bin")
plt.savefig(image_path + "HaHb_bar", dpi=200), plt.close()

plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
plt.scatter(list_bin_steps, list_bin_mean, zorder=2), plt.title("E(Hα-Hβ), Mittelwert"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hα-Hβ)")
plt.savefig(image_path + "HaHb_mean_errorbar", dpi=200), plt.close()

# x_a = np.linspace(min(list_bin_steps), max(list_bin_steps), 20)
# co = np.polyfit(list_bin_steps, list_bin_median, 1)
# linfit = np.poly1d(co)
# plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3)
plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
plt.scatter(list_bin_steps, list_bin_median, zorder=2, color="xkcd:jade green"),
plt.ylim(np.nanmin(list_bin_median) - 0.1, np.nanmax(list_bin_median) + 0.1), plt.title("E(Hα-Hβ), Median"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hα-Hβ)")
plt.savefig(image_path + "HaHb_median", dpi=200), plt.close()

plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
plt.plot(list_bin_steps, list_bin_mean, zorder=2, label="mean", linewidth=1.1), plt.scatter(list_bin_steps, list_bin_median, zorder=3, label="median", color="xkcd:jade green")
plt.title("E(Hα-Hβ), Mittelwert und Median"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hα-Hβ)"), plt.legend()
plt.savefig(image_path + "HaHb_mean-median_errorbar", dpi=200), plt.close()

# co = np.polyfit(list_bin_steps, list_bin_mean, 1)
# linfit = np.poly1d(co)
# plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3)
plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
plt.scatter(list_bin_steps, list_bin_mean, zorder=2)
plt.title("E(Hα-Hβ), Mittelwert"),  plt.ylim(np.nanmin(list_bin_mean)-0.1, np.nanmax(list_bin_mean)+0.1), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hα-Hβ)")
plt.savefig(image_path + "HaHb_mean", dpi=200), plt.close()

# Binning Ha/Hb fg
bin_arr = np.linspace(10, var_dist_kpc, int(var_dist_kpc / 10), dtype=int)
bin_lst = bin_arr.tolist()
direc = {}
for ele in bin_lst:
    direc["kpc_%s" % ele] = []

for ele in range(len(list_fg_HaHb_value)):
    direc["kpc_" + str(int(np.ceil(list_fg_HaHb_dist_kpc[ele] / 10)) * 10)].append(list_fg_HaHb_value[ele])

list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
for ele in range(int(var_dist_kpc / 10)):
    list_bin_median.append(np.median(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_mean.append(np.mean(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_std.append(np.std(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_steps.append((ele + 1) * 10)
    list_bin_count.append(len(direc["kpc_" + str((ele + 1) * 10)]))

plt.bar(list_bin_steps, list_bin_count, width=4)
plt.title("Anzahl Galaxien pro Bin, E(Hα-Hβ), fg"), plt.xlabel("Bins"), plt.ylabel("Anzahl Galaxien pro Bin")
plt.savefig(image_path + "fg_HaHb_bar", dpi=200), plt.close()

plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
plt.scatter(list_bin_steps, list_bin_mean, zorder=2), plt.title("E(Hα-Hβ), Mittelwert, fg"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hα-Hβ)")
plt.savefig(image_path + "fg_HaHb_mean_errorbar", dpi=200), plt.close()

# x_a = np.linspace(min(list_bin_steps), max(list_bin_steps), 20)
# co = np.polyfit(list_bin_steps, list_bin_median, 1)
# linfit = np.poly1d(co)
# plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3)
plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
plt.scatter(list_bin_steps, list_bin_median, zorder=2, color="xkcd:jade green"),
plt.ylim(np.nanmin(list_bin_median) - 0.1, np.nanmax(list_bin_median) + 0.1), plt.title("E(Hα-Hβ), Median, fg"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hα-Hβ)")
plt.savefig(image_path + "fg_HaHb_median", dpi=200), plt.close()

plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
plt.plot(list_bin_steps, list_bin_mean, zorder=2, label="mean", linewidth=1.1), plt.scatter(list_bin_steps, list_bin_median, zorder=3, label="median", color="xkcd:jade green")
plt.title("E(Hα-Hβ), Mittelwert und Median, fg"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hα-Hβ)"), plt.legend()
plt.savefig(image_path + "fg_HaHb_mean-median_errorbar", dpi=200), plt.close()

# co = np.polyfit(list_bin_steps, list_bin_mean, 1)
# linfit = np.poly1d(co)
# plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3)
plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
plt.scatter(list_bin_steps, list_bin_mean, zorder=2)
plt.title("E(Hα-Hβ), Mittelwert, fg"),  plt.ylim(np.nanmin(list_bin_mean)-0.1, np.nanmax(list_bin_mean)+0.1), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hα-Hβ)")
plt.savefig(image_path + "fg_HaHb_mean", dpi=200), plt.close()

# Binning Hg/Hb
bin_arr = np.linspace(10, var_dist_kpc, int(var_dist_kpc / 10), dtype=int)
bin_lst = bin_arr.tolist()
direc = {}
for ele in bin_lst:
    direc["kpc_%s" % ele] = []

for ele in range(len(list_HgHb_value)):
    direc["kpc_" + str(int(np.ceil(list_HgHb_dist_kpc[ele] / 10)) * 10)].append(list_HgHb_value[ele])

list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
for ele in range(int(var_dist_kpc / 10)):
    list_bin_median.append(np.median(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_mean.append(np.mean(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_std.append(np.std(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_steps.append((ele + 1) * 10)
    list_bin_count.append(len(direc["kpc_" + str((ele + 1) * 10)]))

plt.bar(list_bin_steps, list_bin_count, width=4)
plt.title("Anzahl Galaxien pro Bin, E(Hγ-Hβ)"), plt.xlabel("Bins"), plt.ylabel("Anzahl Galaxien pro Bin")
plt.savefig(image_path + "HgHb_bar", dpi=200), plt.close()

plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1), plt.scatter(list_bin_steps, list_bin_mean)
plt.title("E(Hγ-Hβ), Mittelwert"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hγ-Hβ)")
plt.savefig(image_path + "HgHb_mean_errorbar", dpi=200), plt.close()

# x_a = np.linspace(min(list_bin_steps), max(list_bin_steps), 20)
# co = np.polyfit(list_bin_steps, list_bin_median, 1)
# linfit = np.poly1d(co)
# plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3)
plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
plt.scatter(list_bin_steps, list_bin_median, color="xkcd:jade green"),
plt.ylim(np.nanmin(list_bin_median) - 0.1, np.nanmax(list_bin_median) + 0.1), plt.title("E(Hγ-Hβ), Median"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hγ-Hβ)")
plt.savefig(image_path + "HgHb_median", dpi=200), plt.close()

plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
plt.plot(list_bin_steps, list_bin_mean, zorder=2, label="mean", linewidth=1.1), plt.scatter(list_bin_steps, list_bin_median, zorder=3, label="median", color="xkcd:jade green")
plt.title("E(Hγ-Hβ), Mittelwert und Median"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hγ-Hβ)"), plt.legend()
plt.savefig(image_path + "HgHb_mean-median_errorbar", dpi=200), plt.close()

# co = np.polyfit(list_bin_steps, list_bin_mean, 1)
# linfit = np.poly1d(co)
# plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3)
plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
plt.scatter(list_bin_steps, list_bin_mean, zorder=2)
plt.ylim(np.nanmin(list_bin_mean) - 0.1, np.nanmax(list_bin_mean) + 0.1), plt.title("E(Hγ-Hβ), Mittelwert"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hγ-Hβ)")
plt.savefig(image_path + "HgHb_mean", dpi=200), plt.close()

# Binning HgHb fg
bin_arr = np.linspace(10, var_dist_kpc, int(var_dist_kpc / 10), dtype=int)
bin_lst = bin_arr.tolist()
direc = {}
for ele in bin_lst:
    direc["kpc_%s" % ele] = []

for ele in range(len(list_fg_HgHb_value)):
    direc["kpc_" + str(int(np.ceil(list_fg_HgHb_dist_kpc[ele] / 10)) * 10)].append(list_fg_HgHb_value[ele])

list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
for ele in range(int(var_dist_kpc / 10)):
    list_bin_median.append(np.median(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_mean.append(np.mean(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_std.append(np.std(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_steps.append((ele + 1) * 10)
    list_bin_count.append(len(direc["kpc_" + str((ele + 1) * 10)]))

plt.bar(list_bin_steps, list_bin_count, width=4)
plt.title("Anzahl Galaxien pro Bin, E(Hγ-Hβ), fg"), plt.xlabel("Bins"), plt.ylabel("Anzahl Galaxien pro Bin")
plt.savefig(image_path + "fg_HgHb_bar", dpi=200), plt.close()

plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1), plt.scatter(list_bin_steps, list_bin_mean)
plt.title("E(Hγ-Hβ), Mittelwert, fg"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hγ-Hβ)")
plt.savefig(image_path + "fg_HgHb_mean_errorbar", dpi=200), plt.close()

x_a = np.linspace(min(list_bin_steps), max(list_bin_steps), 20)
co = np.polyfit(list_bin_steps, list_bin_median, 1)
linfit = np.poly1d(co)
plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
plt.scatter(list_bin_steps, list_bin_median, color="xkcd:jade green"), plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3)
plt.ylim(np.nanmin(list_bin_median) - 0.1, np.nanmax(list_bin_median) + 0.1), plt.title("E(Hγ-Hβ), Median, fg"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hγ-Hβ)")
plt.savefig(image_path + "fg_HgHb_median", dpi=200), plt.close()

plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
plt.plot(list_bin_steps, list_bin_mean, zorder=2, label="mean", linewidth=1.1), plt.scatter(list_bin_steps, list_bin_median, zorder=3, label="median", color="xkcd:jade green")
plt.title("E(Hγ-Hβ), Mittelwert und Median, fg"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hγ-Hβ)"), plt.legend()
plt.savefig(image_path + "fg_HgHb_mean-median_errorbar", dpi=200), plt.close()

co = np.polyfit(list_bin_steps, list_bin_mean, 1)
linfit = np.poly1d(co)
plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3), plt.scatter(list_bin_steps, list_bin_mean, zorder=2)
plt.ylim(np.nanmin(list_bin_mean) - 0.1, np.nanmax(list_bin_mean) + 0.1), plt.title("E(Hγ-Hβ), Mittelwert, fg"), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hγ-Hβ)")
plt.savefig(image_path + "fg_HgHb_mean", dpi=200), plt.close()
