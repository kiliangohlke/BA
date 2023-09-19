import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpdaf.obj import Cube
import math

# var_dist_kpc_min = 0 var_dist_kpc = 250 var_z_max = 12 var_snr = 15 var_err = 0.5 var_mag = 25
# remove outlier, var_err, test mass

# trend mag var_mag = 27, 28

var_dist_kpc_min = 0
var_dist_kpc = 300
var_z_max = 0.91
var_snr = 5
var_err = 100
var_mag = 100

image_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\Images_Rötung\\"
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
# mass in log(M/M_sun)
mass = fits_info[1].data["MASS_PRO"]
# Dataset of galaxy
dat = fits_info[1].data["DATASET"]

count_pairs, count_hahb, count_hghb = 0, 0, 0

list_HaHb_value, list_HaHb_dist_kpc, list_HgHb_value, list_HgHb_dist_kpc = [], [], [], []
list_fg_HaHb_value, list_fg_HaHb_dist_kpc, list_fg_HgHb_value, list_fg_HgHb_dist_kpc = [], [], [], []
list_add = [14, 202, 207, 417]

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

# HaHb bg
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

# x_a = np.linspace(np.nanmin(list_bin_steps), np.nanmax(list_bin_steps), 20)
# co = np.polyfit(list_bin_steps, list_bin_median, 1)
# linfit = np.poly1d(co)

fig, axs = plt.subplots(2, 1, sharex="col", sharey="col", dpi=200)
fig.subplots_adjust(hspace=0)
fig.text(0.5, 0.03, "Entfernung in Vordergrundebene [kpc]", ha="center", va="center")
fig.text(0.025, 0.5, "E(Hα-Hβ)", ha="center", va="center", rotation="vertical")
axs[0].set_title("E(Hα-Hβ), fg-bg, Median")
plot1 = axs[0].scatter(list_bin_steps, list_bin_median, s=16, zorder=3, linewidth=0.6, label="median bg galaxies", color="xkcd:jade green")
axs[0].spines.bottom.set_visible(False)
axs[0].grid(linestyle="-.", linewidth=0.15, zorder=1)
axs[0].set_ylim(np.nanmin(list_bin_median) - 0.05, np.nanmax(list_bin_median) + 0.05)
axs[0].set_ylim(np.nanmin(list_bin_median) - 0.03, np.nanmax(list_bin_median) + 0.03)
axs[0].axhline(y=np.nanmin(list_bin_median) - 0.03, color="black", linewidth=1.25)
axs[0].hlines(-0.2, -100, 400, linestyle="-.", linewidth=0.4, color="xkcd:stone")
axs[0].set_xlim(-20, 320)

# HaHb fg
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

list_bin_steps2, list_bin_median2 = list_bin_steps.copy(), list_bin_median.copy()
list_bin_steps2.pop(0), list_bin_median2.pop(0)

# x_a = np.linspace(np.nanmin(list_bin_steps), np.nanmax(list_bin_steps), 20)
# co = np.polyfit(list_bin_steps2, list_bin_median2, 1)
# linfit = np.poly1d(co)

plot2 = axs[1].scatter(list_bin_steps, list_bin_median, s=16, zorder=3, linewidth=0.6, label="median fg galaxies", color="indianred")
axs[1].grid(linestyle="-.", linewidth=0.15, zorder=1)
axs[1].spines.top.set_visible(False)
axs[1].legend((plot1, plot2), ("median bg galaxies", "median fg galaxies"), loc="lower right")
axs[1].hlines(-0.2, -100, 400, linestyle="-.", linewidth=0.4, color="xkcd:stone")
axs[1].set_xlim(-20, 320)

plt.show()

# HgHb bg
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

# x_a = np.linspace(np.nanmin(list_bin_steps), np.nanmax(list_bin_steps), 20)
# co = np.polyfit(list_bin_steps, list_bin_median, 1)
# linfit = np.poly1d(co)

fig, axs = plt.subplots(2, 1, sharex="col", sharey="col", dpi=200)
fig.subplots_adjust(hspace=0)
axs[0].set_title("E(Hγ-Hβ), fg-bg, Median")
fig.text(0.5, 0.03, "Entfernung in Vordergrundebene [kpc]", ha="center", va="center")
fig.text(0.03, 0.5, "E(Hγ-Hβ)", ha="center", va="center", rotation="vertical")
plot1 = axs[0].scatter(list_bin_steps, list_bin_median, s=16, zorder=3, linewidth=0.6, label="median bg galaxies", color="xkcd:jade green")
axs[0].grid(linestyle="-.", linewidth=0.15, zorder=1)
axs[0].spines.bottom.set_visible(False)
axs[0].set_ylim(np.nanmin(list_bin_median) - 0.04, np.nanmax(list_bin_median) + 0.01)
axs[0].axhline(y=np.nanmin(list_bin_median) - 0.04, color="black", linewidth=1.25)
axs[0].hlines(0.1, -100, 400, linestyle="-.", linewidth=0.4, color="xkcd:stone")
axs[0].set_xlim(-20, 320)

# HgHb fg

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

list_bin_steps2, list_bin_median2 = list_bin_steps.copy(), list_bin_median.copy()
list_bin_steps2.pop(0), list_bin_median2.pop(0)

# x_a = np.linspace(np.nanmin(list_bin_steps), np.nanmax(list_bin_steps), 20)
# co = np.polyfit(list_bin_steps2, list_bin_median2, 1)
# linfit = np.poly1d(co)
# plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
# plt.plot(x_a, linfit(x_a), color="xkcd:faded green", linewidth=0.9, zorder=3),

plot2 = axs[1].scatter(list_bin_steps, list_bin_median, s=16, zorder=3, linewidth=0.6, label="median fg galaxies", color="indianred")
axs[1].grid(linestyle="-.", linewidth=0.15, zorder=1)
axs[1].spines.top.set_visible(False)
axs[0].legend((plot1, plot2), ("median bg galaxies", "median fg galaxies"), loc="upper right")
axs[1].hlines(0.1, -100, 400, linestyle="-.", linewidth=0.4, color="xkcd:stone")
axs[1].set_xlim(-20, 320)
plt.show()






