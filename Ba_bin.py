import os
import numpy as np
import math
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from mpdaf.obj import Cube
import time

# tc-id, 200kpc, 30snr; 285, 193
# tc-id, 300kpc, 30snr; 275
# tc-id, 300kpc, 30snr, 24 mag; 239

# 5snr, HaHb; 296

# vorher tc-id 21

# var_dist = 10
# var_dist_kpc = 250
# # variable, limit min and max redshift for galaxies
# var_z_max = 0.91
# var_z_min = 0.01
# # variable, lower flux limit for at least 1 galaxy
# var_flux = 0
# fac = 0
# var_confy = 0
# var_snr = 25
# var_mag = 24

# variable, upper limit for projected distance in arcsec
var_dist = 1000
var_dist_kpc = 150
# variable, lower flux limit for at least 1 galaxy
var_flux = 0
fac = 0
var_confy = 0
var_snr = 15
var_mag = 1000

var_mode = "ab"
# variable, limit min and max redshift for galaxies
var_z_min = 0.01
if var_mode == "ab":
    var_z_max = 0.41
else:
    var_z_max = 0.91


folder_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\Images_bin\\"
fit_main_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\UDF_HUDF_Catalogs_FULL\\dr2_main_09.fits"
cube_UDF_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\DATACUBE_UDF-10.fits"
cube_MXDF_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\DATACUBE_MXDF.fits"
cube_MOSAIC_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\DATACUBE_UDF-MOSAIC.fits"
region_path_1 = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\objects_region_1.reg"
region_path_2 = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\objects_region_2.reg"


# function kpc/arcsec at redshift z, conv_func(z)
a, b, c, d, e, f, g, h, o, j, k, l, m, n = [-24.862278979480052, 21.267147446426016, 0.005306164849709134, 19.600027472637663, -11.409572942008737, 4.751202039434343, -1.3170007020351893, 0.2021612052937083, -0.007249014401993737, 0.0016127139142351095, -0.00017642802527037128, 9.986804760095554e-06, -2.1958087984210742e-07, -1.148703972276188e-09]


def conv_func(x):
    return n * x**14 + m * x**13 + l * x**12 + k * x**11 + j * x**10 + o * x**9 + h * x**7 + g * x**6 + f * x**5 + e * x**4 + d * x**3 + a * x**2 + b * x + c

def rem_func(loop):
    var_id = 0
    for ele in range(loop):
        outfile1, outfile2 = open(region_path_1, "w"), open(region_path_2, "w")
        outfile1.write('# Region file format: DS9 version 4.1\n')
        outfile1.write('global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        outfile1.write('fk5\n')

        if var_id == 1:
            max_index_sort = list_size.index(max(list_size))
            list_size.pop(max_index_sort)
            list_id.pop(max_index_sort)
        var_id = 1

        max_list = list_size.index(max(list_size))
        max_index = list_id[max_list]

        print(max_index)
        direc2["list_" + str(max_index)].sort()
        print(direc2["list_" + str(max_index)])
        if max_index > 2:
            print(len(direc2["list_" + str(max_index - 2)]), len(direc2["list_" + str(max_index - 1)]),
                  len(direc2["list_" + str(max_index)]), len(direc2["list_" + str(max_index + 1)]),
                  len(direc2["list_" + str(max_index + 2)]))

        rfg, dfg, zfg = ra[max_index - 1], dec[max_index - 1], z[max_index - 1]
        if ele >= 0:
            outfile1.write('circle(%s,%s,0.9")# color=red width=2 font="helvetica 14 normal roman"\n' % (rfg, dfg))
        print(zfg)

        list_hghb, list_hahb, list_dist_ab, list_dist_gb = [], [], [], []
        for ele3 in direc2["list_" + str(max_index)]:
            rt, dt, zt = ra[ele3], dec[ele3], z[ele3]
            h_a_bg, h_b_bg, h_g_bg = h_a[ele3], h_b[ele3], h_g[ele3]
            ha_snr_bg, hb_snr_bg, hg_snr_bg = h_a_snr[ele3], h_b_snr[ele3], h_g_snr[ele3]
            if ele >= 0:
                outfile1.write('circle(%s,%s,0.8")# color=green width=2 font="helvetica 14 normal roman"\n' % (rt, dt))
            dist = ((rfg - rt) ** 2 + (dfg - dt) ** 2) ** (1 / 2)
            dist_kpc = dist * 3600 * conv_func(zfg)
            if hg_snr_bg > var_snr and hb_snr_bg > var_snr and (dist * 3600 * conv_func(zfg)) >= 0:
                # list_hghb.append(h_g_bg / (h_b_bg * 0.47))
                list_hghb.append(-2.5 * math.log(h_g_bg / (h_b_bg * 0.47)))
                list_dist_gb.append(dist_kpc)
            if ha_snr_bg > var_snr and hb_snr_bg > var_snr and (dist * 3600 * conv_func(zfg)) >= 0:
                # list_hahb.append(h_a_bg / (h_b_bg * 2.86))
                list_hahb.append((-2.5 * math.log(h_a_bg / (h_b_bg * 2.85))))

                list_dist_ab.append(dist_kpc)

        # Balmer-Decrement
        if var_mode == "ab":
            plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
            # plt.scatter(list_dist_ab, list_hahb, s=14), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("Hα/Hβ/2.86")
            plt.scatter(list_dist_ab, list_hahb, s=14), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hα - Hβ)"), plt.title("E(Hα - Hβ), bin")
            plt.savefig(folder_path + str(max_index) + "_hahb.png")
            plt.close()

            bin_arr = np.linspace(10, var_dist_kpc, int(var_dist_kpc / 10), dtype=int)
            bin_lst = bin_arr.tolist()
            direc = {}
            for ele in bin_lst:
                direc["kpc_%s" % ele] = []
            for ele in range(len(list_hahb)):
                direc["kpc_" + str(int(np.ceil(list_dist_ab[ele] / 10)) * 10)].append(list_hahb[ele])

            list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
            for ele in range(int(var_dist_kpc / 10)):
                list_bin_median.append(np.median(direc["kpc_" + str((ele + 1) * 10)]))
                list_bin_mean.append(np.mean(direc["kpc_" + str((ele + 1) * 10)]))
                list_bin_std.append(np.std(direc["kpc_" + str((ele + 1) * 10)]))
                list_bin_steps.append((ele + 1) * 10)
                list_bin_count.append(len(direc["kpc_" + str((ele + 1) * 10)]))

            # plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
            plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
            plt.scatter(list_bin_steps, list_bin_mean, zorder=1, label="Mean", s=32, edgecolors="none"), plt.scatter(list_bin_steps, list_bin_median, color="xkcd:jade green", label="Median", s=40, edgecolors="xkcd:jade green", facecolor="none")
            plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hα - Hβ)"), plt.legend(loc='lower right'), plt.ylim(np.nanmin(list_bin_mean)-0.1, np.nanmax(list_bin_mean)+0.1), plt.title("E(Hα - Hβ), Mittelwert und Median")
            plt.savefig(folder_path + str(max_index) + "_hahb_mean-median.png", dpi=150)
            plt.close()

        if var_mode == "gb":
            plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
            plt.scatter(list_dist_gb, list_hghb, s=14), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hγ - Hβ)"), plt.title("E(Hγ - Hβ), bin")
            plt.savefig(folder_path + str(max_index) + "_hghb.png")
            plt.close()

            bin_arr = np.linspace(10, var_dist_kpc, int(var_dist_kpc / 10), dtype=int)
            bin_lst = bin_arr.tolist()
            direc = {}
            for ele in bin_lst:
                direc["kpc_%s" % ele] = []
            for ele in range(len(list_hghb)):
                direc["kpc_" + str(int(np.ceil(list_dist_gb[ele] / 10)) * 10)].append(list_hghb[ele])

            list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
            for ele in range(int(var_dist_kpc / 10)):
                list_bin_median.append(np.median(direc["kpc_" + str((ele + 1) * 10)]))
                list_bin_mean.append(np.mean(direc["kpc_" + str((ele + 1) * 10)]))
                list_bin_std.append(np.std(direc["kpc_" + str((ele + 1) * 10)]))
                list_bin_steps.append((ele + 1) * 10)
                list_bin_count.append(len(direc["kpc_" + str((ele + 1) * 10)]))

            # plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
            plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
            plt.scatter(list_bin_steps, list_bin_mean, zorder=1, label="Mean", s=32, edgecolors="none"), plt.scatter(list_bin_steps, list_bin_median, color="xkcd:jade green", label="Median", s=40, edgecolors="xkcd:jade green", facecolor="none")
            plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("E(Hγ - Hβ)"), plt.legend(loc='lower right'), plt.ylim(np.nanmin(list_bin_mean)-0.1, np.nanmax(list_bin_mean)+0.1), plt.title("E(Hγ - Hβ), Mittelwert und Median")
            plt.savefig(folder_path + str(max_index) + "_hghb_mean-median.png", dpi=150)
            plt.close()

cube_UDF, cube_MXDF, cube_MOSAIC = Cube(cube_UDF_path, ext=1), Cube(cube_MXDF_path, ext=1), Cube(cube_MOSAIC_path, ext=1)
outfile1, outfile2 = open(region_path_1, "w"), open(region_path_2, "w")
fits_info = fits.open(fit_main_path)

outfile1.write('# Region file format: DS9 version 4.1\n')
outfile1.write('global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
outfile1.write('fk5\n')
outfile2.write('# Region file format: DS9 version 4.1\n')
outfile2.write('global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
outfile2.write('fk5\n')

# right ascension
ra = fits_info[1].data['RA']
# declination
dec = fits_info[1].data['DEC']
# red shift
z = fits_info[1].data['Z']
# redshift confidence
zconf = fits_info[1].data["ZCONF"]
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


# variable, count of galaxy pairs
count, count_HaHb, count_HgHb, count_BV = 0, 0, 0, 0
list_pairs, list_dist = ["i", "y", "dist"], []

bin_arr = np.linspace(1, 2300, 2300, dtype=int)
bin_lst = bin_arr.tolist()
direc2 = {}
for ele in bin_lst:
    direc2["list_%s" % ele] = []
# goes through list table, saves right ascension and declination for every object in var r, d
for i in range(len(ra)):
    ri, di, zi, zconfi, flux_i, mag606_i = ra[i], dec[i], float(z[i]), zconf[i], flux[i], mag606[i]
    print(ri, di, zi, i)

    for y in range(len(ra)):
        # goes through list table, saves right ascension and declination, flux, emission lines for
        ry, dy, zy, zconfy, mag606_y, hay, hby, hgy = ra[y], dec[y], float(z[y]), zconf[y], mag606[y], h_a[y], h_b[y], h_g[y]
        ha_snry, hb_snry, hg_snry = h_a_snr[y], h_b_snr[y], h_g_snr[y]
        ha_ab_y, hb_ab_y, hg_ab_y = ha_ab[y], hb_ab[y], hg_ab[y]

        if i != y:
            # prevents the pairing of galaxies with themselves
            if zi < zy and mag606_i < var_mag:
                if var_z_min < zi < var_z_max and var_z_min < zy < var_z_max:
                    # distance to foreground galaxy
                    dist = ((ry - ri) ** 2 + (dy - di) ** 2) ** (1 / 2)
                    dist_kpc = dist * 3600 * conv_func(zi)

                    if dist_kpc <= var_dist_kpc and dist * 3600 < var_dist:
                        if var_mode == "ab":
                          if var_snr < ha_snry < 1e+15 and var_snr < hb_snry < 1e+15:
                              if hay + ha_ab_y > 0 and hby + hb_ab_y > 0:
                                direc2["list_"+str(i+1)].append(y)
                                print("Pair: "+str(count), "i: " + str(i), "y: " + str(y), "dist-kpc", dist_kpc, dist)
                        if var_mode == "gb":
                          if var_snr < hg_snry < 1e+15 and var_snr < hb_snry < 1e+15:
                              if hgy + hg_ab_y > 0 and hby + hb_ab_y > 0:
                                direc2["list_"+str(i+1)].append(y)
                                print("Pair: "+str(count), "i: " + str(i), "y: " + str(y), "dist-kpc", dist_kpc, dist)

list_id, list_size, list_z = [], [], []
for ele in range(2300):
    list_id.append(ele+1)
    list_size.append(len(direc2["list_"+str(ele+1)]))

rem_func(5)


# close files
outfile1.close()
outfile2.close()
