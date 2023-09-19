import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from mpdaf.obj import Cube
import time
# var_dist_kpc_min = 0 var_dist_kpc = 250 var_z_max = 12 var_snr = 15 var_err = 0.5 var_mag = 25
# remove outlier, var_err, test mass

# trend mag var_mag = 27, 28

var_mode = "single"
var_dist_kpc_min = 0
var_dist_kpc = 300
var_z_max = 12
var_err = 0.3
var_mag = 100
var_mass = 0


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

count_F6F7 = 0
list_F6F7_value, list_F6F7_dist_kpc, list_F7F8_value, list_F7F8_dist_kpc = [], [], [], []
list_fg_F6F7_value, list_fg_F6F7_dist_kpc, list_fg_F7F8_value, list_fg_F7F8_dist_kpc = [], [], [], []
list_bv, list_gal_dist = [], []
list_i = [300, 600, 900, 1200, 1500, 1800]

if var_mode == "single":
    list_HaHb_value, list_HaHb_dist_kpc, list_HgHb_value, list_HgHb_dist_kpc = [], [], [], []
    list_F6F7_value, list_F6F7_dist_kpc, list_F7F8_value, list_F7F8_dist_kpc = [], [], [], []
    # goes through list table, saves right ascension and declination for every object in var r, d
    for i in range(len(ra)):
        ri, di, zi, flux_i = ra[i], dec[i], float(z[i]), flux[i]
        if i in list_i:
            input("Enter to continue")
        for y in range(len(ra)):
            # goes through list table, saves right ascension and declination, flux, emission lines for
            ry, dy, zy, mag606_y, mass_y = ra[y], dec[y], float(z[y]), mag606[y], mass[y]

            if i != y:
                # if zi > zy and 8 < mass_y < 1e+10:
                if zi > zy and mag606_y < var_mag:
                # if zi > zy and mag606_y < var_mag and var_mass < mass_y < 1e+10:
                    if 0.01 < zi < var_z_max and zy > 0.01 and 0.01 < zi - zy < 1:
                        dist = ((ry - ri) ** 2 + (dy - di) ** 2) ** (1 / 2)
                        dist_kpc = dist * 3600 * conv_func(zy)

                        if var_dist_kpc_min < dist_kpc <= var_dist_kpc:
                            print("")
                            print("i: " + str(i+1), "y: " + str(y+1), "dist-kpc", dist_kpc, dist)

                            mag606_i, mag775_i, mag850_i, mag606_err_i, mag775_err_i, mag850_err_i = mag606[i], mag775[i], mag850[i], mag606_err[i], mag775_err[i], mag850_err[i]
                            mag606_y, mag775_y, mag850_y, mag606_err_y, mag775_err_y, mag850_err_y = mag606[y], mag775[y], mag850[y], mag606_err[y], mag775_err[y], mag850_err[y]


                            # Reddening
                            if mag606_err_i < var_err and mag775_err_i < var_err:
                                count_F6F7 = count_F6F7 + 1
                                list_F6F7_value.append(mag606_i - mag775_i)
                                list_F6F7_dist_kpc.append(dist_kpc)

                            if mag606_err_y < var_err and mag775_err_y < var_err:
                                list_fg_F6F7_value.append(mag606_y - mag775_y)
                                list_fg_F6F7_dist_kpc.append(dist_kpc)

                            if mag775_err_i < var_err and mag850_err_i < var_err:
                                list_F7F8_value.append(mag775_i - mag850_i)
                                list_F7F8_dist_kpc.append(dist_kpc)

                            sp_i = 0
                            if str(dat[i]) == "UDF10":
                                sp_i = cube_UDF.aperture((di, ri), 0.5)
                            if str(dat[i]) == "MXDF":
                                sp_i = cube_MXDF.aperture((di, ri), 0.5)
                            if str(dat[i]) == "MOSAIC":
                                sp_i = cube_MOSAIC.aperture((di, ri), 0.5)

                            if abs(np.std(sp_i.data)) < abs(np.median(sp_i.data)):
                                print(sp_i.abmag_filter_name("B")[0], sp_i.abmag_filter_name("V")[0], sp_i.abmag_filter_name("Rc")[0], sp_i.abmag_filter_name("Ic")[0])
                                list_bv.append(sp_i.abmag_filter_name("B")[0] - sp_i.abmag_filter_name("V")[0])
                                # list_v.append(sp_i.abmag_filter_name("V")[0])
                                # list_r.append(sp_i.abmag_filter_name("Rc")[0])
                                # list_i.append(sp_i.abmag_filter_name("Ic")[0])
                                list_gal_dist.append(dist * 3600 * conv_func(zy))


    print(count_F6F7)
    # 606 - 775
    plt.scatter(list_F6F7_dist_kpc, list_F6F7_value, s=2.5)
    plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F606W - F775W [mag]")
    plt.savefig(image_path + "606-775_scatter", dpi=200), plt.close()

    bin_arr = np.linspace(10, var_dist_kpc, int(var_dist_kpc / 10), dtype=int)
    bin_lst = bin_arr.tolist()
    direc = {}
    for ele in bin_lst:
        direc["kpc_%s" % ele] = []

    for ele in range(len(list_F6F7_value)):
        direc["kpc_" + str(int(np.ceil(list_F6F7_dist_kpc[ele] / 10)) * 10)].append(list_F6F7_value[ele])
    list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
    for ele in range(int(var_dist_kpc / 10)):
        list_bin_median.append(np.median(direc["kpc_" + str((ele + 1) * 10)]))
        list_bin_mean.append(np.mean(direc["kpc_" + str((ele + 1) * 10)]))
        list_bin_std.append(np.std(direc["kpc_" + str((ele + 1) * 10)]))
        list_bin_steps.append((ele + 1) * 10)
        list_bin_count.append(len(direc["kpc_" + str((ele + 1) * 10)]))

    plt.bar(list_bin_steps, list_bin_count, width=4)
    plt.title("Anzahl Galaxien pro Bin, F606W - F775W"), plt.xlabel("Bins"), plt.ylabel("Anzahl Galaxien pro Bin")
    plt.savefig(image_path + "606-775_bar", dpi=200), plt.close()

    plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1), plt.scatter(list_bin_steps, list_bin_mean, zorder=2)
    plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F606W - F775W [mag]")
    plt.savefig(image_path + "606-775_mean_errorbar", dpi=200), plt.close()

    x_a = np.linspace(min(list_bin_steps), np.nanmax(list_bin_steps), 20)
    co = np.polyfit(list_bin_steps, list_bin_median, 1)
    linfit = np.poly1d(co)
    plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
    plt.scatter(list_bin_steps, list_bin_median, zorder=2, color="xkcd:jade green")
    # plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3), plt.xlabel("Entfernung in Vordergrundebene [kpc]"),\
    plt.ylabel("F606W - F775W [mag]"), plt.ylim(np.nanmin(list_bin_median)-0.075, np.nanmax(list_bin_median)+0.075)
    plt.savefig(image_path + "606-775_median", dpi=200), plt.close()

    plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
    plt.plot(list_bin_steps, list_bin_mean, zorder=2, label="mean", linewidth=1.1), plt.scatter(list_bin_steps, list_bin_median, zorder=3, label="median", color="xkcd:jade green")
    plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F606W - F775W [mag]"), plt.legend()
    plt.savefig(image_path + "606-775_mean-median_errorbar", dpi=200), plt.close()

    co = np.polyfit(list_bin_steps, list_bin_mean, 1)
    linfit = np.poly1d(co)
    plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
    # plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3),
    plt.scatter(list_bin_steps, list_bin_mean, zorder=2)
    plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F606W - F775W [mag]"), plt.ylim(np.nanmin(list_bin_mean)-0.075, np.nanmax(list_bin_mean)+0.075)
    plt.savefig(image_path + "606-775_mean", dpi=200), plt.close()

    # # 606 - 775 fg
    # plt.scatter(list_fg_F6F7_dist_kpc, list_fg_F6F7_value, s=2.5)
    # plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F606W - F775W [mag]"), plt.title("fg")
    # plt.savefig(image_path + "fg_606-775_scatter", dpi=200), plt.close()
    #
    # bin_arr = np.linspace(10, var_dist_kpc, int(var_dist_kpc / 10), dtype=int)
    # bin_lst = bin_arr.tolist()
    # direc = {}
    # for ele in bin_lst:
    #     direc["kpc_%s" % ele] = []
    #
    # for ele in range(len(list_fg_F6F7_value)):
    #     direc["kpc_" + str(int(np.ceil(list_fg_F6F7_dist_kpc[ele] / 10)) * 10)].append(list_fg_F6F7_value[ele])
    # list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
    # for ele in range(int(var_dist_kpc / 10)):
    #     list_bin_median.append(np.median(direc["kpc_" + str((ele + 1) * 10)]))
    #     list_bin_mean.append(np.mean(direc["kpc_" + str((ele + 1) * 10)]))
    #     list_bin_std.append(np.std(direc["kpc_" + str((ele + 1) * 10)]))
    #     list_bin_steps.append((ele + 1) * 10)
    #     list_bin_count.append(len(direc["kpc_" + str((ele + 1) * 10)]))
    #
    # plt.bar(list_bin_steps, list_bin_count, width=4)
    # plt.title("Anzahl Galaxien pro Bin, F606W - F775W"), plt.xlabel("Bins"), plt.ylabel("Anzahl Galaxien pro Bin"), plt.title("fg")
    # plt.savefig(image_path + "fg_606-775_bar", dpi=200), plt.close()
    #
    # plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1), plt.scatter(list_bin_steps, list_bin_mean, zorder=2)
    # plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F606W - F775W [mag]")
    # plt.savefig(image_path + "fg_606-775_mean_errorbar", dpi=200), plt.close()
    #
    # x_a = np.linspace(min(list_bin_steps), np.nanmax(list_bin_steps), 20)
    # co = np.polyfit(list_bin_steps, list_bin_median, 1)
    # linfit = np.poly1d(co)
    # plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
    # plt.scatter(list_bin_steps, list_bin_median, zorder=2, color="xkcd:jade green")
    # plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3), plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F606W - F775W [mag]"), plt.ylim(np.nanmin(list_bin_median) - 0.075, np.nanmax(list_bin_median) + 0.075)
    # plt.savefig(image_path + "fg_606-775_median", dpi=200), plt.close()
    #
    # plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
    # plt.plot(list_bin_steps, list_bin_mean, zorder=2, label="mean", linewidth=1.1), plt.scatter(list_bin_steps, list_bin_median, zorder=3, label="median", color="xkcd:jade green")
    # plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F606W - F775W [mag]"), plt.legend(), plt.title("fg")
    # plt.savefig(image_path + "fg_606-775_mean-median_errorbar", dpi=200), plt.close()
    #
    # co = np.polyfit(list_bin_steps, list_bin_mean, 1)
    # linfit = np.poly1d(co)
    # plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
    # plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3), plt.scatter(list_bin_steps, list_bin_mean, zorder=2)
    # plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F606W - F775W [mag]"), plt.ylim(np.nanmin(list_bin_mean) - 0.075, np.nanmax(list_bin_mean) + 0.075), plt.title("fg")
    # plt.savefig(image_path + "fg_606-775_mean", dpi=200), plt.close()
    #
    # 775 - 850
    # plt.scatter(list_F7F8_dist_kpc, list_F7F8_value, s=2.5)
    # plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F775W - F850LP [mag]")
    # plt.savefig(image_path + "775-850_scatter", dpi=200), plt.close()
    #
    # bin_arr = np.linspace(10, var_dist_kpc, int(var_dist_kpc / 10), dtype=int)
    # bin_lst = bin_arr.tolist()
    # direc = {}
    # for ele in bin_lst:
    #     direc["kpc_%s" % ele] = []
    #
    # for ele in range(len(list_F7F8_value)):
    #     direc["kpc_" + str(int(np.ceil(list_F7F8_dist_kpc[ele] / 10)) * 10)].append(list_F7F8_value[ele])
    #
    # list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
    # for ele in range(int(var_dist_kpc / 10)):
    #     list_bin_median.append(np.median(direc["kpc_" + str((ele + 1) * 10)]))
    #     list_bin_mean.append(np.mean(direc["kpc_" + str((ele + 1) * 10)]))
    #     list_bin_std.append(np.std(direc["kpc_" + str((ele + 1) * 10)]))
    #     list_bin_steps.append((ele + 1) * 10)
    #     list_bin_count.append(len(direc["kpc_" + str((ele + 1) * 10)]))
    #
    # plt.bar(list_bin_steps, list_bin_count, width=4)
    # plt.title("Anzahl Galaxien pro Bin, F775W - F850LP"), plt.xlabel("Bins"), plt.ylabel("Anzahl Galaxien pro Bin")
    # plt.savefig(image_path + "775-850_bar", dpi=200), plt.close()
    #
    # plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1), plt.scatter(list_bin_steps, list_bin_mean, zorder=2),
    # plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F775W - F850LP [mag]")
    # plt.savefig(image_path + "775-850_mean_errorbar", dpi=200), plt.close()
    #
    # x_a = np.linspace(np.nanmin(list_bin_steps), np.nanmax(list_bin_steps), 20)
    # co = np.polyfit(list_bin_steps, list_bin_median, 1)
    # linfit = np.poly1d(co)
    # plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
    # plt.scatter(list_bin_steps, list_bin_median, zorder=2, color="xkcd:jade green"), plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3)
    # plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F775W - F850LP [mag]"), plt.ylim(np.nanmin(list_bin_median)-0.075, np.nanmax(list_bin_median)+0.075)
    # plt.savefig(image_path + "775-850_median", dpi=200), plt.close()
    #
    # plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
    # plt.plot(list_bin_steps, list_bin_mean, zorder=2, label="mean", linewidth=1.1), plt.scatter(list_bin_steps, list_bin_median, zorder=3, label="median", color="xkcd:jade green")
    # plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F775W - F850LP [mag]"), plt.legend()
    # plt.savefig(image_path + "775-850_mean-median_errorbar", dpi=200), plt.close()
    #
    # co = np.polyfit(list_bin_steps, list_bin_mean, 1)
    # linfit = np.poly1d(co)
    # plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
    # plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3), plt.scatter(list_bin_steps, list_bin_mean, zorder=2)
    # plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F775W - F850LP [mag]"), plt.ylim(np.nanmin(list_bin_mean)-0.075, np.nanmax(list_bin_mean)+0.075)
    # plt.savefig(image_path + "775-850_mean", dpi=200), plt.close()

    # B-V
    plt.scatter(list_gal_dist, list_bv, s=2.5)
    plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("B - V [mag]")
    plt.savefig(image_path + "B-V", dpi=200), plt.close()

    bin_arr = np.linspace(10, var_dist_kpc, int(var_dist_kpc / 10), dtype=int)
    bin_lst = bin_arr.tolist()
    direc = {}
    for ele in bin_lst:
        direc["kpc_%s" % ele] = []

    for ele in range(len(list_bv)):
        direc["kpc_" + str(int(np.ceil(list_gal_dist[ele] / 10)) * 10)].append(list_bv[ele])
    list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
    for ele in range(int(var_dist_kpc / 10)):
        list_bin_median.append(np.median(direc["kpc_" + str((ele + 1) * 10)]))
        list_bin_mean.append(np.mean(direc["kpc_" + str((ele + 1) * 10)]))
        list_bin_std.append(np.std(direc["kpc_" + str((ele + 1) * 10)]))
        list_bin_steps.append((ele + 1) * 10)
        list_bin_count.append(len(direc["kpc_" + str((ele + 1) * 10)]))

    plt.bar(list_bin_steps, list_bin_count, width=4)
    plt.title("Anzahl Galaxien pro Bin, B - V"), plt.xlabel("Bins"), plt.ylabel("Anzahl Galaxien pro Bin")
    plt.savefig(image_path + "B-V_bar", dpi=200), plt.close()

    plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83,
                 capsize=4, zorder=1), plt.scatter(list_bin_steps, list_bin_mean, zorder=2)
    plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("B - V [mag]")
    plt.savefig(image_path + "B-V_mean_errorbar", dpi=200), plt.close()

    x_a = np.linspace(min(list_bin_steps), np.nanmax(list_bin_steps), 20)
    co = np.polyfit(list_bin_steps, list_bin_median, 1)
    linfit = np.poly1d(co)
    plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
    plt.scatter(list_bin_steps, list_bin_median, zorder=2, color="xkcd:jade green")
    # plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3)
    plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("B - V [mag]"), plt.ylim(np.nanmin(list_bin_median) - 0.075, np.nanmax(list_bin_median) + 0.075)
    plt.savefig(image_path + "B-V_median", dpi=200), plt.close()

    plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
    plt.plot(list_bin_steps, list_bin_mean, zorder=2, label="mean", linewidth=1.1), plt.scatter(list_bin_steps,list_bin_median, zorder=3, label="median", color="xkcd:jade green")
    plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("B - V [mag]"), plt.legend()
    plt.savefig(image_path + "B-V_mean-median_errorbar", dpi=200), plt.close()

    co = np.polyfit(list_bin_steps, list_bin_mean, 1)
    linfit = np.poly1d(co)
    plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
    # plt.plot(x_a, linfit(x_a), color="indianred", linewidth=0.94, zorder=3),
    plt.scatter(list_bin_steps, list_bin_mean, zorder=2)
    plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("B - V [mag]"), plt.ylim(np.nanmin(list_bin_mean) - 0.075, np.nanmax(list_bin_mean) + 0.075)
    plt.savefig(image_path + "B-V_mean", dpi=200), plt.close()


list_mag, list_mass, list_color, list_lin_color = [26, 28, 30, 32], [7.8, 8.5, 9.2, 10], ["xkcd:earth", "xkcd:dusty purple", "tab:red", "tab:orange"], ["xkcd:faded blue", "xkcd:faded green", "xkcd:faded red", "xkcd:faded orange"],
list_mag_range, list_mass_range = [24, 28, 28, 32], [8.7, 12, 6, 8.7]
count_ele = 0

if var_mode == "multi":
    # for ele4 in list_mass:
    # for ele4 in list_mag:
    for ele5 in range(int(len(list_mag_range)/2)):
        ele_range = ele5 * 2
        count_ele = count_ele + 1

        list_HaHb_value, list_HaHb_dist_kpc, list_HgHb_value, list_HgHb_dist_kpc = [], [], [], []
        list_F6F7_value, list_F6F7_dist_kpc, list_F6F8_value, list_F6F8_dist_kpc = [], [], [], []
        # goes through list table, saves right ascension and declination for every object in var r, d
        for i in range(len(ra)):
            ri, di, zi, flux_i = ra[i], dec[i], float(z[i]), flux[i]
            for y in range(len(ra)):
                # goes through list table, saves right ascension and declination, flux, emission lines for
                ry, dy, zy, mag606_y, mass_y = ra[y], dec[y], float(z[y]), mag606[y], mass[y]

                if i != y:
                    # if zi > zy and ele4 < mass_y < 1e+10:
                    if zi > zy and list_mag_range[ele_range] < mag606_y < list_mag_range[ele_range + 1]:
                    # if zi > zy and list_mass_range[ele_range] < mass_y < list_mass_range[ele_range + 1]:
                        if 0.01 < zi < var_z_max and zy > 0.01 and 0.01 < zi - zy < 1:
                            dist = ((ry - ri) ** 2 + (dy - di) ** 2) ** (1 / 2)
                            dist_kpc = dist * 3600 * conv_func(zy)

                            if var_dist_kpc_min < dist_kpc <= var_dist_kpc:
                                print("")
                                print("i: " + str(i+1), "y: " + str(y+1), "dist-kpc", dist_kpc, dist)

                                h_a_i, h_b_i, h_g_i = h_a[i], h_b[i], h_g[i],
                                ha_snr_i, hb_snr_i, hg_snr_i = h_a_snr[i], h_b_snr[i], h_g_snr[i]
                                mag606_i, mag775_i, mag850_i, mag606_err_i, mag775_err_i, mag850_err_i = mag606[i], mag775[i], mag850[i], mag606_err[i], mag775_err[i], mag850_err[i]
                                ha_ab_i, hb_ab_i, hg_ab_i = ha_ab[i], hb_ab[i], hg_ab[i]

                                if mag606_err_i < var_err and mag775_err_i < var_err:
                                    list_F6F7_value.append(mag606_i - mag775_i)
                                    list_F6F7_dist_kpc.append(dist_kpc)


        bin_arr = np.linspace(10, var_dist_kpc, int(var_dist_kpc / 10), dtype=int)
        bin_lst = bin_arr.tolist()
        direc = {}
        for ele in bin_lst:
            direc["kpc_%s" % ele] = []

        for ele in range(len(list_F6F7_value)):
            direc["kpc_" + str(int(np.ceil(list_F6F7_dist_kpc[ele] / 10)) * 10)].append(list_F6F7_value[ele])
        list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
        for ele in range(int(var_dist_kpc / 10)):
            list_bin_median.append(np.median(direc["kpc_" + str((ele + 1) * 10)]))
            list_bin_mean.append(np.mean(direc["kpc_" + str((ele + 1) * 10)]))
            list_bin_std.append(np.std(direc["kpc_" + str((ele + 1) * 10)]))
            list_bin_steps.append((ele + 1) * 10)
            list_bin_count.append(len(direc["kpc_" + str((ele + 1) * 10)]))

        if count_ele-1 == 0:
            fig, axs = plt.subplots(2, 1, sharex="col", sharey="col", dpi=200)
            fig.subplots_adjust(hspace=0)
            fig.text(0.5, 0.03, "Entfernung in Vordergrundebene [kpc]", ha="center", va="center")
            fig.text(0.03, 0.5, "F606W - F775W [mag]", ha="center", va="center", rotation="vertical")
            axs[count_ele - 1].spines.bottom.set_visible(False)
            axs[count_ele - 1].set_ylim(np.nanmin(list_bin_median) - 0.03, np.nanmax(list_bin_median) + 0.03)
            axs[count_ele - 1].axhline(y=np.nanmin(list_bin_median) - 0.03, color="black", linewidth=1.25)
            plot1 = axs[count_ele - 1].scatter(list_bin_steps, list_bin_median, s=16, zorder=3, linewidth=0.6, label="fg in F606W : (" + str(list_mag_range[ele_range]) + " to " + str(list_mag_range[ele_range + 1]) + ") Mag", color=list_color[count_ele - 1])
        else:
            axs[count_ele - 1].spines.top.set_visible(False)
            plot2 = axs[count_ele-1].scatter(list_bin_steps, list_bin_median, s=16, zorder=3, linewidth=0.6, label="fg in F606W : (" + str(list_mag_range[ele_range]) + " to " + str(list_mag_range[ele_range + 1]) + ") Mag", color=list_color[count_ele - 1])
        axs[count_ele-1].grid(linestyle="-.", linewidth=0.15, zorder=1)
        axs[count_ele-1].set_xlim(-5, 310)

        # plt.scatter(list_bin_steps, list_bin_median, zorder=2, color=list_color[count_ele - 1], s=5, label="Mass fg : (" + str(list_mass_range[ele_range]) + " to " + str(list_mass_range[ele_range + 1]) + ") log(M/M$_{\u2609}$)")
        # plt.scatter(list_bin_steps, list_bin_median, zorder=2, color=list_color[count_ele - 1], s=5, label="fg in F606W : (" + str(list_mag_range[ele_range]) + " to " + str(list_mag_range[ele_range + 1]) + ") Mag")
        # plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F606W - F775W [mag]"), plt.ylim(min(list_bin_median) - 0.075, max(list_bin_median) + 0.075)

    axs[0].legend((plot1, plot2), ("fg in F606W : (24 bis 28) Mag", "fg in F606W : (28 bis 32) Mag"), loc="upper right")
    # plt.savefig(image_path + "Multi", dpi=200)
    plt.show()