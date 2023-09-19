import time

import numpy as np
import math
from astropy.io import fits
import matplotlib.pyplot as plt

var_bg_count = 1
var_dist_min = 10
var_dist_kpc = 250
var_dist_arc = 100
var_snr = 15
var_mag = 24
var_mass = 0
var_mode = "gb"

if var_mode == "ab":
    var_z_max = 0.41
else:
    var_z_max = 0.91

image_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\Images\\w7\\"
fit_main_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\UDF_HUDF_Catalogs_FULL\\dr2_main_09.fits"
region_path_1 = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\objects_region_1.reg"

# function kpc/arcsec at redshift z, conv_func(z)
a, b, c, d, e, f, g, h, o, j, k, l, m, n = [-24.862278979480052, 21.267147446426016, 0.005306164849709134, 19.600027472637663, -11.409572942008737, 4.751202039434343, -1.3170007020351893, 0.2021612052937083, -0.007249014401993737, 0.0016127139142351095, -0.00017642802527037128, 9.986804760095554e-06, -2.1958087984210742e-07, -1.148703972276188e-09]

def conv_func(x):
    return n * x**14 + m * x**13 + l * x**12 + k * x**11 + j * x**10 + o * x**9 + h * x**7 + g * x**6 + f * x**5 + e * x**4 + d * x**3 + a * x**2 + b * x + c

outfile1 = open(region_path_1, "w")
fits_info = fits.open(fit_main_path)

outfile1.write('# Region file format: DS9 version 4.1\n')
outfile1.write('global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
outfile1.write('fk5\n')

ra = fits_info[1].data['RA']
# declination
dec = fits_info[1].data['DEC']
# red shift
z = fits_info[1].data['Z']
# redshift confidence
zconf = fits_info[1].data["ZCONF"]
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

list_val_dist, list_val_Ba, list_val_Ba_med, list_val_dist_med = [], [], [], []
count_hahb, count_hghb = 0, 0

for i in range(len(ra)):
    ri, di, zi, mag606_i, mass_i = ra[i], dec[i], z[i], mag606[i], mass[i]
    hai, hbi, hgi, ha_snri, hb_snri, hg_snri = h_a[i], h_b[i], h_g[i], h_a_snr[i], h_b_snr[i], h_g_snr[i]
    list_dist, list_id, list_close, list_med_dist, list_med_Ba = [], [], [], [], []
    ha_ab_i, hb_ab_i, hg_ab_i = ha_ab[i], hb_ab[i], hg_ab[i]

    if var_mode == "ab":
        if var_snr < ha_snri < 1e+15 and var_snr < hb_snri < 1e+15:
            count_hahb = count_hahb + 1
            for y in range(len(ra)):
                ry, dy, zy, mag606_y, mass_y = ra[y], dec[y], z[y], mag606[y], mass[y]
                dist_kpc = ((ri - ry) ** 2 + (di - dy) ** 2) ** (1 / 2) * 3600 * conv_func(zy)

                if i != y and zi > zy and mag606_y < var_mag and 0.01 < zy < var_z_max and var_dist_min < dist_kpc < var_dist_kpc:
                    if hai + ha_ab_i > 0 and hbi + hb_ab_i > 0:
                        list_dist.append(dist_kpc)
                        list_id.append(y)

    if var_mode == "gb":
        if var_snr < hg_snri < 1e+15 and var_snr < hb_snri < 1e+15:
            count_hghb = count_hghb + 1
            for y in range(len(ra)):
                ry, dy, zy, mag606_y, mass_y = ra[y], dec[y], z[y], mag606[y], mass[y]
                dist_kpc = ((ri - ry) ** 2 + (di - dy) ** 2) ** (1 / 2) * 3600 * conv_func(zy)

                if i != y and zi > zy and mag606_y < var_mag and 0.01 < zy < var_z_max and var_dist_min < dist_kpc < var_dist_kpc:
                    if hgi + hg_ab_i > 0 and hbi + hb_ab_i > 0:
                        list_dist.append(dist_kpc)
                        list_id.append(y)
    try:
        list_close.append(list_id[list_dist.index(min(list_dist))])
        list_id.pop(list_dist.index(min(list_dist)))
        list_dist.pop(list_dist.index(min(list_dist)))
    except ValueError:
        pass

    for ele in list_close:
        rc, dc, zc, ha_c, hb_c, hg_c = ra[ele], dec[ele], z[ele], h_a[i], h_b[i], h_g[i]
        dist_kpc = ((ri - rc)**2 + (di - dc)**2)**(1/2) * 3600 * conv_func(zc)
        outfile1.write('circle(%s,%s,0.35")# color=green width=2 font="helvetica 14 normal roman" text={%s}\n' % (rc, dc, i + 1))

        if var_mode == "ab":
            list_val_dist.append(dist_kpc)
            # list_val_Ba.append(ha_c / (hb_c * 2.86))
            # list_med_Ba.append(ha_c / (hb_c * 2.86))
            list_val_Ba.append(-2.5 * math.log(ha_c / (hb_c * 2.85)))
            list_med_Ba.append(-2.5 * math.log(ha_c / (hb_c * 2.85)))
        if var_mode == "gb":
            list_val_dist.append(dist_kpc)
            # list_val_Ba.append(hg_c / (hb_c * 0.47))
            # list_med_Ba.append(hg_c / (hb_c * 0.47))
            list_val_Ba.append(-2.5 * math.log(hg_c / (hb_c * 0.47)))
            list_med_Ba.append(-2.5 * math.log(hg_c / (hb_c * 0.47)))
        list_med_dist.append(dist_kpc)

    list_val_Ba_med.append(np.median(list_med_Ba))
    list_val_dist_med.append(np.median(list_med_dist))

print(count_hahb, count_hghb)

# plt.scatter(list_val_dist, list_val_Ba, s=20)
# plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
# plt.xlabel("Entfernung in Vordergrundebene [kpc]")
# if var_mode == "ab":
#     # plt.ylabel("Hα/Hβ/2.85")
#     plt.ylabel("E(Hα - Hβ)")
#     plt.title("E(Hα - Hβ), nahe fg-bg")
# if var_mode == "gb":
#     # plt.ylabel("Hγ/Hβ/0.47")
#     plt.ylabel("E(Hγ - Hβ)")
#     plt.title("E(Hγ - Hβ), nahe fg-bg")
# plt.show()

count_pop = 0
for ele in list_val_Ba:
    # ab
    # if abs(ele - np.median(list_val_Ba)) / 3 > np.std(list_val_Ba):
    # gb
    if abs(ele - np.median(list_val_Ba)) / 2 > np.std(list_val_Ba):
        pop_index = list_val_Ba.index(ele)
        list_val_Ba.pop(pop_index)
        list_val_dist.pop(pop_index)
        count_pop = count_pop + 1
print(count_pop, "count_pop")
# ab
# list_even = np.linspace(3, 21, 7)
# gb
list_even = np.linspace(2.3, 20.3, 7)
list_even_step = [int(ele3 * 10) for ele3 in list_even]
print(list_even_step)

bin_arr = np.linspace(1, 2300, 2300, dtype=int)
bin_lst = bin_arr.tolist()
direc = {}
for ele in bin_lst:
    direc["kpc_%s" % ele] = []
for ele in range(len(list_val_Ba)):
    var_loop = 0
    for ele2 in list_even_step:
        if list_val_dist[ele] < ele2 and var_loop == 0:
            direc["kpc_" + str(ele2)].append(list_val_Ba[ele])
            var_loop = 1


list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
for ele in list_even_step:
    list_bin_median.append(np.median(direc["kpc_" + str(ele)]))
    list_bin_mean.append(np.mean(direc["kpc_" + str(ele)]))
    list_bin_std.append(np.std(direc["kpc_" + str(ele)]))
    list_bin_steps.append(ele)


# plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
# plt.plot(list_bin_steps, list_bin_mean, zorder=2, label="mean", linewidth=1.1), plt.scatter(list_bin_steps, list_bin_median, zorder=3, label="median", color="xkcd:jade green")
# plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.legend()
# if var_mode == "ab":
#     plt.ylabel("Hα/Hβ/0.47")
# if var_mode == "gb":
#     plt.ylabel("Hγ/Hβ/0.47")
# plt.show()

list_steps_mid = [ele5 - 15 for ele5 in list_bin_steps]
list_xerr_width = [ele6 * 0 + 15 for ele6 in list_bin_steps]

# plt.rcParams.update({'font.size': 16})
# plt.figure(figsize=(10, 6))
plt.scatter(list_val_dist, list_val_Ba, s=20)
plt.hlines(0, -100, 400, linestyle="-.", linewidth=0.4, color="indianred")
plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
plt.errorbar(list_steps_mid, list_bin_mean, yerr=list_bin_std, xerr=list_xerr_width, fmt="none", ecolor="xkcd:dark", elinewidth=0.9, capsize=2, zorder=1)
plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylim(np.nanmin(list_bin_mean)-0.3, np.nanmax(list_bin_mean)+0.3),
plt.xlim(-10, 190)
# plt.xlim(-10, 200)
if var_mode == "ab":
    # plt.ylabel("Hα/Hβ/2.85")
    plt.ylabel("E(Hα - Hβ)")
    # plt.title("E(Hα - Hβ), Mittelwert und Median")
if var_mode == "gb":
    # plt.ylabel("Hγ/Hβ/0.47")
    plt.ylabel("E(Hγ - Hβ)")
    # plt.title("E(Hγ - Hβ), Mittelwert und Median")
plt.show()
# plt.savefig(image_path + "HgHb_bg-near_mean-hline", dpi=150)

outfile1.close()