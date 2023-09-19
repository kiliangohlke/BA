import time
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from mpdaf.obj import Cube

# strong Trend
# var_dist_kpc = 250
# var_err = 100
# var_mag = 100
# var_mass = 0
# var_zdiff_max = 1, 0.1

var_dist_kpc = 250
var_err = 0.3
var_mag = 25
var_mass = 0
var_z_max = 12
var_zdiff_max = 1

# function kpc/arcsec at redshift z, conv_func(z)
a, b, c, d, e, f, g, h, o, j, k, l, m, n = [-24.862278979480052, 21.267147446426016, 0.005306164849709134, 19.600027472637663, -11.409572942008737, 4.751202039434343, -1.3170007020351893, 0.2021612052937083, -0.007249014401993737, 0.0016127139142351095, -0.00017642802527037128, 9.986804760095554e-06, -2.1958087984210742e-07, -1.148703972276188e-09]

def conv_func(x):
    return n * x**14 + m * x**13 + l * x**12 + k * x**11 + j * x**10 + o * x**9 + h * x**7 + g * x**6 + f * x**5 + e * x**4 + d * x**3 + a * x**2 + b * x + c

image_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\Images_BV_bg-near\\"
fit_main_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\UDF_HUDF_Catalogs_FULL\\dr2_main_09.fits"
fits_info = fits.open(fit_main_path)
cube_UDF_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\DATACUBE_UDF-10.fits"
cube_MXDF_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\DATACUBE_MXDF.fits"
cube_MOSAIC_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\DATACUBE_UDF-MOSAIC.fits"
region_path_1 = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\objects_region_1.reg"
outfile1 = open(region_path_1, "w")
outfile1.write('# Region file format: DS9 version 4.1\n')
outfile1.write('global color=green dashlist=8 3 width=2 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
outfile1.write('fk5\n')

cube_UDF, cube_MXDF, cube_MOSAIC = Cube(cube_UDF_path, ext=1), Cube(cube_MXDF_path, ext=1), Cube(cube_MOSAIC_path, ext=1)

# right ascension
ra = fits_info[1].data['RA']
# declination
dec = fits_info[1].data['DEC']
# red shift
z = fits_info[1].data['Z']
# redshift confidence
zconf = fits_info[1].data["ZCONF"]
# Dataset of galaxy
dat = fits_info[1].data["DATASET"]
# mag in F606W
mag606 = fits_info[1].data["MAG_F606W"]
# mag error in F606W
mag606_err = fits_info[1].data["MAGERR_F606W"]
# mag in F775W
mag775 = fits_info[1].data["MAG_F775W"]
# mag error in F775W
mag775_err = fits_info[1].data["MAGERR_F775W"]
# mass in log(M/M_sun)
mass = fits_info[1].data["MASS_PRO"]

bin_arr = np.linspace(1, 2300, 2300, dtype=int)
bin_lst, direc, direc2 = bin_arr.tolist(), {}, {}
for ele in bin_lst:
    direc["list_%s" % ele] = []
    direc2["list_%s" % ele] = []

list_hahb_val, list_hahb_dist, list_hghb_val, list_hghb_dist, list_hahb_val2, list_hahb_dist2, list_hghb_val2, list_hghb_dist2, list_hahb_mag, list_hghb_mag = [], [], [], [], [], [], [], [], [], []

for i in range(len(ra)):
    ri, di, zi = ra[i], dec[i], z[i]
    mag606_i, mag775_i, mag606err_i, mag775err_i = mag606[i], mag775[i], mag606_err[i], mag775_err[i]
    list_dist, list_pairs = [], []
    for y in range(len(ra)):
        ry, dy, zy, mag606_y, mass_y = ra[y], dec[y], z[y], mag606[y], mag775[y]
        if i != y and zi > zy and 0.01 < zi < var_z_max and 0.01 < zy < var_z_max:
            dist = ((ri - ry) ** 2 + (di - dy) ** 2) ** (1 / 2) * 3600 * conv_func(zy)
            dist_arc = ((ri - ry) ** 2 + (di - dy) ** 2) ** (1 / 2)
            print("Pair", i, y, dist)

            if 0 < dist < var_dist_kpc and 0.01 < zi - zy < var_zdiff_max and var_mass < mass_y < 1e+15 and mag606_y < var_mag:
                # if 0.001 < mag606_i < var_mag and 0.001 < mag775_i < var_mag and mag606err_i < var_err and mag775err_i < var_err and mass_y > var_mass :
                list_pairs.append(y)
                list_dist.append(dist)

    try:
        min_index = list_dist.index(min(list_dist))
        # min_index = list_pairs_mag_ab.index(min(list_pairs_mag_ab))
        # mag606_i, mag775_i, mag606err_i, mag775err_i = mag606[list_pairs[min_index]], mag775[list_pairs[min_index]], mag606_err[list_pairs[min_index]], mag775_err[list_pairs[min_index]]
        r2, d2, z2 = ra[list_pairs[min_index]], dec[list_pairs[min_index]], z[list_pairs[min_index]]
        dist = ((ri - r2) ** 2 + (di - d2) ** 2) ** (1 / 2) * 3600 * conv_func(z2)
        if mag606err_i < var_err and mag775err_i < var_err:
            list_hahb_dist.append(dist)
            list_hahb_val.append(mag606_i - mag775_i)


        outfile1.write('circle(%s,%s,0.5")# color=red width=2 font="helvetica 14 normal roman" text={%s}\n' % (
            ri, di, i + 1))
        outfile1.write('circle(%s,%s,0.9")# color=green width=2 font="helvetica 14 normal roman" text={%s}\n' % (
            r2, d2, i + 1))
    except ValueError:
        pass

count_pop = 0
for ele in list_hahb_val:
    if abs(ele - np.median(list_hahb_val)) / 5 > np.std(list_hahb_val):
        pop_index = list_hahb_val.index(ele)
        list_hahb_val.pop(pop_index)
        list_hahb_dist.pop(pop_index)
        count_pop = count_pop + 1
print(count_pop, "count_pop")

plt.scatter(list_hahb_dist, list_hahb_val, s=10)
plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F606W - F775W [mag]")
plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
plt.savefig(image_path + "606-775_bg-near_scatter", dpi=150), plt.close()
# plt.show()

bin_arr = np.linspace(10, var_dist_kpc, int(var_dist_kpc / 10), dtype=int)
bin_lst = bin_arr.tolist()
direc = {}
for ele in bin_lst:
    direc["kpc_%s" % ele] = []
for ele in range(len(list_hahb_val)):
    direc["kpc_" + str(int(np.ceil(list_hahb_dist[ele] / 10)) * 10)].append(list_hahb_val[ele])

list_bin_mean, list_bin_median, list_bin_std, list_bin_steps, list_bin_count = [], [], [], [], []
for ele in range(int(var_dist_kpc / 10)):
    list_bin_median.append(np.median(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_mean.append(np.mean(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_std.append(np.std(direc["kpc_" + str((ele + 1) * 10)]))
    list_bin_steps.append((ele + 1) * 10)
    list_bin_count.append(len(direc["kpc_" + str((ele + 1) * 10)]))

plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83,capsize=4, zorder=1)
plt.grid(linestyle="-.", linewidth=0.2, zorder=1)
plt.plot(list_bin_steps, list_bin_mean, zorder=2, label="mean", linewidth=1.1), plt.scatter(list_bin_steps, list_bin_median, zorder=3, label="median", color="xkcd:jade green")
plt.xlabel("Entfernung in Vordergrundebene [kpc]"), plt.ylabel("F606W - F775W [mag]"), plt.legend(),
# plt.ylim(np.nanmin(list_bin_mean) - 0.35, np.nanmax(list_bin_mean) + 0.35)
plt.savefig(image_path + "606-775_bg-near_mean-median-err", dpi=150), plt.close()
# plt.show()

# # plt.errorbar(list_bin_steps, list_bin_mean, yerr=list_bin_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
# plt.scatter(list_bin_steps, list_bin_mean, zorder=1, label="Mean", s=24, edgecolors="none"), plt.scatter(list_bin_steps, list_bin_median, color="xkcd:jade green", label="Median", s=28, edgecolors="xkcd:jade green", facecolor="none")
# plt.show()