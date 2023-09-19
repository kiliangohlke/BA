import os
import numpy as np
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from mpdaf.obj import Cube
import time

# variable, upper limit for projected distance
var_dist = 30
var_dist_kpc = 180
var_dist_kpc_min = 120
# variable, limit min and max redshift for galaxies
var_z_max = 12
var_z_min = 0.6
# variable, lower flux limit for at least 1 galaxy
var_flux = 0
fac = 0
var_ang = 1
var_mag_i = 26
var_mag_y = 26
var_mass_y = 0
var_c = 513

var_mode = "div"

var_rej = 2
rej_out = "y"
mask_lines = "y"
var_obs = "fg"


txt_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\Spectra-med\\"
txt_path2 = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\Spectra_Restframe\\v3_b\\"
# txt_path2 = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\Spectra_Restframe\\"
# txt_path2 = "D:\\Spectra_Restframe\\"
fit_main_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\UDF_HUDF_Catalogs_FULL\\dr2_main_09.fits"
region_path_4 = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\objects_region_4.reg"

outfile4 = open(region_path_4, 'w')
outfile4.write('# Region file format: DS9 version 4.1\n')
outfile4.write(
    'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
outfile4.write('fk5\n')

# function kpc/arcsec at redshift z, conv_func(z)
a, b, c, d, e, f, g, h, o, j, k, l, m, n = [-24.862278979480052, 21.267147446426016, 0.005306164849709134, 19.600027472637663, -11.409572942008737, 4.751202039434343, -1.3170007020351893, 0.2021612052937083, -0.007249014401993737, 0.0016127139142351095, -0.00017642802527037128, 9.986804760095554e-06, -2.1958087984210742e-07, -1.148703972276188e-09]

def conv_func(x):
    return n * x**14 + m * x**13 + l * x**12 + k * x**11 + j * x**10 + o * x**9 + h * x**7 + g * x**6 + f * x**5 + e * x**4 + d * x**3 + a * x**2 + b * x + c

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
# H-beta emission flux line
h_b = fits_info[1].data["HBETA_EMI_FLUX"]
# H_gamma emission flux line
h_g = fits_info[1].data["HGAMMA_EMI_FLUX"]
# Dataset of galaxy
dat = fits_info[1].data["DATASET"]
# confidence redshift
zconf = fits_info[1].data["ZCONF"]
# mag in F606W
mag606 = fits_info[1].data["MAG_F606W"]
# log(M) / M_sun
mass = fits_info[1].data["MASS_PRO"]

count, count_HaHb, count_HgHb = 0, 0, 0

# 6564, 4862, 4341, 2801, 2799, 1548, 1550, 2600, 2586, 2383, 2375, 2344, 1855, 1807, 1670, 1524, 1393, 1301, 1258, 1206, 1608,1906, 1403, 2325, 2396, 2625, 2611
if mask_lines == "y":
    list_lines_sort, list_sort, list_lines = [], np.linspace(-10, 10, 21), [6564, 4862, 4341, 2801, 2799, 1548, 1550, 2600, 2586, 2383, 2375, 2344, 1855, 1807, 1670, 1524, 1393, 1301, 1258, 1206, 1608,1906, 1403, 2325, 2396, 2625, 2611, 1333, 2365, 2793, 1215]
else:
    list_lines_sort, list_sort, list_lines = [], np.linspace(-3, 3, 7), []
for eles in list_lines:
    for eles2 in list_sort:
        list_lines_sort.append(eles + eles2)
print(list_lines_sort)

bin_arr = np.linspace(1, 10000, 10000, dtype=int)
bin_lst = bin_arr.tolist()
direc = {}
for ele in bin_lst:
    direc["list_%s" % ele] = []

for i in range(len(ra)):
    ri, di, zi, flux_i, mag606_i = ra[i], dec[i], float(z[i]), flux[i], mag606[i]

    for y in range(len(ra)):
        # goes through list table, saves right ascension and declination, flux, emission lines for
        ry, dy, zy, zconfy, mag606_y, mass_y = ra[y], dec[y], float(z[y]), zconf[y], mag606[y], mass[y]

        if i != y:
            # prevents the pairing of galaxies with themselves
            if zi > zy and zconfy > 1 and mag606_i < var_mag_i and mag606_y < var_mag_y and mass_y > var_mass_y:
                if var_z_min < zy < var_z_max and 0.1 < zi - zy < 10.0:
                    # distance to foreground galaxy
                    dist = ((ry - ri) ** 2 + (dy - di) ** 2) ** (1 / 2)
                    dist_kpc = dist * 3600 * conv_func(zy)


                    if var_dist_kpc_min < dist_kpc < var_dist_kpc and dist * 3600 < 100:

                        var_fg = 0
                        for elez in range(len(ra)):
                            rz, dz, zz = ra[elez], dec[elez], float(z[elez])
                            dist3 = (((rz - ri) ** 2 + (dz - di) ** 2) ** (1 / 2)) * 3600
                            if elez != i and dist3 < var_ang and zz < zi:
                                var_fg = 1

                        if var_fg == 0:
                            outfile4.write('circle(%s,%s,0.5") # color=red width=2 font="helvetica 10 normal roman" \n' % (ri, di))
                            outfile4.write('circle(%s,%s,0.5") # color=green width=2 font="helvetica 10 normal roman" \n' % (ry, dy))
                            count = count + 1
                            print("")
                            print("Pair: "+str(count), "i: " + str(i+1), "y: " + str(y+1), "dist-kpc", dist_kpc, zy)

                            # with open(txt_path2 + str(i) + "sub_v2.txt", "r") as file:
                            # with open(txt_path2 + str(i) + "_sub.txt", "r") as file:
                            # with open(txt_path2 + str(i) + ".txt", "r") as file:
                            if var_mode == "div":
                                with open(txt_path2 + str(i) + "_div_v3.txt", "r") as file:
                                    lines = []
                                    for line in file:
                                        lines.append(float(line.rstrip()))
                                    spec_std, spec_mean, sp_st = float(lines[0]), float(lines[1]), int(lines[2])
                                    lines.pop(0), lines.pop(0), lines.pop(0)
                                    for ele in range(len(lines)):
                                        ele_step = ele * 1.25
                                        if round((sp_st + ele_step) / (1 + zi)) in list_lines_sort:
                                            # print(round((sp_st + ele_step) / (1 + zi)))
                                            pass
                                        else:
                                            if var_obs == "fg":
                                                direc["list_" + str(round((sp_st + ele_step) / (1 + zy)))].append(
                                                    lines[ele])
                                            if var_obs == "bg":
                                                direc["list_" + str(round((sp_st + ele_step) / (1 + zi)))].append(
                                                    lines[ele])
                            if var_mode == "sub":
                                with open(txt_path2 + str(i) + "_sub_v3.txt", "r") as file:
                                    lines = []
                                    for line in file:
                                        lines.append(float(line.rstrip()))
                                    spec_std, spec_mean, sp_st = float(lines[0]), float(lines[1]), int(lines[2])
                                    lines.pop(0), lines.pop(0), lines.pop(0)
                                    for ele in range(len(lines)):
                                        ele_step = ele * 1.25
                                        if round((sp_st + ele_step) / (1 + zi)) in list_lines_sort:
                                            # print(round((sp_st + ele_step) / (1 + zi)))
                                            pass
                                        else:
                                            if var_obs == "fg":
                                                direc["list_" + str(round((sp_st + ele_step) / (1 + zy)))].append(lines[ele])
                                            if var_obs == "bg":
                                                direc["list_" + str(round((sp_st + ele_step) / (1 + zi)))].append(lines[ele])
dmed, direc_place, med_direc = [], [], 0
list_direc_mean, list_direc_median, list_direc_std, list_direc_bin, list_direc_count, list_direc_sum, = [], [], [], [], [], []

if rej_out == "y":
    for ele in range(800, 6200):
        print(ele)
        dmed = []
        direc_place = direc["list_"+str(ele + 1)]
        med_direc = np.median(direc_place)

        for ele2 in direc_place:
            dmed.append(abs(ele2 - med_direc))

        med_dmed = np.median(dmed)
        direc_place_copy = direc_place.copy()

        for ele2 in direc_place_copy:
            if abs(ele2 - med_direc) > var_rej * med_dmed:
                index = direc_place.index(ele2)
                direc_place.pop(index)

for ele in range(10000):
    list_direc_median.append(np.median(direc["list_"+str(ele + 1)]))
    list_direc_mean.append(np.mean(direc["list_"+str(ele + 1)]))
    list_direc_std.append(np.std(direc["list_" + str(ele + 1)]))
    list_direc_bin.append(ele + 1)
    list_direc_count.append(len(direc["list_"+str(ele + 1)]))
    list_direc_sum.append(sum(direc["list_" + str(ele + 1)]))

with open(txt_path+"Info, ("+str(var_c)+").txt", "w") as f:
    f.write("max dist deg: "+str(var_dist)+"\n")
    f.write("min dist kpc: " + str(var_dist_kpc_min) + "\n")
    f.write("max dist kpc: " + str(var_dist_kpc) + "\n")
    f.write("redshift range: "+str(var_z_min)+"-"+str(var_z_max)+"\n")
    f.write("min flux: "+str(var_flux)+"\n")
    f.write("number of pairs: "+str(count)+"\n")
    f.write("fac snr: "+str(fac)+"\n")
    f.write("ang separation: "+str(var_ang)+"\n")
    f.write("mag606_i: "+str(var_mag_i)+"\n")
    f.write("mag606_y: "+str(var_mag_y)+"\n")

with open(txt_path+"Spec-med, bin, ("+str(var_c)+").txt", "w") as f:
    for ele5 in range(len(list_direc_bin)):
        f.write(str(list_direc_bin[ele5]) + "\n")
with open(txt_path+"Spec-med, median, ("+str(var_c)+").txt", "w") as f:
    for ele5 in range(len(list_direc_bin)):
        f.write(str(list_direc_median[ele5]) + "\n")
with open(txt_path+"Spec-med, mean, ("+str(var_c)+").txt", "w") as f:
    for ele5 in range(len(list_direc_bin)):
        f.write(str(list_direc_mean[ele5]) + "\n")
with open(txt_path+"Spec-med, std, ("+str(var_c)+").txt", "w") as f:
    for ele5 in range(len(list_direc_bin)):
        f.write(str(list_direc_std[ele5]) + "\n")
with open(txt_path+"Spec-med, count, ("+str(var_c)+").txt", "w") as f:
    for ele5 in range(len(list_direc_bin)):
        f.write(str(list_direc_count[ele5]) + "\n")
with open(txt_path+"Spec-med, sum, ("+str(var_c)+").txt", "w") as f:
    for ele5 in range(len(list_direc_bin)):
        f.write(str(list_direc_sum[ele5]) + "\n")
with open(txt_path+"Spec-med, sum-rel, ("+str(var_c)+").txt", "w") as f:
    for ele5 in range(len(list_direc_bin)):
        f.write(str(list_direc_sum[ele5] - list_direc_count[ele5]) + "\n")
