import math
import time

import numpy as np
from astropy import units as u
from astropy.io import fits
import matplotlib.pyplot as plt
from mpdaf.obj import Cube
from sympy.solvers import solve
from sympy import Symbol


var_dist_kpc = 60
var_dist_kpc_min = 0
# variable, lower limit for redshift difference
var_z_diff = 0.00001
# variable, limit min and max redshift for galaxies
var_z_max = 0.42
var_z_min = 0
# variable, lower flux limit for at least 1 galaxy
var_flux = 0
var_mag = 25
var_mass = 0
var_plot = "n"
var_mode = "scan"

fac = 3
var_steps = 15
list_3 = np.linspace(-20, 20, 10)

# function kpc/arcsec at redshift z, conv_func(z)
a, b, c, d, e, f, g, h, o, j, k, l, m, n = [-24.862278979480052, 21.267147446426016, 0.005306164849709134, 19.600027472637663, -11.409572942008737, 4.751202039434343, -1.3170007020351893, 0.2021612052937083, -0.007249014401993737, 0.0016127139142351095, -0.00017642802527037128, 9.986804760095554e-06, -2.1958087984210742e-07, -1.148703972276188e-09]

def conv_func(x):
    return n * x ** 14 + m * x ** 13 + l * x ** 12 + k * x ** 11 + j * x ** 10 + o * x ** 9 + h * x ** 7 + g * x ** 6 + f * x ** 5 + e * x ** 4 + d * x ** 3 + a * x ** 2 + b * x + c

image_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\Images, Ba-Spec\\"
fit_main_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\UDF_HUDF_Catalogs_FULL\\dr2_main_09.fits"
cube_UDF_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\DATACUBE_UDF-10.fits"
cube_MXDF_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\DATACUBE_MXDF.fits"
cube_MOSAIC_path = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\AMUSED\\DATACUBE_UDF-MOSAIC.fits"
region_path_3 = "C:\\Users\\Kilian Gohlke\\Documents\\Uni-Potsdam\\Bachelorarbeit\\Python\\objects_region_3.reg"

cube_UDF, cube_MXDF, cube_MOSAIC = Cube(cube_UDF_path, ext=1), Cube(cube_MXDF_path, ext=1), Cube(cube_MOSAIC_path, ext=1)

outfile3 = open(region_path_3, 'w')
outfile3.write('# Region file format: DS9 version 4.1\n')
outfile3.write(
    'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
outfile3.write('fk5\n')

fits_info = fits.open(fit_main_path)

# right ascension
ra = fits_info[1].data['RA']
# declination
dec = fits_info[1].data['DEC']
# red shift
z = fits_info[1].data['Z']
# max flux
flux = fits_info[1].data["FLUX_MAX"]
# Dataset of galaxy
dat = fits_info[1].data["DATASET"]
# mag in F606W
mag606 = fits_info[1].data["MAG_F606W"]
# mass in log(M/M_sun)
mass = fits_info[1].data["MASS_PRO"]
# confidence redshift
zconf = fits_info[1].data["ZCONF"]

list_hghb, list_hghb_dist, list_hahb,list_hahb_dist = [], [], [], []
list_co_gb, list_xa_gb, list_co_ab, list_xa_ab, list_co_blre, list_co_BV = [], [], [], [], [], []
list_id = [192, 284]

# for i in range(len(ra)):
#     ri, di, zi, dat_i, mag_i, mass_i, zconfi = ra[i], dec[i], float(z[i]), dat[i], mag606[i], mass[i], zconf[i]
#     if zi < 0.01:
#         print(zi, mass_i, zconfi)
#
#         outfile3.write('circle(%s,%s,0.48") # color=yellow width=2 dash=1 font="helvetica 10 normal roman" text={%s} \n' % (ri, di, i))

# for ele9 in range(len(ra)):
#     if ele9 == 192:
#         zele = z[ele9]
#         for y in range(len(ra)):
#             ry, dy, zy, mag_y = ra[y], dec[y], float(z[y]), mag606[y]
#             if zy < zele:
#                 outfile3.write('circle(%s,%s,0.48") # color=yellow width=2 font="helvetica 10 normal roman" text={%s} \n' % (ry, dy, y))


# goes through list table, saves right ascension and declination for every object in var r, d
for i in range(len(ra)):
    ri, di, zi, dat_i, mag_i = ra[i], dec[i], float(z[i]), dat[i], mag606[i]

    for y in range(len(ra)):
        ry, dy, zy, mag_y = ra[y], dec[y], float(z[y]), mag606[y]

        if i != y and zi < zy and 0.01 < zi < var_z_max and 0.01 < zy < var_z_max and i in list_id and y in list_id:
        # if i != y and zi < zy and 0.01 < zi < var_z_max and 0.01 < zy < var_z_max:
            dist = ((ry - ri)**2 + (dy - di)**2)**(1/2)
            dist_kpc = ((ry - ri)**2 + (dy - di)**2)**(1/2) * 3600 * conv_func(zi)
            print(dist_kpc, zy)

            if var_dist_kpc_min < dist_kpc < var_dist_kpc and mag_i < var_mag and mag_y < var_mag:
                outfile3.write('circle(%s,%s,0.6") # color=white width=2 font="helvetica 10 normal roman" \n' % (ri, di))
                # outfile3.write('circle(%s,%s,0.8") # color=blue width=2 font="helvetica 10 normal roman" text={%s} \n' % (ry, dy, y))

                m1 = (dy - di) / (ry - ri)
                b1 = di - ri * m1

                list_hghb, list_hghb_dist, list_hahb, list_hahb_dist, list_HgHb_std, list_HaHb_std, list_angl = [], [], [], [], [], [], []

                for ele in list_3:
                    # find line with angle in list_3
                    m2 = (m1 - math.tan(ele * math.pi / 180)) / (1 + math.tan(ele * math.pi / 180) * m1)
                    b2 = di - m2 * ri
                    p = Symbol("p")
                    # rs is point where distance on new line is dist between centers of galaxies
                    r_list = solve((p - ri) ** 2 + ((m2 * p + b2) - di) ** 2 - (dist * 1.5) ** 2)
                    rs, rg = r_list[0], r_list[1]

                    if abs(((rs - ry) ** 2 + ((m2 * rs + b2) - dy) ** 2)) ** (1 / 2) < abs(((rg - ry) ** 2 + ((m2 * rg + b2) - dy) ** 2)) ** (1 / 2):
                        for ele3 in range(var_steps):
                            xc = ri - ele3 * (ri - rs) / var_steps
                            yc = m2 * xc + b2
                            list_angl.append(xc)
                            list_angl.append(yc)
                            list_angl.append(ele3)
                    else:
                        for ele3 in range(var_steps):
                            xc = ri - ele3 * (ri - rg) / var_steps
                            yc = m2 * xc + b2
                            list_angl.append(xc)
                            list_angl.append(yc)
                            list_angl.append(ele3)

                for ele2 in range(var_steps):
                    list_HgHb_id, list_HaHb_id, list_BV_id, list_VR_id, list_RI_id, list_BlRe_id = [], [], [], [], [], []
                    for ele in range(int(len(list_angl) / 3)):
                        ele_cord = ele * 3
                        if list_angl[ele_cord + 2] == ele2 and list_angl[ele_cord + 2] != 0 and list_angl[ele_cord + 2] != 1 and list_angl[ele_cord + 2] != 2 and var_mode == "scan":
                            cube = 0
                            if str(dat[i]) == "UDF10":
                                cube = cube_UDF
                            if str(dat[i]) == "MXDF":
                                cube = cube_MXDF
                            if str(dat[i]) == "MOSAIC":
                                cube = cube_MOSAIC
                            rcord, dcord = list_angl[ele_cord], list_angl[ele_cord + 1]
                            dist_i, dist_y = ((rcord - ri)**2 + (dcord - di)**2)**(1/2), ((rcord - ry)**2 + (dcord - dy)**2)**(1/2)
                            if dist_i > dist_y:
                                outfile3.write('circle(%s,%s,0.2") # color=red width=2 font="helvetica 10 normal roman" \n' % (rcord, dcord))
                                try:
                                    spec = cube.aperture((dcord, rcord), 0.2)
                                except ValueError:
                                    cube = cube_MOSAIC
                                    spec = cube.aperture((dcord, rcord), 0.2)

                                dist2 = ((ri - rcord) ** 2 + (di - dcord) ** 2) ** (1 / 2)
                                ha, hg, hb = (1 + zy) * 6564, (1 + zy) * 4341, (1 + zy) * 4862
                                hg_flu, hb_flu, ha_flu = 0, 0, 0

                                if zy < 0.91:
                                    sp_hg = spec.subspec(lmin=(hg - 20), lmax=(hg + 20), unit=u.angstrom)
                                    sp_hg_mean, sp_hg_std = np.mean(sp_hg.data), np.std(sp_hg.data)
                                    if sp_hg_mean + fac * sp_hg_std < sp_hg.data.max():
                                        var_min_hg, var_max_hg, hg_maskmin, hg_maskmax, hg_lmin, hg_lmax = 0, 0, hg - 5, hg + 5, 41 - 5, 41 + 5
                                        for ele6 in range(15):
                                            if var_min_hg == 0:
                                                if sp_hg.subspec(hg - ele6 - 1) <= np.median(sp_hg.data):
                                                    hg_lmin, hg_maskmin, var_min_hg = 15 - ele6 - 1, hg - ele6 - 2, 1
                                            if var_max_hg == 0:
                                                if sp_hg.subspec(hg + ele6 + 1) <= np.median(sp_hg.data):
                                                    hg_lmax, hg_maskmax, var_max_hg = 15 + ele6 + 1, hg + ele6 + 1, 1

                                        sp_hg.mask_region(hg_maskmin, hg_maskmax)
                                        fit_hg = sp_hg.poly_spec(6)
                                        sp_hg.unmask()
                                        sp_hg_wico = sp_hg - fit_hg
                                        hg_flu = sp_hg_wico.sum()[0] * sp_hg_wico.get_step(unit=sp_hg_wico.wave.unit)
                                        if var_plot == "y":
                                            sp_hg.plot(linewidth=0.5)
                                            fit_hg.plot(color="red")
                                            sp_hg_wico[int(hg_lmin):int(hg_lmax)].plot(color="orange", linewidth=2.0)
                                            sp_hg_wico.plot(color="green")
                                            plt.title("hg")
                                            plt.show()

                                    sp_hb = spec.subspec(lmin=(hb - 20), lmax=(hb + 20), unit=u.angstrom)
                                    sp_hb_mean, sp_hb_std = np.mean(sp_hb.data), np.std(sp_hb.data)
                                    if sp_hb_mean + fac * sp_hb_std < sp_hb.data.max():
                                        var_min_hb, var_max_hb, hb_maskmin, hb_maskmax, hb_lmin, hb_lmax = 0, 0, hb - 5, hb + 5, 41 - 5, 41 + 5
                                        for ele6 in range(15):
                                            if var_min_hb == 0:
                                                if sp_hb.subspec(hb - ele6 - 1) <= np.median(sp_hb.data):
                                                    hb_lmin, hb_maskmin, var_min_hb = 15 - ele6 - 1, hb - ele6 - 1, 1
                                            if var_max_hb == 0:
                                                if sp_hb.subspec(hb + ele6 + 1) <= np.median(sp_hb.data):
                                                    hb_lmax, hb_maskmax, var_max_hb = 15 + ele6 + 1, hb + ele6 + 1, 1
                                        sp_hb.mask_region(hb_maskmin, hb_maskmax)
                                        fit_hb = sp_hb.poly_spec(6)
                                        sp_hb.unmask()
                                        sp_hb_wico = sp_hb - fit_hb
                                        hb_flu = sp_hb_wico.sum()[0] * sp_hb_wico.get_step(unit=sp_hb_wico.wave.unit)
                                        if var_plot == "y":
                                            sp_hb.plot()
                                            fit_hb.plot(color="red")
                                            sp_hb_wico[int(hb_lmin):int(hb_lmax)].plot(color="orange", linewidth=2.0)
                                            sp_hb_wico.plot(color="green")
                                            plt.title("hb")
                                        plt.show()

                                    if hg_flu > 0 and hb_flu > 0:
                                        if 0 < hg_flu / hb_flu < 5:
                                            # list_HgHb_id.append(hg_flu / (hb_flu * 0.47))
                                            list_HgHb_id.append(-2.5 * math.log(hg_flu / (hb_flu * 0.47)))
                                            print("hg/hb", hg_flu / hb_flu)

                                if zy < 0.41:
                                    sp_ha = spec.subspec(lmin=(ha - 20), lmax=(ha + 20), unit=u.angstrom)
                                    sp_ha_mean, sp_ha_std = np.mean(sp_ha.data), np.std(sp_ha.data)

                                    if sp_ha_mean + fac * sp_ha_std < sp_ha.data.max():
                                        var_min_ha, var_max_ha, ha_maskmin, ha_maskmax, ha_lmin, ha_lmax = 0, 0, ha - 5, ha + 5, 41 - 5, 41 + 5
                                        for ele6 in range(15):
                                            if var_min_ha == 0:
                                                if sp_ha.subspec(ha - ele6 - 1) <= np.median(sp_ha.data):
                                                    ha_lmin, ha_maskmin, var_min_ha = 15 - ele6 - 1, ha - ele6 - 1, 1
                                            if var_max_ha == 0:
                                                if sp_ha.subspec(ha + ele6 + 1) <= np.median(sp_ha.data):
                                                    ha_lmax, ha_maskmax, var_max_ha = 15 + ele6 + 1, ha + ele6 + 1, 1
                                        sp_ha.mask_region(ha_maskmin, ha_maskmax)
                                        fit_ha = sp_ha.poly_spec(5)
                                        sp_ha.unmask()
                                        sp_ha_wico = sp_ha - fit_ha
                                        ha_flu = sp_ha_wico[int(ha_lmin):int(ha_lmax)].sum()[0] * sp_ha_wico.get_step(unit=sp_ha_wico.wave.unit)
                                        if var_plot == "y":
                                            sp_ha.plot()
                                            fit_ha.plot(color="red")
                                            sp_ha_wico[int(ha_lmin):int(ha_lmax)].plot(color="orange", linewidth=2.0)
                                            sp_ha_wico.plot(color="green")
                                            plt.title("ha")
                                            plt.show()

                                    if ha_flu > 0 and hb_flu > 0:
                                        if 0 < ha_flu / hb_flu < 15:
                                            # list_HaHb_id.append(ha_flu / (hb_flu * 2.86))
                                            list_HaHb_id.append(-2.5 * math.log(ha_flu / (hb_flu * 2.85)))

                                            print("ha/hb", ha_flu / hb_flu)

                    if list_HgHb_id and len(list_HgHb_id) != 1:
                        list_hghb.append(np.median(list_HgHb_id))
                        list_HgHb_std.append(np.std(list_HgHb_id))
                        list_hghb_dist.append(float(dist2 * 3600 * conv_func(zi)))
                    if list_HaHb_id and len(list_HaHb_id) != 1:
                        list_hahb.append(np.median(list_HaHb_id))
                        list_HaHb_std.append(np.std(list_HaHb_id))
                        list_hahb_dist.append(float(dist2 * 3600 * conv_func(zi)))

                if len(list_hghb) > 0:
                    plt.errorbar(list_hghb_dist, list_hghb, yerr=list_HgHb_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
                    plt.grid(linestyle="-.", linewidth=0.25, zorder=1)
                    co = np.polyfit(list_hghb_dist, list_hghb, 1)
                    linfit = np.poly1d(co)
                    x_a = np.linspace(min(list_hghb_dist), max(list_hghb_dist), 50)
                    list_co_gb.append(co[0])
                    plt.scatter(list_hghb_dist, list_hghb)
                    # plt.plot(x_a, linfit(x_a), color="dimgrey", linewidth=0.94, linestyle="dashed", zorder=3)
                    plt.xlabel("distance to foreground galaxy in kpc"), plt.ylabel("E(Hγ-Hβ)"), plt.title("E(Hγ-Hβ), scan")
                    plt.savefig(image_path + "HgHb-" + str(i + 1) + "-" + str(y + 1) + ".png")
                    plt.close()
                if len(list_hahb) > 3:
                    plt.errorbar(list_hahb_dist, list_hahb, yerr=list_HaHb_std, fmt="none", ecolor="indianred", elinewidth=0.83, capsize=4, zorder=1)
                    plt.grid(linestyle="-.", linewidth=0.25, zorder=1)
                    co = np.polyfit(list_hahb_dist, list_hahb, 1)
                    linfit = np.poly1d(co)
                    x_a = np.linspace(min(list_hahb_dist), max(list_hahb_dist), 50)
                    list_co_ab.append(co[0])
                    plt.scatter(list_hahb_dist, list_hahb)
                    # plt.plot(x_a, linfit(x_a), color="dimgrey", linewidth=0.94, linestyle="dashed", zorder=3),
                    # plt.xlabel("distance to foreground galaxy in kpc"), plt.ylabel("Ha/Hb/2.85"), plt.title("Ha/Hb-dist(kpc)")
                    plt.xlabel("distance to foreground galaxy in kpc"), plt.ylabel("E(Hα-Hβ)"), plt.title("E(Hα-Hβ), scan")
                    plt.savefig(image_path + "HaHb-" + str(i + 1) + "-" + str(y + 1) + ".png")
                    plt.close()


