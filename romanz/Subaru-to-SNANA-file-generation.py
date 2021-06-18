"""Subaru-to-SNANA-file-generation.py

A script to convert form a set of Subaru PFS ETC output files to an SNANA
Spectrograph definition file.

Benjamin Rose
Jun 14, 2021
initially used Python 3.9, but should be compatible with >3.7
"""
# from astropy.table.groups import TableGroups
# import numpy as np
# import astropy
from astropy.table import Table, Column
from astropy.io import ascii, fits

# import matplotlib.pyplot as plt

# What parts
MAG = ["18", "28"]
EXP_TIME = [".16", "1", "5", "10"]
FILE = "ref.snc.{}-{}.dat"
Tables = {}

# esasier to keep these header formats of MAG and EXP_TIME hand coded
MAGREF_LIST = "22 28"
TEXPOSE_LIST = "600 3600 18000 36000"


# Read ECT data
###############
for mag in MAG:
    for et in EXP_TIME:
        Tables[mag + et] = ascii.read(
            FILE.format(mag, et),
            names=(
                "arm_num",
                "pixel",
                "wavelength",
                "S/N",
                "signal",
                "noise_var",
                "noise_var_1",
                "input_spectra",
                "conversion_factor",
                "sampling_factor",
                "background",
            ),
        )


# Verifying something that can mess up SNANA.
# I was told to just do it.
for x, y in zip(
    Tables[MAG[0] + EXP_TIME[0]]["wavelength"],
    Tables[MAG[1] + EXP_TIME[1]]["wavelength"],
):
    if x - y != 0:
        print(x - y)


# Reformat Data
###############
T = Table(
    names=(
        "lam min",
        "lam max",
        "lam res",
        "snr1t1",
        "snr2t1",
        "snr1t2",
        "snr2t2",
        "snr1t3",
        "snr2t3",
        "snr1t4",
        "snr2t4",
    )
)

c = []  # will be SNANA line key "SPECBIN:"
for x in range(0, len(Tables[MAG[0] + EXP_TIME[0]]["wavelength"]) - 1):
    T.add_row(
        [
            Tables[MAG[0] + EXP_TIME[0]]["wavelength"][x],
            Tables[MAG[0] + EXP_TIME[0]]["wavelength"][x + 1],
            0,
            # Itterate: for each expsoure time, print mag[0] then mag[1]
            # 10 mins
            Tables[MAG[0] + EXP_TIME[0]]["S/N"][x],
            Tables[MAG[1] + EXP_TIME[0]]["S/N"][x],
            # 1 hr
            Tables[MAG[0] + EXP_TIME[1]]["S/N"][x],
            Tables[MAG[1] + EXP_TIME[1]]["S/N"][x],
            # 5 hr
            Tables[MAG[0] + EXP_TIME[2]]["S/N"][x],
            Tables[MAG[1] + EXP_TIME[2]]["S/N"][x],
            # 10 hr
            Tables[MAG[0] + EXP_TIME[3]]["S/N"][x],
            Tables[MAG[1] + EXP_TIME[3]]["S/N"][x],
        ]
    )
    c.append("SPECBIN:")
c1 = Column(data=c, name="#")
print(T)


# Save Data
##########
# FITs Table
T.write("test_spec_subaru.fits", format="fits", overwrite=True)
# update header, https://docs.astropy.org/en/stable/io/fits/index.html#save-file-changes
with fits.open("test_spec_subaru.fits", mode="update") as tab:
    hdr = tab[0].header
    hdr.set("Instrument:", "Subaru+PSF")
    hdr.set("MAGREF_LIST:", MAGREF_LIST)
    hdr.set("TEXPOSE_LIST:", TEXPOSE_LIST)

# SNANA .dat file
T.add_column(c1, index=0)
with open("SUBARUPFS_table_20210614.dat", "w") as f:
    f.write(
        f"""INSTRUMENT: SUBARU+PFS
MAGREF_LIST: {MAGREF_LIST} #used to define SNR1 and SNR2
TEXPOSE_LIST: {TEXPOSE_LIST} #seconds
#     LAM LAM LAM
#     MIN MAX RES SNR1 SNR2 SNR1 SNR2 SNR1 SNR2 SNR1 SNR2
"""
    )
    # There has to be a better way, epically since this is the order of `T`!
    # But I inherited this code block and it works.
    for x in range(0, len(T["#"])):
        f.write(
            T["#"][x]
            + " "
            + str(T["lam min"][x])
            + " "
            + str(T["lam max"][x])
            + " "
            + str(T["lam res"][x])
            + " "
            + str(T["snr1t1"][x])
            + " "
            + str(T["snr2t1"][x])
            + " "
            + str(T["snr1t2"][x])
            + " "
            + str(T["snr2t2"][x])
            + " "
            + str(T["snr1t3"][x])
            + " "
            + str(T["snr2t3"][x])
            + " "
            + str(T["snr1t4"][x])
            + " "
            + str(T["snr2t4"][x])
            + " \n"
        )
