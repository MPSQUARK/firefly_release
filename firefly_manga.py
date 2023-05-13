"""
Created on Tue May 12 11:00:35 2020

.. moduleauthor:: Justus Neumann <jusneuma.astro__at__gmail.com>
.. contributions:: Daniel Goddard <daniel.goddard__at__port.ac.uk>

General purpose:
................

Reads in a MaNGA datacube and analyses each spectrum from Voronoi binned spectra.
The bin number is optional and should only be given if only a single bin is meant to be fitted.

"""

import sys
import os
import time
import multiprocessing
from multiprocessing import Pool
import numpy as np
from astropy.io import fits
from src.config import Config
from src.setup import Setup
from src.models import StellarPopulationModel

sys.path.append(os.path.join(os.getcwd(), "python"))
os.environ["FF_DIR"] = os.getcwd()
os.environ["STELLARPOPMODELS_DIR"] = os.path.join(
    os.environ["FF_DIR"], "stellar_population_models"
)

t0 = time.time()

config = Config().create("config/manga.yaml")

# paths to input files
logcube_dir = "example_data/manga/"
maps_dir = "example_data/manga/"
output_dir = "output/manga/"
dap_file = "example_data/manga/dapall-v3_1_1-3.1.0.fits"  # For DR15 use dapall-v2_4_3-2.2.1.fits

# define plate number, IFU number, bin number
plate = 8080
ifu = 12701
bin_number = 0  #'all' if entire IFU, otherwise bin number

# define MaNGA Product Launch (MPL) to use
# mpl-11 = DR17 final data realease (recommended)
# mpl-7 = DR15
mpl = "mpl-11"


def f(i):
    galaxy_bin_number = i
    print(f"Fitting bin number {i}")
    output_file = f"{direc}/spFly-{plate}-{ifu}-bin{i}"
    spec = Setup(
        maps,
        config.milky_way_reddening,
        config.hpf_mode,
        config.n_angstrom_masked,
    )
    spec.openMANGASpectrum(
        logcube, 
        dap_file, 
        galaxy_bin_number, 
        plate, 
        ifu, 
        config.emlines, 
        mpl
    )

    # prepare model templates
    model = StellarPopulationModel(
        spec,
        output_file,
        config.cosmo,
        models=config.model_key,
        model_libs=config.model_lib,
        imfs=config.imfs,
        age_limits= config.age_limits,
        downgrade_models=True,
        data_wave_medium=config.data_wave_medium,
        Z_limits=config.z_limits,
        suffix=config.suffix,
        use_downgraded_models=False,
        dust_law=config.dust_law,
        max_ebv=config.max_ebv,
        num_dust_vals=config.num_dust_vals,
        dust_smoothing_length=config.dust_smoothing_length,
        max_iterations=config.max_iterations,
        pdf_sampling=config.pdf_sampling,
        flux_units=config.flux_units,
    )
    # initiate fit
    model.fit_models_to_data()

def f_mp(unique_bin_number):  # multiprocessing function
    if __name__ == "__main__":
        pool = Pool(
            processes=multiprocessing.cpu_count()
        )  # Define number of cores to use
        pool.map(f, unique_bin_number)
        pool.close()
        pool.join()

# ------------------------------
# START Running firefly

print("\nStarting firefly ...")
print(f"Plate = {plate}, IFU = {ifu}")

# Set MAPS and LOGCUBE paths.
# For DR17 use defaul, for DR15 use '[...]VOR10-GAU-MILESHC.fits.gz'
logcube = os.path.join(
    logcube_dir,
    str(plate),
    str(ifu),
    f"manga-{plate}-{ifu}-LOGCUBE-VOR10-MILESHC-MASTARSSP.fits.gz",
)
maps = os.path.join(
    maps_dir,
    str(plate),
    str(ifu),
    f"manga-{plate}-{ifu}-MAPS-VOR10-MILESHC-MASTARSSP.fits.gz",
)

# Create output path if it doesn't exist.
direc = os.path.join(output_dir, str(plate), str(ifu))
if not os.path.exists(direc):
    os.makedirs(direc)

# Read in MAPS file as this contains part of the information.
maps_header = fits.open(maps)
unique_bin_number = list(np.unique(maps_header["BINID"].data)[1:])
print(f"Number of bins = {len(unique_bin_number)}")

if bin_number == "all":
    N = f_mp(unique_bin_number)

else:
    f(bin_number)

print("Time to complete: ", (time.time()) - t0)
