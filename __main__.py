import sys
import os
import time
import traceback
import numpy as np
from astropy.io import fits
from src.setup import Setup
from src.models import StellarPopulationModel
from src.config import Config
from src.io import IO

params = sys.argv[1:]

os.environ["STELLARPOPMODELS_DIR"] = os.path.join(
    os.getcwd(), "stellar_population_models"
)

t0 = time.time()

config = Config()

if len(params) == 0:
    config.create("config/default.yaml")
    data = np.loadtxt(config.file, unpack=True)
    wavelength = data[0, :]
    flux = data[1, :]
    error = data[2, :]
elif params[1] == "sdss":
    config.create("config/sdss.yaml")
    hdul = fits.open(config.file)
    # update config
    config.redshift = hdul[2].data["Z"][0]
    config.ra = hdul[0].header["RA"]
    config.dec = hdul[0].header["DEC"]
    config.vdisp = hdul[2].data["VDISP"][0]
    
    wavelength = 10 ** hdul[1].data["loglam"]
    flux = hdul[1].data["flux"]
    error = hdul[1].data["ivar"] ** (-0.5)
else: 
    raise Exception("Unknown Config")

restframe_wavelength = wavelength / (1 + config.redshift)
# instrumental resolution
r_instrument = np.zeros(len(wavelength))
for wi, w in enumerate(wavelength):
    r_instrument[wi] = 2000    

config.verify()

print("\nStarting firefly ...")

# set output folder and output filename in firefly directory
# and write output file
outputFolder = os.path.join(os.getcwd(), "output")
output_file_name, _ = os.path.splitext(os.path.basename(config.file))
output_file = os.path.join(outputFolder, f"spFly-{output_file_name}.fits")

IO.warn_overwrite(output_file)
IO.ensure_dir(outputFolder)

print(f"\nOutput file: {output_file}\n")

tables = IO.create_fits_tables(output_file, config)

did_not_converge = 0.0
try:
    # define input object to pass data on to firefly modules and initiate run
    spec = Setup(config).openSingleSpectrum(wavelength,flux,error,r_instrument)
    
    # prepare model templates
    model = StellarPopulationModel(spec,output_file,config)

    # initiate fit
    model.fit_models_to_data()
    tables.append(model.tbhdu)
    
except ValueError:
    tables.append(model.create_dummy_hdu())
    did_not_converge += 1
    print("did not converge")
except Exception:
    traceback.print_exc()
finally:
    if did_not_converge < 1:
        complete_hdus = fits.HDUList(tables)
        if os.path.isfile(output_file):
            os.remove(output_file)
        complete_hdus.writeto(output_file)

    print("\nDone... total time:", int(time.time() - t0), "seconds.\n")
