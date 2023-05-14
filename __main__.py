import os
import time
import traceback
import numpy as np
from astropy.io import fits
from src.setup import Setup
from src.models import StellarPopulationModel
from src.io import IO

os.environ["FF_DIR"] = os.getcwd()
os.environ["STELLARPOPMODELS_DIR"] = os.path.join(
    os.environ["FF_DIR"], "stellar_population_models"
)

time_start = time.time()

config, data = IO.read_config_and_input()

restframe_wavelength = data.get('wavelength') / (1 + config.redshift)

# instrumental resolution
r_instrument = np.full(len(data.get('wavelength')), 2000)

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
    spec = Setup(config).openSingleSpectrum(data.get('wavelength'),data.get('flux'),data.get('error'),r_instrument)
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

    print(f"\nDone... total time: {time.time() - time_start:.2f} seconds.\n")
