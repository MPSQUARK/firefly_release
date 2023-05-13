import sys
import os
import time
import traceback
import numpy as np
from astropy.io.fits import Header, PrimaryHDU, HDUList
import astropy.cosmology as co
from src.setup import Setup
from src.models import StellarPopulationModel
from src.config import Config

#sys.path.append(os.path.join(os.getcwd(), "python"))
os.environ["FF_DIR"] = os.getcwd()
os.environ["STELLARPOPMODELS_DIR"] = os.path.join(os.environ["FF_DIR"], "stellar_population_models")

t0=time.time()

config = Config().create('config/default.yaml').verify()

data = np.loadtxt(config.file, unpack=True)

wavelength = data[0,:]
flux = data[1,:]
error = data[2,:]
restframe_wavelength = wavelength/(1+config.redshift)
#instrumental resolution
r_instrument = np.zeros(len(wavelength))
for wi,w in enumerate(wavelength):
    r_instrument[wi] = 2000


print('\nStarting firefly ...')

Z_min=config.z_limits[0]
Z_max=config.z_limits[1]

#set output folder and output filename in firefly directory
#and write output file
outputFolder = os.path.join( os.environ['FF_DIR'], 'output')
output_file_name, _ = os.path.splitext(os.path.basename(config.file))
output_file = os.path.join(outputFolder, f'spFly-{output_file_name}.fits')


if os.path.isfile(output_file):
    print('\nWarning: This object has already been processed, the file will be over-witten.')
    answer = input('** Do you want to continue? (Y/N)')
    if answer.upper() == 'N':
        sys.exit()
    os.remove(output_file)
if os.path.isdir(outputFolder) is False:
    os.mkdir(outputFolder)

print(f'\nOutput file: {output_file}\n')

prihdr = Header()
prihdr['FILE']     = os.path.basename(output_file)
prihdr['MODELS']   = config.model_key
prihdr['FITTER']   = "FIREFLY"
prihdr['AGEMIN']   = str(config.age_limits[0])
prihdr['AGEMAX']   = str(config.age_limits[1])
prihdr['ZMIN']     = str(Z_min)
prihdr['ZMAX']	   = str(Z_max)
prihdr['redshift'] = config.redshift
prihdr['HIERARCH age_universe']	= np.round(config.age_of_universe(),3)
prihdu = PrimaryHDU(header=prihdr)
tables = [prihdu]

#define input object to pass data on to firefly modules and initiate run
spec= Setup(config.file,config.milky_way_reddening,config.hpf_mode,config.n_angstrom_masked)

spec.openSingleSpectrum(wavelength, flux, error, config.redshift, config.ra, config.dec, config.vdisp, config.emlines, r_instrument)

did_not_converge = 0.

try:
    #prepare model templates
    model = StellarPopulationModel(
        spec,
        output_file,
        config.cosmo,
        models = config.model_key,
        model_libs = config.model_lib, 
        imfs = config.imfs,
        age_limits = config.age_limits,
        downgrade_models = True,
        data_wave_medium = config.data_wave_medium,
        Z_limits = config.z_limits,
        suffix=config.suffix,
        use_downgraded_models = False,
        dust_law=config.dust_law,
        max_ebv=config.max_ebv,
        num_dust_vals=config.num_dust_vals,
        dust_smoothing_length=config.dust_smoothing_length,
        max_iterations=config.max_itterations,
        pdf_sampling=config.pdf_sampling,
        flux_units=config.flux_units,
        )

    #initiate fit
    model.fit_models_to_data()
    tables.append(model.tbhdu)
except ValueError:
    tables.append( model.create_dummy_hdu() )
    did_not_converge +=1
    print('did not converge')
except Exception:
    traceback.print_exc()
finally:
    if did_not_converge < 1 :
        complete_hdus = HDUList(tables)
        if os.path.isfile(output_file):
            os.remove(output_file)
        complete_hdus.writeto(output_file)

    print("\nDone... total time:", int(time.time()-t0) ,"seconds.\n")
