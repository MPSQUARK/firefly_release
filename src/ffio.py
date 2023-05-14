import os
import sys
import numpy as np
from astropy.io import fits
from src.config import Config

class IO():
    
    OVERWRITEWARNING = "\nWarning: This object has already been processed, \
        the file will be over-witten.\n** Do you want to continue? (Y/N)"
    
    @staticmethod
    def ensure_dir(path :str) -> None:
        '''Create directory if it does not exist'''
        if not os.path.isdir(path):
            os.mkdir(path)
        
        
    @staticmethod
    def warn_overwrite(file :str) -> None:
        '''Warn if file exists and ask for permission to overwrite'''
        if not os.path.isfile(file):
            return 
        if input(IO.OVERWRITEWARNING).upper() == "N":
            sys.exit()
        os.remove(file)
        
    @staticmethod
    def create_fits_tables(output_file, config : Config):
        prihdr = fits.Header()
        prihdr["FILE"] = os.path.basename(output_file)
        prihdr["MODELS"] = config.model_key
        prihdr["FITTER"] = "FIREFLY"
        prihdr["AGEMIN"] = str(config.age_limits[0])
        prihdr["AGEMAX"] = str(config.age_limits[1])
        prihdr["ZMIN"] = str(config.z_limits[0])
        prihdr["ZMAX"] = str(config.z_limits[1])
        prihdr["redshift"] = config.redshift
        prihdr["HIERARCH age_universe"] = np.round(config.age_of_universe(), 3)
        prihdu = fits.PrimaryHDU(header=prihdr)
        return [prihdu]
    
    @staticmethod
    def read_config_and_input() -> tuple[Config, dict]:
        config = Config()
        data = IO.read_input(config)
        return (config, data)
    
    @staticmethod
    def read_input(config : Config) -> dict:
        params = sys.argv[1:]
        if len(params) == 0:
            config.create("config/default.yaml")
            return IO.read_dat_input(config.file)
        elif params[1] == "sdss":
            config.create("config/sdss.yaml")
            return IO.read_sdss_input(config)
        raise Exception("Unknown Config")
            
    @staticmethod
    def read_dat_input(filepath) -> dict:
        file_contents = np.loadtxt(filepath, unpack=True)
        return {
            "wavelength": file_contents[0, :], 
            "flux": file_contents[1, :], 
            "error": file_contents[2, :]
            }
        
    @staticmethod
    def read_sdss_input(config : Config) -> dict:
        hdul = fits.open(config.file)
        
        # update config
        config.redshift = hdul[2].data["Z"][0]
        config.ra = hdul[0].header["RA"]
        config.dec = hdul[0].header["DEC"]
        config.vdisp = hdul[2].data["VDISP"][0]
        
        # read data
        data = {
            "wavelength": 10 ** hdul[1].data["loglam"],
            "flux": hdul[1].data["flux"],
            "error": hdul[1].data["ivar"] ** (-0.5)
            }
        
        return data
        