import os
import sys
from astropy.io.fits import Header, PrimaryHDU
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
        prihdr = Header()
        prihdr["FILE"] = os.path.basename(output_file)
        prihdr["MODELS"] = config.model_key
        prihdr["FITTER"] = "FIREFLY"
        prihdr["AGEMIN"] = str(config.age_limits[0])
        prihdr["AGEMAX"] = str(config.age_limits[1])
        prihdr["ZMIN"] = str(config.z_limits[0])
        prihdr["ZMAX"] = str(config.z_limits[1])
        prihdr["redshift"] = config.redshift
        prihdr["HIERARCH age_universe"] = np.round(config.age_of_universe(), 3)
        prihdu = PrimaryHDU(header=prihdr)
        return [prihdu]