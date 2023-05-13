'''Configuration module for firefly'''

from os import path
from typing import Union
import yaml
from astropy.cosmology import FLRW, Planck15

class Config():
    '''Configuration class for firefly'''

    _cosmo : FLRW
    _file : str
    _suffix : str
    _redshift : float
    _ra : float
    _dec : float
    _vdisp : float
    _n_angstrom_masked : int
    _emlines : list[str]
    _model_key : str
    _model_lib : list[str]
    _imfs : list[str]
    _age_limits : tuple[float, float]
    _z_limits : tuple[float, float]
    _data_wave_medium : str
    _flux_units : str
    _write_results : bool
    _milky_way_reddening : bool
    _hpf_mode : str
    _dust_law : str
    _max_ebv : float
    _num_dust_vals : int
    _dust_smoothing_length : int
    _max_itterations : int
    _pdf_sampling : int

    @property
    def cosmo(self) -> FLRW:
        '''Cosmology object from astropy.cosmology'''
        return self._cosmo
    @cosmo.setter
    def cosmo(self, value: FLRW) -> None:
        self._cosmo = value

    @property
    def file(self) -> str:
        '''Path to the input file'''
        return self._file
    @file.setter
    def file(self, value: str) -> None:
        self._file = value

    @property
    def suffix(self) -> str:
        '''Suffix to append to the output file name'''
        return self._suffix
    @suffix.setter
    def suffix(self, value: str) -> None:
        self._suffix = value

    @property
    def redshift(self) -> float:
        '''Redshift of the object'''
        return self._redshift
    @redshift.setter
    def redshift(self, value: float) -> None:
        self._redshift = value

    @property
    def ra(self) -> float:
        '''Right ascension of the object'''
        return self._ra
    @ra.setter
    def ra(self, value: float) -> None:
        self._ra = value

    @property
    def dec(self) -> float:
        '''Declination of the object'''
        return self._dec
    @dec.setter
    def dec(self, value: float) -> None:
        self._dec = value

    @property
    def vdisp(self) -> float:
        '''Velocity dispersion of the object'''
        return self._vdisp
    @vdisp.setter
    def vdisp(self, value: float) -> None:
        self._vdisp = value

    @property
    def n_angstrom_masked(self) -> int:
        '''Number of angstroms to mask on either side of the emission lines'''
        return self._n_angstrom_masked
    @n_angstrom_masked.setter
    def n_angstrom_masked(self, value: int) -> None:
        self._n_angstrom_masked = value

    @property
    def emlines(self) -> list[str]:
        '''List of emission lines to mask'''
        return self._emlines
    @emlines.setter
    def emlines(self, value: list[str]) -> None:
        self._emlines = value

    @property
    def model_key(self) -> str:
        '''Which model to use'''
        return self._model_key
    @model_key.setter
    def model_key(self, value: str) -> None:
        self._model_key = value

    @property
    def model_lib(self) -> list[str]:
        '''Which model library to use'''
        return self._model_lib
    @model_lib.setter
    def model_lib(self, value: list[str]) -> None:
        self._model_lib = value

    @property
    def imfs(self) -> list[str]:
        '''List of IMFs to use'''
        return self._imfs
    @imfs.setter
    def imfs(self, value: list[str]) -> None:
        self._imfs = value

    @property
    def age_limits(self) -> tuple[float, float]:
        '''Age limits to use'''
        return self._age_limits
    @age_limits.setter
    def age_limits(self, value: tuple[float, float]) -> None:
        self._age_limits = value

    @property
    def z_limits(self) -> tuple[float, float]:
        '''Redshift limits to use'''
        return self._z_limits
    @z_limits.setter
    def z_limits(self, value: tuple[float, float]) -> None:
        self._z_limits = value

    @property
    def data_wave_medium(self) -> str:
        '''Which medium the data is in'''
        return self._data_wave_medium
    @data_wave_medium.setter
    def data_wave_medium(self, value: str) -> None:
        self._data_wave_medium = value

    @property
    def flux_units(self) -> str:
        '''Units of the flux'''
        return self._flux_units
    @flux_units.setter
    def flux_units(self, value: str) -> None:
        self._flux_units = value

    @property
    def write_results(self) -> bool:
        '''Whether or not to write the results to a file'''
        return self._write_results
    @write_results.setter
    def write_results(self, value: bool) -> None:
        self._write_results = value

    @property
    def milky_way_reddening(self) -> bool:
        '''Whether or not to apply Milky Way reddening'''
        return self._milky_way_reddening
    @milky_way_reddening.setter
    def milky_way_reddening(self, value: bool) -> None:
        self._milky_way_reddening = value

    @property
    def hpf_mode(self) -> str:
        '''Which mode to use for the high pass filter'''
        return self._hpf_mode
    @hpf_mode.setter
    def hpf_mode(self, value: str) -> None:
        self._hpf_mode = value

    @property
    def dust_law(self) -> str:
        '''Which dust law to use'''
        return self._dust_law
    @dust_law.setter
    def dust_law(self, value: str) -> None:
        self._dust_law = value

    @property
    def max_ebv(self) -> float:
        '''Maximum E(B-V) to use'''
        return self._max_ebv
    @max_ebv.setter
    def max_ebv(self, value: float) -> None:
        self._max_ebv = value

    @property
    def num_dust_vals(self) -> int:
        '''Number of dust values to use'''
        return self._num_dust_vals
    @num_dust_vals.setter
    def num_dust_vals(self, value: int) -> None:
        self._num_dust_vals = value

    @property
    def dust_smoothing_length(self) -> int:
        '''Smoothing length for the dust'''
        return self._dust_smoothing_length
    @dust_smoothing_length.setter
    def dust_smoothing_length(self, value: int) -> None:
        self._dust_smoothing_length = value

    @property
    def max_itterations(self) -> int:
        '''Maximum number of itterations to use'''
        return self._max_itterations
    @max_itterations.setter
    def max_itterations(self, value: int) -> None:
        self._max_itterations = value

    @property
    def pdf_sampling(self) -> int:
        '''Number of samples to use for the PDF'''
        return self._pdf_sampling
    @pdf_sampling.setter
    def pdf_sampling(self, value: int) -> None:
        self._pdf_sampling = value

    def age_of_universe(self) -> float:
        '''Age of the universe
        Notes:
        ---
        cosmo & redshift must be set
        '''
        return self.cosmo.age(self.redshift).value

    def age_limit_to_value(self, value : Union[str,float]):
        '''
        Convert the age limit to a value
        Notes:
        ---
        cosmo & redshift must be set
        '''
        if isinstance(value, str):
            if value.lower() == 'aou':
                return self.age_of_universe()
            raise ValueError(f'Unknown age limit: {value}')
        return value

    def comology_model(self, value : str) -> FLRW:
        '''Get the cosmology model'''
        if value.lower() == 'planck15':
            return Planck15
        raise ValueError(f'Unknown cosmology model: {value}')

    def create(self, config_file_path : str):
        '''Create the config object from a yaml file'''
        assert path.exists(config_file_path), f'Config file {config_file_path} does not exist'

        with open(config_file_path, 'r', encoding='utf-8') as config_file:
            config = yaml.safe_load(config_file)

        self.cosmo = self.comology_model(config.get('Cosmo'))
        self.file = config.get('File')
        self.suffix = config.get('Suffix')
        self.redshift = config.get('Redshift')
        self.ra = config.get('Right_Ascension')
        self.dec = config.get('Declination')
        self.vdisp = config.get('Velocity_Dispersion')
        self.n_angstrom_masked = config.get('N_Angstrom_masked', 20)
        self.model_key = config.get('Model_Key')
        self.model_lib = config.get('Model_Library')
        self.imfs = config.get('IMFs')
        self.age_limits = (
            self.age_limit_to_value(config.get('Age_Min')),
            self.age_limit_to_value(config.get('Age_Max'))
            )
        self.z_limits = (
            config.get('Z_Min'),
            config.get('Z_Max')
        )
        self.data_wave_medium = config.get('Wave_Medium')
        self.flux_units = config.get('Flux_Units')
        self.write_results = config.get('Write_Results')
        self.milky_way_reddening = config.get('Milkyway_Reddening', True)
        self.hpf_mode = config.get('HPF_Mode', 'on')
        self.dust_law = config.get('Dust_Law')
        self.max_ebv = config.get('Max_EBV')
        self.num_dust_vals = config.get('Num_Dust_Vals')
        self.dust_smoothing_length = config.get('Dust_Smoothing_Length')
        self.max_itterations = config.get('Max_Itterations')
        self.pdf_sampling = config.get('PDF_Sampling')
        self.emlines = config.get('Emission_Lines')
        
        return self

    def verify(self):
        # Check if all the required values are set
        assert self.cosmo is not None, 'Cosmo not set'
        assert self.file is not None, 'File not set'
        # OPTIONAL: assert self.suffix is not None, 'Suffix not set'
        assert self.redshift is not None, 'Redshift not set'
        assert self.ra is not None, 'Right Ascension not set'
        assert self.dec is not None, 'Declination not set'
        assert self.vdisp is not None, 'Velocity Dispersion not set'
        assert self.n_angstrom_masked is not None, 'N Angstrom Masked not set'
        assert self.model_key is not None, 'Model Key not set'
        assert self.model_lib is not None, 'Model Library not set'
        assert self.imfs is not None, 'IMFs not set'
        assert self.age_limits is not None, 'Age Limits not set'
        assert self.z_limits is not None, 'Z Limits not set'
        assert self.data_wave_medium is not None, 'Wave Medium not set'
        assert self.flux_units is not None, 'Flux Units not set'
        assert self.write_results is not None, 'Write Results not set'
        assert self.milky_way_reddening is not None, 'Milkyway Reddening not set'
        assert self.hpf_mode is not None, 'HPF Mode not set'
        assert self.dust_law is not None, 'Dust Law not set'
        assert self.max_ebv is not None, 'Max E(B-V) not set'
        assert self.num_dust_vals is not None, 'Num Dust Vals not set'
        assert self.dust_smoothing_length is not None, 'Dust Smoothing Length not set'
        assert self.max_itterations is not None, 'Max Itterations not set'
        assert self.pdf_sampling is not None, 'PDF Sampling not set'
        assert self.emlines is not None, 'Emission Lines not set'
        
        # Check if the values are valid
        # assert self.file != '', 'File path is empty'
        # assert path.exists(self.file), f'File {self.file} does not exist'
        # assert self.model_key in {'m11', 'MaStar'}, 'Invalid model key'
        # if self.model_key == 'm11':
        #     assert self.model_lib in {'MILES', 'STELIB', 'ELODIE', 'MARCS'}, 'Invalid M11 model library'
        # elif self.model_key == 'MaStar':
        #     assert self.model_lib in {'gold'}, 'Invalid MaStar model library'
        # for imf in self.imfs:
        #     assert imf in {'kr', 'ss'}, 'Invalid IMF'
        # assert self.age_limits[0] < self.age_limits[1], 'Min age is greater than max age'
        # assert self.data_wave_medium in {'vacuum', 'air'}, 'Invalid wave medium'
        # assert self.hpf_mode in {'on', 'hpf_only'}
        # assert self.dust_law in {'calzetti', 'allen', 'prevot'}, 'Invalid dust law'
        
        return self
    