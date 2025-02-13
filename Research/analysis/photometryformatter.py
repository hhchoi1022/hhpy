
#%%
from astropy.io import ascii
from astropy.table import Table
import numpy as np
#%%
class PhotometryFormatter:
    
    def __init__(self,
                 filename : str, format : str = 'fixed_width'):
        self.filename = filename
        self.original_data = self._load_table(format = format)
        self.formatted_data = None
        self.format_history = list()

    def formatting(self):
        original_data = self.original_data.copy()
        # Set/Modify the columns to the standard format
        formatted_tbl = Table()
        default_data = {key: np.array([value[0]] * len(original_data), dtype = value[1]) for key, value in self.default_format.items()}
        # for key, value in self.default_format.items():
        #     formatted_tbl[key] = [value[0]] * len(original_data)
        #     formatted_tbl[key] = value[1]
        formatted_tbl = Table(default_data)
                    
        for key in original_data.keys():
            canonical_key = self._normalize_required_keys(key)
            if canonical_key is not None:
                formatted_tbl[canonical_key] = original_data[key]
            
        if "status" in original_data.keys():
            formatted_tbl['detected'] = [True if status == 'detected' else False for status in original_data['status']]
        self.formatted_data = formatted_tbl
    
    def change_value(self, key : str, value):
        if not self.formatted_data:
            raise ValueError('No formatted data. Run PhotometryFormatter.formatting()')
        if not key in self.formatted_data.keys():
            raise ValueError(f'Key "{key}" not in formatted data')
        self.formatted_data[key] = value
    
    def change_magsys_depending_on_filter(self,
                                          filtername_list : list = ['UBVRI', 'ugriz'],
                                          magsys_list : list = ['Vega', 'AB']):
        if isinstance(filtername_list,str):
            filtername_list = [filtername_list]
        if isinstance(magsys_list,str):
            magsys_list = [magsys_list]
        if not self.formatted_data:
            raise ValueError('No formatted data. Run PhotometryFormatter.formatting()')
        
        filter_to_magsys = {}
        for filter_names, maysys_names in zip(filtername_list, magsys_list):
            for filter_name in filter_names:
                filter_to_magsys[filter_name] = maysys_names

        matched_magsys = [filter_to_magsys.get(f, None) for f in self.formatted_data['filter']]
        self.formatted_data['magsys'] = matched_magsys
        
    def convert_magsys(self, to_magsys :str = 'AB'):
        from convert_AB_Vega import ABVegaMagnitude
        if not self.formatted_data:
            raise ValueError('No formatted data. Run PhotometryFormatter.formatting()')
        if not 'magsys' in self.formatted_data.keys():
            raise ValueError('No magsys key in the formatted data')
        
        maglist = []
        magsyslist = []
        for value in self.formatted_data:
            if bool(value['magsys']) & bool(value['mag']) & bool(value['filter']):
                try:
                    mag = ABVegaMagnitude(magnitude = value['mag'], magsys = value['magsys'], filter_ = value['filter'])
                    if to_magsys == 'AB':
                        maglist.append(mag.AB)
                        magsyslist.append('AB')
                    elif to_magsys == 'Vega':
                        maglist.append(mag.Vega)
                        magsyslist.append('Vega')
                    else:
                        raise ValueError(f'Not defiened magsys: {to_magsys}')
                except:
                    maglist.append(value['mag'])
                    magsyslist.append(value['magsys'])
            else:
                if bool(value['mag']):
                    maglist.append(value['mag'])
                else:
                    maglist.append(np.nan)
                magsyslist.append('Undefined')
        self.formatted_data['mag'] = maglist    
        self.formatted_data['magsys'] = magsyslist        
    
    def write(self, filename : str, format : str = 'ascii.fixed_width'):
        if not self.formatted_data:
            raise ValueError('No formatted data. Run PhotometryFormatter.formatting()')
        self.formatted_data.write(filename, format = format, overwrite = True)
        
    def _normalize_required_keys(self, key: str):

        # Iterate through the dictionary to find a match
        for canonical_key, variants in self.required_key_variants.items():
            if key.lower() in variants:
                print(canonical_key)
                return canonical_key
        return None

    def _load_table(self, format : str = 'fixed_width'):
        tbl = ascii.read(self.filename, format = format)
        return tbl

    @property
    def default_format(self):
        default_format = dict()
        default_format['obsdate'] = (np.nan,'float64')
        default_format['mag'] = (np.nan,'float64')
        default_format['e_mag'] = (np.nan,'float64')
        default_format['magsys'] = (None,'str')
        default_format['filter'] = (None,'str')
        default_format['depth_5sig'] = (np.nan,'float64')
        default_format['zp'] = (np.nan,'float64')
        default_format['observatory'] = (None,'str')
        default_format['detected'] = (None,'bool')
        default_format['A_filter'] = (0.0,'float64')
        return default_format
    
    @property
    def required_key_variants(self):
        # Define key variants, if a word is duplicated in the same variant, posit the word with the highest priority first
        required_key_variants_lower = {
            'obsdate': ['obsdate', 'MJD'],
            'mag': ['mag', 'magnitude'],
            'e_mag': ['e_mag', 'magnitude error'],
            'magsys': ['magsys', 'magnitude system'],
            'filter': ['filter', 'filter_', 'band'],
            'depth_5sig': ['depth_5sig', '5sigma depth', 'depth'],
            'zp': ['zp, zero point', 'zeropoint'],
            'observatory': ['observatory'],
            'detected': ['detected'],
            'A_filter': ['magsys']
        }
        # Sort each list in the dictionary by string length (descending order)
        sorted_required_key_variants = {
            key: sorted(variants, key=len, reverse=True)
            for key, variants in required_key_variants_lower.items()
        }
        return sorted_required_key_variants
# %%
P = PhotometryFormatter('./data/SN2021aefx/Ashall2022_all.dat')

# %%
P.formatting()
P.change_magsys_depending_on_filter(['UBVRI','ugriz',[ 'B_S', 'UVM2', 'UVW1', 'UVW2', 'U_S', 'V_S']], ['Vega','AB','AB'])
P.convert_magsys(to_magsys = 'AB')
P.write('./data/SN2021aefx/Ashall2022.dat')
# %%
