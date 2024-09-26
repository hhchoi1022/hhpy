import os
import requests
import json
from collections import OrderedDict
import time
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
from datetime import date
import datetime
from astropy.io import ascii
import time
#%%
class tns_search:
    def __init__(self, 
                 photdata : bool = True,
                 specdata : bool = False):
        self._photkey = int(photdata)
        self._speckey = int(specdata)
        TNS                 = "www.wis-tns.org"
        #TNS                 = "sandbox.wis-tns.org"
        self._url_tns_api         = "https://" + TNS + "/api/get"

        #TNS account info
        self._TNS_BOT_ID          = "117913"
        self._TNS_BOT_NAME        = "HH_bot"
        self._TNS_API_KEY         = "9da2ee1367ba2c229d1c5d6c4aebeca425a5d67c"

        # current working directory
        cwd                 = os.getcwd()

        # directory for downloaded files
        self._download_dir        = os.path.join(cwd, "downloaded_files")

        # external http errors
        self._ext_http_errors     = [403, 500, 503]
        self._err_msg             = ["Forbidden", "Internal Server Error: Something is broken", "Service Unavailable"]

    def _set_bot_tns_marker(self):
        tns_marker = 'tns_marker{"tns_id": "' + str(self._TNS_BOT_ID) + '", "type": "bot", "name": "' + self._TNS_BOT_NAME + '"}'
        return tns_marker

    def _format_to_json(self, source):
        parsed = json.loads(source, object_pairs_hook = OrderedDict)
        result = json.dumps(parsed, indent = 4)
        return result

    def _is_string_json(self, string):
        try:
            json_object = json.loads(string)
        except Exception:
            return False
        return json_object

    def _print_status_code(self, response):
        json_string = self._is_string_json(response.text)
        if json_string != False:
            print ("status code ---> [ " + str(json_string['id_code']) + " - '" + json_string['id_message'] + "' ]\n")
        else:
            status_code = response.status_code
            if status_code == 200:
                status_msg = 'OK'
            elif status_code in self._ext_http_errors:
                status_msg = self._err_msg[self._ext_http_errors.index(status_code)]
            else:
                status_msg = 'Undocumented error'
            print ("status code ---> [ " + str(status_code) + " - '" + status_msg + "' ]\n")

    def _search(self, search_obj):
        search_url = self._url_tns_api + "/search"
        tns_marker = self._set_bot_tns_marker()
        headers = {'User-Agent': tns_marker}
        json_file = OrderedDict(search_obj)
        search_data = {'api_key': self._TNS_API_KEY, 'data': json.dumps(json_file)}
        response = requests.post(search_url, headers = headers, data = search_data)
        return response

    def _get(self, get_obj):
        get_url = self._url_tns_api + "/object"
        tns_marker = self._set_bot_tns_marker()
        headers = {'User-Agent': tns_marker}
        json_file = OrderedDict(get_obj)
        get_data = {'api_key': self._TNS_API_KEY, 'data': json.dumps(json_file)}
        response = requests.post(get_url, headers = headers, data = get_data)
        return response

    def get_obj(self, objname):
        get_obj = [("objname", f"{objname}"), ("objid", ""), ("photometry", f"{self._photkey}"), ("spectra", f"{self._speckey}")]
        get_response = self._get(get_obj)
        json_data = json.loads(get_response.text)
        return json_data
    
    def search_obj(self, ra,dec,radius='25',units='arcmin'):
        search_obj          = [("ra", f"{ra}"), ("dec", f"{dec}"), ("radius", f"{radius}"), ("units", f"{units}"), 
                            ("objname", ""), ("objname_exact_match", 0), ("internal_name", ""), 
                            ("internal_name_exact_match", 0),("objname", ""), ("public_timestamp", "",),
                            ("photometry", f"{self._photkey}"), ("spectra", f"{self._speckey}")]
        search_response = self._search(search_obj)
        json_data = json.loads(search_response.text)
        return json_data
    
    def exist_AT(self, ra, dec, radius='25'):
        # Find transients in the region(ra,dec,radius)
        # Output : pandas table
        ra = str(ra)
        dec = str(dec)
        objdata = self._search_obj(ra,dec,radius=radius)
        if ':' in ra:
            radec = SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
        else:
            radec = SkyCoord(ra,dec,unit=(u.deg,u.deg))
        ATtable = pd.DataFrame(index=range(0,len(objdata['data']['reply'])),columns=['TNSName','DetDate','DetMag','RA','Dec','Separation[arcmin]','reporter'])
        for i, ATinfo in enumerate(objdata['data']['reply']):
            time.sleep(3)
            objname = ATinfo['objname']
            objinfo = self._get_obj(objname)
            repoinfo = objinfo['data']['reply']['reporting_group']['group_name']
            nameinfo = objinfo['data']['reply']['name_prefix']+' '+objinfo['data']['reply']['objname']
            dateinfo = objinfo['data']['reply']['discoverydate']
            maginfo = objinfo['data']['reply']['discoverymag']
            rainfo = objinfo['data']['reply']['ra']
            decinfo = objinfo['data']['reply']['dec']
            radecinfo = SkyCoord(rainfo,decinfo,unit=(u.hourangle,u.deg))
            sepinfo = SkyCoord.separation(radecinfo,radec).arcminute
            
            ATtable['TNSName'][i] = nameinfo
            ATtable['DetDate'][i] = dateinfo
            ATtable['DetMag'][i] = maginfo
            ATtable['RA'][i] = rainfo
            ATtable['Dec'][i] = decinfo
            ATtable['Separation[arcmin]'][i] = sepinfo
            ATtable['reporter'][i] = repoinfo
        return ATtable