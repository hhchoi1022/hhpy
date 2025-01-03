#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Avgust 2021 by Nikola Knezevic ASTRO DATA
Edited on January 2022 by Hyeonho Choi

Developed and tested on:

- Linux 20.04 LTS
- Python 3.6 (Spyder 4)
"""

#%%
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
#%%#----------------------------------------------------------------------------------


TNS                 = "www.wis-tns.org"
#TNS                 = "sandbox.wis-tns.org"
url_tns_api         = "https://" + TNS + "/api/get"

#TNS account info
TNS_BOT_ID          = "117913"
TNS_BOT_NAME        = "HH_bot"
TNS_API_KEY         = "9da2ee1367ba2c229d1c5d6c4aebeca425a5d67c"

# current working directory
cwd                 = os.getcwd()

# directory for downloaded files
download_dir        = os.path.join(cwd, "downloaded_files")

# external http errors
ext_http_errors     = [403, 500, 503]
err_msg             = ["Forbidden", "Internal Server Error: Something is broken", "Service Unavailable"]


#%%#----------------------------------------------------------------------------------
def set_bot_tns_marker():
    tns_marker = 'tns_marker{"tns_id": "' + str(TNS_BOT_ID) + '", "type": "bot", "name": "' + TNS_BOT_NAME + '"}'
    return tns_marker

def is_string_json(string):
    try:
        json_object = json.loads(string)
    except Exception:
        return False
    return json_object

def print_status_code(response):
    json_string = is_string_json(response.text)
    if json_string != False:
        print ("status code ---> [ " + str(json_string['id_code']) + " - '" + json_string['id_message'] + "' ]\n")
    else:
        status_code = response.status_code
        if status_code == 200:
            status_msg = 'OK'
        elif status_code in ext_http_errors:
            status_msg = err_msg[ext_http_errors.index(status_code)]
        else:
            status_msg = 'Undocumented error'
        print ("status code ---> [ " + str(status_code) + " - '" + status_msg + "' ]\n")

def search(search_obj):
    search_url = url_tns_api + "/search"
    tns_marker = set_bot_tns_marker()
    headers = {'User-Agent': tns_marker}
    json_file = OrderedDict(search_obj)
    search_data = {'api_key': TNS_API_KEY, 'data': json.dumps(json_file)}
    response = requests.post(search_url, headers = headers, data = search_data)
    return response

def get(get_obj):
    get_url = url_tns_api + "/object"
    tns_marker = set_bot_tns_marker()
    headers = {'User-Agent': tns_marker}
    json_file = OrderedDict(get_obj)
    get_data = {'api_key': TNS_API_KEY, 'data': json.dumps(json_file)}
    response = requests.post(get_url, headers = headers, data = get_data)
    return response

'''
def get_file():
    filename = os.path.basename(file_tns_url)
    tns_marker = set_bot_tns_marker()
    headers = {'User-Agent': tns_marker}
    api_data = {'api_key': TNS_API_KEY}
    print ("Downloading file '" + filename + "' from the TNS...\n")
    response = requests.post(file_tns_url, headers = headers, data = api_data, stream = True)    
    print_status_code(response)
    path = os.path.join(download_dir, filename)
    if response.status_code == 200:
        with open(path, 'wb') as f:
            for chunk in response:
                f.write(chunk)
        print ("File was successfully downloaded.\n")
    else:
        print ("File was not downloaded.\n")
'''

def print_response(response, json_file, counter):
    response_code = str(response.status_code) if json_file == False else str(json_file['id_code'])
    stats = 'Test #' + str(counter) + '| return code: ' + response_code + \
            ' | Total Rate-Limit: ' + str(response.headers.get('x-rate-limit-limit')) + \
            ' | Remaining: ' + str(response.headers.get('x-rate-limit-remaining')) + \
            ' | Reset: ' + str(response.headers.get('x-rate-limit-reset'))
    if(response.headers.get('x-cone-rate-limit-limit') != None):
        stats += ' || Cone Rate-Limit: ' + str(response.headers.get('x-cone-rate-limit-limit')) + \
                 ' | Cone Remaining: ' + str(response.headers.get('x-cone-rate-limit-remaining')) + \
                 ' | Cone Reset: ' + str(response.headers.get('x-cone-rate-limit-reset'))
    print (stats)

def get_reset_time(response):
    # If any of the '...-remaining' values is zero, return the reset time
    for name in response.headers:
        value = response.headers.get(name)
        if name.endswith('-remaining') and value == '0':
            return int(response.headers.get(name.replace('remaining', 'reset')))
    return None

def rate_limit_handling():
    counter = 0
    while True:
        counter = counter + 1
        response = search()
        json_file = is_string_json(response.text)
        print_response(response, json_file, counter)
        # Checking if rate-limit reached (...-remaining = 0)
        reset = get_reset_time(response)
        # A general verification if not some error 
        if (response.status_code == 200):
            if reset != None:
                # Sleeping for reset + 1 sec
                print("Sleep for " + str(reset + 1) + " sec") 
                time.sleep(reset + 1)
        	    # Can continue to submit requests...
                print ("Continue to submit requests...")
                for i in range(3):
                    counter = counter + 1
                    response = search()
                    json_file = is_string_json(response.text)
                    print_response(response, json_file, counter)
                print ("etc...\n") 
                break
        else:
            print_status_code(response)       
            break
#%%#----------------------------------------------------------------------------------

def search_obj(ra,dec,radius='25',units='arcmin'):
    search_obj          = [("ra", f"{ra}"), ("dec", f"{dec}"), ("radius", "25"), ("units", "arcmin"), 
                           ("objname", ""), ("objname_exact_match", 0), ("internal_name", ""), 
                           ("internal_name_exact_match", 0),("objname", ""), ("public_timestamp", "")]
    search_response = search(search_obj)
    json_data = json.loads(search_response.text)
    return json_data





def exist_AT(ra,dec,radius='25'):
    # Find transients in the region(ra,dec,radius)
    # Output : pandas table
    ra = str(ra)
    dec = str(dec)
    objdata = search_obj(ra,dec,radius=radius)
    if ':' in ra:
        radec = SkyCoord(ra,dec,unit=(u.hourangle,u.deg))
    else:
        radec = SkyCoord(ra,dec,unit=(u.deg,u.deg))
    ATtable = pd.DataFrame(index=range(0,len(objdata['data']['reply'])),columns=['TNSName','DetDate','DetMag','RA','Dec','Separation[arcmin]','reporter'])
    for i, ATinfo in enumerate(objdata['data']['reply']):
        time.sleep(3)
        objname = ATinfo['objname']
        objinfo = get_obj(objname)
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

#%%#
os.makedirs('/Users/hhchoi1022/Gitrepo/makereference/TNS')
os.chdir('/Users/hhchoi1022/Gitrepo/makereference/TNS')

for name, ra, dec in alltargetlist['obj','ra','dec']:
    print(f'Searching Transient in {name}')
    ATtable = exist_AT(ra,dec,radius = 140)
    ATtable.to_csv(f'TNS_{name}.csv')
    time.sleep(3)



"""
# EXAMPLE 3 (get file from TNS)
file_tns_url        = "https://" + TNS + "/system/files/uploaded/"\
                      "Padova-Asiago/tns_2017A_2457777.69_Ekar_AFOSC_Padova-Asiago.txt"
get_file()
"""

"""
# EXAMPLE 4 (test rate-limit search)
search_obj          = [("ra", ""), ("dec", ""), ("radius", ""), ("units", ""), 
                       ("objname", "2021rak"), ("objname_exact_match", 0), ("internal_name", ""), 
                       ("internal_name_exact_match", 0), ("objid", ""), ("public_timestamp", "")]
rate_limit_handling()
"""

"""
# EXAMPLE 5 (test rate-limit cone search)
search_obj          = [("ra", "15:57:28"), ("dec", "+30:03:39"), ("radius", "5"), ("units", "arcsec"), 
                       ("objname", ""), ("objname_exact_match", 0), ("internal_name", ""), 
                       ("internal_name_exact_match", 0), ("objid", ""), ("public_timestamp", "")]
rate_limit_handling()
"""

#----------------------------------------------------------------------------------

#%%
import os
import json
import time
import requests
import pandas as pd
from collections import OrderedDict
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u


class TNSQuerier:
    """
    A class for querying the TNS server and saving results to CSV files.
    """

    def __init__(self, api_key, bot_id, bot_name, download_dir):
        """
        Initialize the TNSQuerier class.

        Parameters:
            api_key (str): The API key for TNS.
            bot_id (str): The bot ID.
            bot_name (str): The bot name.
            download_dir (str): Directory to save the downloaded files.
        """
        self.api_key = api_key
        self.bot_id = bot_id
        self.bot_name = bot_name
        self.base_url = "https://www.wis-tns.org/api/get"
        self.download_dir = download_dir
        self.tns_marker = self._set_bot_tns_marker()
        os.makedirs(download_dir, exist_ok=True)

    def _set_bot_tns_marker(self):
        """Set the bot marker for TNS API."""
        return f'tns_marker{{"tns_id": "{self.bot_id}", "type": "bot", "name": "{self.bot_name}"}}'

    def search_by_objname(self, objname : str):
        objinfo_dict = [("objname", f"{objname}"), ("objid", ""), ("photometry", "0"), ("spectra", "0")]
        get_url = self.base_url + "/object"
        headers = {'User-Agent': self.tns_marker}
        query_json = OrderedDict(objinfo_dict)
        search_data = {'api_key': self.api_key, 'data': json.dumps(query_json)}
        response = requests.post(get_url, headers = headers, data = search_data)
        json_data = json.loads(response.text)
        return json_data

    def search_by_coordinates(self, 
                              ra : float, # float 
                              dec : float, # float
                              radius : str = 30 # arcmin
                              ):
        """
        Search for objects in a region around the given RA and Dec.

        Parameters:
            ra (float): Right Ascension in degrees.
            dec (float): Declination in degrees.
            radius (str): Search radius (default: 25 arcmin).
            units (str): Units of the radius (default: arcmin).

        Returns:
            dict: JSON response from the TNS server.
        """
        objinfo_dict = [
            ("ra", f"{ra}"), ("dec", f"{dec}"), ("radius", radius), ("units", 'arcmin'),
            ("objname", ""), ("objname_exact_match", 0), ("internal_name", ""),
            ("internal_name_exact_match", 0), ("objname", ""), ("public_timestamp", "")
        ]
        search_url = f"{self.base_url}/search"
        headers = {"User-Agent": self.tns_marker}
        query_json = OrderedDict(objinfo_dict)
        search_data = {"api_key": self.api_key, "data": json.dumps(query_json)}
        response = requests.post(search_url, headers=headers, data=search_data)
        json_data = json.loads(response.text)
        return json_data

    def search_by_date(self, 
                       since_date : str = None):
        """
        Query recent transients from TNS within the specified date range.

        Parameters:
            start_date (str): Start date in ISO format (YYYY-MM-DD).
            end_date (str): End date in ISO format (default is the current date).

        Returns:
            pd.DataFrame: DataFrame containing objname, ra, and dec.
        """
        if since_date is None:
            since_date = (Time.now() - 5 * u.day).datetime.strftime("%Y-%m-%d")

        objinfo_dict = [
            ("public_timestamp", f"{since_date}")
        ]
        search_url = f"{self.base_url}/search"
        headers = {"User-Agent": self.tns_marker}
        query_json = OrderedDict(objinfo_dict)
        search_data = {"api_key": self.api_key, "data": json.dumps(query_json)}
        response = requests.post(search_url, headers=headers, data=search_data)
        json_data = json.loads(response.text)
        data = response.get("data", {}).get("reply", [])
        
        # Extract relevant information
        transients = [
            {
                "objname": obj["objname"],
                "ra": float(obj["ra"]),
                "dec": float(obj["dec"])
            }
            for obj in data if "ra" in obj and "dec" in obj
        ]
        return pd.DataFrame(transients)
    

    def save_to_csv(self, df, filename):
        """
        Save a DataFrame to a CSV file.

        Parameters:
            df (pd.DataFrame): DataFrame to save.
            filename (str): Name of the output CSV file.

        Returns:
            str: Path to the saved CSV file.
        """
        filepath = os.path.join(self.download_dir, filename)
        df.to_csv(filepath, index=False)
        print(f"Saved data to {filepath}")
        return filepath


# Example usage
if __name__ == "__main__":
    # Initialize with your API key and bot credentials
    api_key = "9da2ee1367ba2c229d1c5d6c4aebeca425a5d67c"
    bot_id = "117913"
    bot_name = "HH_bot"
    download_dir = "./downloaded_files"
    tns_querier = TNSQuerier(api_key, bot_id, bot_name, download_dir)
# %%
