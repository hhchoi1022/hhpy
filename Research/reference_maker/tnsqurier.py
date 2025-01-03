#%%
import os
import json
import time
import requests
from collections import OrderedDict
from datetime import datetime
from astropy.table import Table
from astropy.io import ascii
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import vstack
from astropy.time import Time


class TNSQuerier:
    """
    A class for querying the TNS server and handling results.
    """

    def __init__(self, api_key, bot_id, bot_name):
        """
        Initialize the TNSQuerier class.

        Parameters:
            api_key (str): The API key for TNS.
            bot_id (str): The bot ID for TNS.
            bot_name (str): The bot name.
            user_id (str): User ID (if querying as a user).
            user_name (str): User name (if querying as a user).
        """
        self.api_key = api_key
        self.bot_id = bot_id
        self.bot_name = bot_name

        self.base_url = "https://www.wis-tns.org"
        self.search_url = f"{self.base_url}/search"
        self.headers = {'User-Agent': self._set_tns_marker()}

    def _set_tns_marker(self):
        """
        Set the TNS marker for API requests.

        Parameters:
            bot (bool): If True, use bot credentials. If False, use user credentials.

        Returns:
            str: TNS marker string.
        """
        return f'tns_marker{{"tns_id": "{self.bot_id}", "type": "bot", "name": "{self.bot_name}"}}'
    
    
    # function for searching TNS with specified url parameters
    def search_tns(self, url_parameters, 
                   alert_save_dir : str = './tns_alerts'
                   #merge_to_single_file : bool = True,
                   ):
        #--------------------------------------------------------------------
        # extract keywords and values from url parameters
        keywords = list(url_parameters.keys())
        values = list(url_parameters.values())
        #--------------------------------------------------------------------
        # flag for checking if url is with correct keywords
        wrong_url = False
        # check if keywords are correct
        for i in range(len(keywords)):
            if keywords[i] not in self._all_url_params:
                print ("Unknown url keyword '"+keywords[i]+"'\n")
                wrong_url = True
        # check flag
        if wrong_url == True:
            print ("TNS search url is not in the correct format.\n")
        #--------------------------------------------------------------------
        # else, if everything is correct
        else:
            # current date and time
            current_datetime = datetime.now()
            current_date_time = current_datetime.strftime("%Y%m%d_%H%M%S")
            # current working directory
            cwd = os.getcwd()        
            if "format" in keywords:
                extension = "." + url_parameters["format"]
            else:
                extension = ".txt"
            tns_search_file = "tns_search_data_" + current_date_time + extension
            tns_search_file_path = os.path.join(cwd, alert_save_dir, tns_search_file)       
            if not os.path.exists(os.path.join(cwd, alert_save_dir)):
                os.makedirs(os.path.join(cwd, alert_save_dir))     
            #--------------------------------------------------------------------
            # build TNS search url
            url_par = ['&' + x + '=' + y for x, y in zip(keywords, values)]
            tns_search_url = self.search_url + '?' + "".join(url_par)
            #--------------------------------------------------------------------
            # page number
            page_num = 0
            # searched data
            searched_data = []
            # go trough every page
            while True:
                # url for download
                url = tns_search_url + "&page=" + str(page_num)        
                headers = self.headers
                # downloading file using request module
                response = requests.post(url, headers=headers, stream=True)
                # chek if response status code is not 200, or if returned data is empty
                if (response.status_code != 200) or (len((response.text).splitlines()) <= 1):
                    if response.status_code != 200:
                        print ("Sending download search request for page num " + str(page_num + 1) + "...")
                        self._print_response(response, page_num + 1)
                    break            
                print ("Sending download search request for page num " + str(page_num + 1) + "...")
                # print status code of the response
                self._print_response(response, page_num + 1)
                # get data
                data = (response.text).splitlines()
                # add to searched data
                if page_num == 0:
                    searched_data.append(data)
                else:
                    searched_data.append(data[1 : ])
                # check reset time
                reset = self._get_reset_time(response)
                if reset != None:
                    # Sleeping for reset + 1 sec
                    print("\nSleep for " + str(reset + 1) + " sec and then continue...\n") 
                    time.sleep(reset + 1)
                # increase page num
                page_num = page_num + 1
            #--------------------------------------------------------------------
            # if there is searched data, write to file
            if searched_data != []:            
                searched_data = [j for i in searched_data for j in i]            
                #if merge_to_single_file == 1:
                f = open(tns_search_file_path, 'w')
                for el in searched_data:
                    f.write(el + '\n')
                f.close()
                if len(searched_data) > 2:
                    print ("\nTNS searched data returned " + str(len(searched_data) - 1) + " rows. File '" + \
                        tns_search_file + "' is successfully created.\n")
                    return tns_search_file_path
                else: 
                    print ("\nTNS searched data returned 1 row. File '" + tns_search_file + "' is successfully created.\n")            
                    return tns_search_file_path
            else:
                print ("TNS searched data returned empty list. No file(s) created.\n")
    
    def update_alerthisotry(self, new_alert_path, historyfile_path : str = 'targets.ascii'):
        # Read the new alert file
        new_alert_tbl = ascii.read(new_alert_path)
        new_alert_ids = new_alert_tbl['ID']
        
        # Extract relevant fields and create a new table
        ra = new_alert_tbl['RA']
        dec = new_alert_tbl['DEC']
        coord = self._get_skycoord(ra, dec)  # Assuming this method is already defined
        new_tbl = Table()
        new_tbl['TNS_ID'] = new_alert_ids
        new_tbl['Name'] = new_alert_tbl['Name']
        new_tbl['RA'] = np.round(coord.ra.value, 5)
        new_tbl['DEC'] = np.round(coord.dec.value, 5)
        new_tbl['Discovery_UTC'] = new_alert_tbl['Discovery Date (UT)']
        
        # Check if the history file exists
        try:
            history_tbl = ascii.read(historyfile_path)
        except FileNotFoundError:
            # If no history file exists, save the new table as the history
            new_tbl.write(historyfile_path, format='ascii', overwrite=True)
            print(f"Created new history file: {historyfile_path}")
            return
        
        # Combine the new data with the history, avoiding duplicates
        history_ids = history_tbl['TNS_ID']
        unique_mask = ~np.isin(new_alert_ids, history_ids)  # Filter out IDs already in history
        unique_new_tbl = new_tbl[unique_mask]
        
        if len(unique_new_tbl) > 0:
            updated_tbl = vstack([history_tbl, unique_new_tbl])  # Merge tables
            updated_tbl.sort('Discovery_UTC', reverse=True)
            updated_tbl.write(historyfile_path, format='ascii', overwrite=True)
            print(f"Updated history file: {historyfile_path} with {len(unique_new_tbl)} new entries.")
        else:
            print("No new unique entries to add to the history file.")

    def _response_status(self, response):
        # external http errors
        ext_http_errors       = [403, 500, 503]
        err_msg               = ["Forbidden", "Internal Server Error: Something is broken", "Service Unavailable"]
        def is_string_json(string):
            try:
                json_object = json.loads(string)
            except Exception:
                return False
            return json_object

        json_string = is_string_json(response.text)
        if json_string != False:
            status = "[ " + str(json_string['id_code']) + " - '" + json_string['id_message'] + "' ]"
        else:
            status_code = response.status_code
            if status_code == 200:
                status_msg = 'OK'
            elif status_code in ext_http_errors:
                status_msg = err_msg[ext_http_errors.index(status_code)]
            else:
                status_msg = 'Undocumented error'
            status = "[ " + str(status_code) + " - '" + status_msg + "' ]"
        return status

    def _print_response(self, response, page_num):
        status = self._response_status(response)
        if response.status_code == 200:     
            stats = 'Page number ' + str(page_num) + ' | return code: ' + status + \
                    ' | Total Rate-Limit: ' + str(response.headers.get('x-rate-limit-limit')) + \
                    ' | Remaining limit: ' + str(response.headers.get('x-rate-limit-remaining')) + \
                    ' | Limit reset after ' + str(response.headers.get('x-rate-limit-reset') + ' sec')
        
        else:       
            stats = 'Page number ' + str(page_num) + ' | return code: ' + status        
        print (stats)
                    
    def _send_request(self, parameters):
        """
        Send a POST request to TNS.

        Parameters:
            parameters (dict): Request parameters.

        Returns:
            Response: The HTTP response object.
        """
        search_data = {"api_key": self.api_key, "data": json.dumps(OrderedDict(parameters))}
        return requests.post(self.search_url, headers=self.headers, data=search_data)

    def _get_reset_time(self, response):
        """
        Get the rate-limit reset time from the response headers.

        Parameters:
            response (Response): The HTTP response object.

        Returns:
            int or None: Reset time in seconds, or None if no reset is required.
        """
        # If any of the '...-remaining' values is zero, return the reset time
        for name in response.headers:
            value = response.headers.get(name)
            if name.endswith('-remaining') and value == '0':
                return int(response.headers.get(name.replace('remaining', 'reset')))
        return None     
    
    def _get_skycoord(self,
                     ra : str or float,
                     dec: str or float,
                     ra_unit : str = 'hourangle',
                     dec_unit : str = 'deg'
                     ):
        """
        Parameters
        ==========
        ra : str | float = Right Ascension, if str, it should be in hms format (e.g., "10:20:30"), if float, it should be in decimal degrees
        dec : str | float = Declination, if str, it should be in dms format (e.g., "+20:30:40"), if float, it should be in decimal degrees
        
        Return
        ======
        coord : SkyCoord = SkyCoord object
        """
        
        u_ra = u.hourangle if ra_unit == 'hourangle' else u.deg
        u_dec = u.deg if dec_unit == 'deg' else u.deg
        coord = SkyCoord(ra, dec, unit=(u_ra, u_dec))
        return coord
            
    @property
    def _all_url_params(self):
        URL_PARAMETERS   = ["discovered_period_value", "discovered_period_units", "unclassified_at", "classified_sne", "include_frb",
                            "name", "name_like", "isTNS_AT", "public", "ra", "decl", "radius", "coords_unit", "reporting_groupid[]",
                            "groupid[]", "classifier_groupid[]", "objtype[]", "at_type[]", "date_start[date]",
                            "date_end[date]",  "discovery_mag_min", "discovery_mag_max", "internal_name", "discoverer", "classifier",
                            "spectra_count", "redshift_min", "redshift_max", "hostname", "ext_catid", "ra_range_min", "ra_range_max",
                            "decl_range_min", "decl_range_max", "discovery_instrument[]", "classification_instrument[]",
                            "associated_groups[]", "official_discovery", "official_classification", "at_rep_remarks", "class_rep_remarks",
                            "frb_repeat", "frb_repeater_of_objid", "frb_measured_redshift", "frb_dm_range_min", "frb_dm_range_max",
                            "frb_rm_range_min", "frb_rm_range_max", "frb_snr_range_min", "frb_snr_range_max", "frb_flux_range_min",
                            "frb_flux_range_max", "format", "num_page"]
        return URL_PARAMETERS
    
    
#%%
# Example usage
if __name__ == "__main__":
    now = Time.now() 
    two_days_ago = now - u.days(2)
    # Example TNS parameters
    api_key = "9da2ee1367ba2c229d1c5d6c4aebeca425a5d67c"
    bot_id = "117913"
    bot_name = "HH_bot"
    tns = TNSQuerier(api_key, bot_id, bot_name)
    file_ = tns.search_tns(url_parameters = {"date_start[date]": two_days_ago, "date_end[date]": now.isot, "num_page": '3', "format": 'csv'}, alert_save_dir = '/home/7dt/7dt_too/backend/data/target_hisotry')
    file2_ = tns.update_alerthisotry(file_, historyfile_path = '/home/7dt/7dt_too/backend/data/target_hisotry/targets.ascii')
# %%
