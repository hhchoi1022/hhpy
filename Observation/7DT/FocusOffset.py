#%%
from tcspy.utils.target.db_target.GS_connector import GoogleSheet
gs = GoogleSheet()
name_sheet = 'Focustest'
tbl = gs.get_sheet_data(sheet_name = name_sheet, format_ = 'Table')
# %%
import datetime
all_obs_datetime = [datetime.datetime.strptime(date, '%y/%m/%d %H:%M') for date in tbl['obs_time']]
all_obs_same_datetime = [datetime.datetime(year = 2023, month = 10, day = 10, hour = date.hour, minute = date.minute) for 
date in all_obs_datetime]
tbl['obs_date'] = all_obs_same_datetime
# %%
tbl_filters = tbl.group_by('filter')
# %%
import matplotlib.pyplot as plt
for i in range(len(tbl_filters.groups)):
    filter_ = tbl_filters.groups.keys[i][0]
    tbl = tbl_filters.groups[i]
    plt.figure(dpi = 300)
    plt.title(filter_)
    plt.scatter(tbl['obs_time'], np.array(tbl['focuspos']).astype('float'))
    plt.xticks(rotation = 45)
