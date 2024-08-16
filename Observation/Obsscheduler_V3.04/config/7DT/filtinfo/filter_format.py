#%%
import json
#%%

def get_filtformat(filtformat_file = './filtformat.json'):
    with open(filtformat_file, 'r') as f:
        filtformat_dict = json.load(f)
    return filtformat_dict
def get_seformat(filter_, seformat_file = './seformat.json'):
    with open(seformat_file, 'r') as f:
        seformat_dict = json.load(f)
    seformat_dict['Items']['$values'][0]['Filter']['_name'] = filter_
    return seformat_dict
def get_filtinfo(filtinfo_file = './filtinfo.config'):
    with open(filtinfo_file, 'r') as f:
        filtinfo_dict = json.load(f)
    return filtinfo_dict
#%%
filtformat_file = './filtformat.json'
seformat_file = './seformat.json'
filtinfo_file = './filtinfo.config'
result_dict = dict()
filtinfo_dict = get_filtinfo(filtinfo_file)
for name_telescope in filtinfo_dict.keys():
    filtformat_dict = get_filtformat(filtformat_file)
    all_filterset = filtinfo_dict[name_telescope]
    for filter_ in all_filterset:
        se_dict = get_seformat(filter_ = filter_, seformat_file = seformat_file)
        filtformat_dict['Items']['$values'][1]['Items']['$values'].append(se_dict)
    result_dict[name_telescope] = filtformat_dict
for name_telescope in filtinfo_dict.keys():
    filtformat_dict = result_dict[name_telescope]
    path = f'./{name_telescope}'
    filepath = f'{path}/filter.json'
    os.makedirs(path, exist_ok = True)
    with open(filepath, 'w') as f:
        json.dump(filtformat_dict, f, indent = 2)
for name_telescope in filtinfo_dict.keys():
    filtformat_dict = result_dict[name_telescope]
    str_allfilter = str()
    for i in range(len(filtformat_dict['Items']['$values'][1]['Items']['$values'])):
        str_allfilter += str(filtformat_dict['Items']['$values'][1]['Items']['$values'][i]['Items']['$values'][0]['Filter']['_name']) + ' '
    print(f'{name_telescope} = {str_allfilter}')

# %%
