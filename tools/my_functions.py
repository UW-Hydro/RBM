#!/usr/local/anaconda/bin/python

# -------------------------------------------------------------------- #
def read_config(config_file, default_config=None):
    """
    Return a dictionary with subdictionaries of all configFile options/values
    """

    from netCDF4 import Dataset
    try:
        from cyordereddict import OrderedDict
    except:
        from collections import OrderedDict
    try:
        from configparser import SafeConfigParser
    except:
        from ConfigParser import SafeConfigParser
    import configobj

    config = SafeConfigParser()
    config.optionxform = str
    config.read(config_file)
    sections = config.sections()
    dict1 = OrderedDict()
    for section in sections:
        options = config.options(section)
        dict2 = OrderedDict()
        for option in options:
            dict2[option] = config_type(config.get(section, option))
        dict1[section] = dict2

    if default_config is not None:
        for name, section in dict1.items():
            if name in default_config.keys():
                for option, key in default_config[name].items():
                    if option not in section.keys():
                        dict1[name][option] = key

    return dict1
# -------------------------------------------------------------------- #

# -------------------------------------------------------------------- #
def config_type(value):
    """
    Parse the type of the configuration file option.
    First see the value is a bool, then try float, finally return a string.
    """
    val_list = [x.strip() for x in value.split(',')]
    if len(val_list) == 1:
        value = val_list[0]
        if value in ['true', 'True', 'TRUE', 'T']:
            return True
        elif value in ['false', 'False', 'FALSE', 'F']:
            return False
        elif value in ['none', 'None', 'NONE', '']:
            return None
        else:
            try:
                return int(value)
            except:
                pass
            try:
                return float(value)
            except:
                return value
    else:
        try:
            return list(map(int, val_list))
        except:
            pass
        try:
            return list(map(float, val_list))
        except:
            return val_list
# -------------------------------------------------------------------- #

#============================================================#
#============================================================#

def modify_hydraulics_at_reservoir(lat, lon, depth, width, year_operated, da_depth, da_width, da_velocity, da_flow, min_velocity):
    ''' This function modifies flow hydraulics at a reservoir grid cell
    Input:
        lat, lon: lat and lon of the grid cell of the reservoir [float]
        depth, width: depth and width for this reservoir [ft]
        year_operated: the calendar year the beginning of which reservoir operation started [int]
        da_depth, da_width, da_flow: xray.DataArray of flow depth, width, velocity and flow discharge [ft and s] (must have the same time indices)

    Return:
        da_depth, da_width, da_velocity: modified xray.DataArray
    '''

    import pandas as pd
    import datetime as dt
   
    #=== Extract time index of data ===#
    dates = da_depth['time'].to_series().index
    start_date = dates[0]  # start datetime of data
    end_date = dates[-1]  # end datetime of data

    #=== Calculate new depth, width and velocity ===#
    if year_operated<=start_date.year:  # if reservoir operation started before simulation period,
                                        # modify hydraulics for the whole period
        s_depth = pd.Series(depth, index=dates)
        s_width = pd.Series(width, index=dates)
        s_velocity = da_flow.loc[:,lat,lon] / s_depth / s_width
        s_velocity[s_velocity<min_velocity] = min_velocity
    else:  # if reservoir operation started after simulation period,
           # modify hydraulics only from the year starting operation
        #=== Modify depth ===#
        s_depth = da_depth.loc[:,lat,lon].to_series()  # original depth
        dates_to_modify = s_depth.truncate(before=dt.datetime(year_operated,1,1)).index
                            # dates after reservoir operation began
        s_depth.loc[dates_to_modify] = depth  # modify depth
        #=== Modify width ===#
        s_width = da_width.loc[:,lat,lon].to_series()  # originalwidth 
        s_width.loc[dates_to_modify] = width  # modify width 
        #=== Modify velocity ===#
        s_velocity = da_flow.loc[:,lat,lon] / s_depth / s_width
        s_velocity[s_velocity<min_velocity] = min_velocity

    #=== Modify DataArray ===#
    da_depth.loc[:,lat,lon] = s_depth
    da_width.loc[:,lat,lon] = s_width
    da_velocity.loc[:,lat,lon] = s_velocity
    
    return da_depth, da_width, da_velocity 






