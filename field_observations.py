import numpy as np
import scipy.stats as stats
from scipy import spatial
import matplotlib.pyplot as plt
import matplotlib.colors as clb
from mpl_toolkits.basemap import Basemap, cm, interp
import scipy.interpolate as intrplt

import jsbach_analysis as jsbach
import clm_analysis as clm
import ctem_analysis as ctem
import blaze_analysis as blaze
import orchidee_analysis as orchidee
import inferno_analysis as inferno
import spitfire_analysis as spitfire
import globfirm_analysis as globfirm

import gfed_analysis as gfed

import spatial_comparison as spt

observ_latlon = [[-37.09,145.08], [-37.09,145.08]]
observ_FC = [2,4]

def find_indices(model):
    lons_list, lats_list = spt.get_lons_lats(model)
    lats, lons = np.meshgrid(lats_list,lons_list)
    latlon_data = zip(lats.ravel(), lons.ravel())
    model_grid = spatial.KDTree(latlon_data)
    query = model_grid.query(observ_latlon)
    indices = query[1]
    for index in indices:
        latlon = model_grid.data[index]
        print lons_list.tolist().index(latlon[1])
        print lats_list.index(latlon[0])
    
    
    model_shape = model_grid.data.shape
    return observ_indices, model_shape

def get_observ_grid(model):
    indices, shape = find_indices(model)
    grid = np.empty(shape)
    for i in range(len(indices)):
        index = indices[i]
        value = observ_FC[i]
        if grid[index] == 0:
            grid[index] = [value]
        else:
            grid[index].append(value)
    return grid
    
#print get_observ_grid('ctem')















