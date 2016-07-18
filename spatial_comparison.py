import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
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


def interp_regions():
    regions_GFED = gfed.grid_GFED["basis_regions"]
    
    lats_gfed = gfed.data_GFED["lat"][:,0]
    #lats_gfed = np.sort(lats_gfed)
    lons_gfed = gfed.data_GFED["lon"][0,:]
    #lons_gfed = np.sort(lons_gfed)
    #lons_gfed, lats_gfed = np.meshgrid(lons_gfed, lats_gfed)
    
    lats = jsbach.grid_JSBACH["latitude"]
    lons = jsbach.grid_JSBACH["longitude"]
    #lons, lats = np.meshgrid(lons, lats)
    
    #regions = interp(regions_GFED, 
    #                    lons_gfed,lats_gfed,
    #                    lons, lats,checkbounds=False, masked=False, order=0)
    
    regions = intrplt.griddata([lons_gfed,lats_gfed], regions_GFED,
                            [lons,lats], method='nearest')

    return regions


def plot_regions():
    map_data = gfed.grid_GFED["basis_regions"]
    map_data = interp_regions()
    
    lats = gfed.data_GFED["lat"][:,0]
    lons = gfed.data_GFED["lon"][0,:]
    lons, lats = np.meshgrid(lons, lats)

    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-60, 
        urcrnrlon=180,urcrnrlat=80,projection='mill')
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs = m.contourf(lons,lats, map_data, 100, cmap=plt.cm.YlOrRd, latlon=True)
    cb = m.colorbar(cs, "bottom", size="5%", pad="2%")
    cb.set_label("Burnt Area (billions of $m^2$)")
    plt.title("Mean Burnt Area, LPJ BLAZE")
    plt.show()

#plot_regions()



def load_var_grid(year, year_period, model, var='FC'):
    """
    Year is in absolute terms, e.g. 1997.
    
    Takes argument model, which is a string that can be one of the
    following: 'gfed', 'jsbach', 'clm', 'ctem', 'blaze', 'orchidee',
    'inferno'
    
    Also takes argument var, which is a string that can take either
    'FC', 'BA', or 'emis' as an input, which will return the
    spatial histogram for fuel consumption, burnt area, or emissions
    respectively.
    """
    year_adj = year-1700
    year_ctem = year-1861
    year_gfed = year-1997
    if var == 'FC':
        if model == 'gfed':
            grid = gfed.get_grid_fuel_consumption(year_gfed,year_period*12,
                                                     gfed.data_GFED)
        elif model == 'jsbach':
            grid = jsbach.get_grid_fuel_consumption(year_adj,year_period*12,
                                     jsbach.emis_JSBACH,jsbach.BA_JSBACH)
        elif model == 'clm':
            grid = clm.get_grid_fuel_consumption(year_adj,year_period*12,
                               clm.emis_CLM,clm.BA_CLM,clm.time_data)
        elif model == 'ctem':
            grid = ctem.get_grid_fuel_consumption(year_ctem,year_period*12,
                         ctem.emis_CTEM,ctem.BA_CTEM,ctem.landCover_CTEM)
        elif model == 'blaze':
            grid = blaze.get_grid_fuel_consumption(year_adj,year_period*12,
                               blaze.emis_BLAZE,blaze.BA_BLAZE,blaze.time_data)
        elif model == 'orchidee':
            grid = orchidee.get_grid_fuel_consumption(year_adj,year_period*12,
                     orchidee.emis_ORCHIDEE,orchidee.BA_ORCHIDEE,
                     orchidee.landCover_ORCHIDEE,orchidee.time_data)
        elif model == 'inferno':
            grid = inferno.get_grid_fuel_consumption(year_adj,year_period*12,
                     inferno.emis_Inferno,inferno.BA_Inferno,
                     inferno.landmask_Inferno,inferno.landCover_Inferno)
    elif var == 'emis':
        if model == 'gfed':
            grid = gfed.get_grid_emissions(year_gfed,year_period*12,
                                       gfed.data_GFED, gfed.grid_GFED)
        elif model == 'jsbach':
            grid = jsbach.get_grid_emissions(year_adj,year_period*12,
                                         jsbach.emis_JSBACH, jsbach.grid_JSBACH)
        elif model == 'clm':
            grid = clm.get_grid_emissions(year_adj,year_period*12,
                                clm.emis_CLM,clm.grid_CLM,clm.time_data)
        elif model == 'ctem':
            grid = ctem.get_grid_emissions(year_ctem,year_period*12,
                         ctem.emis_CTEM,ctem.grid_CTEM,ctem.landCover_CTEM)    
        elif model == 'blaze':
            grid = blaze.get_grid_emissions(year_adj,year_period*12,
                               blaze.emis_BLAZE,blaze.grid_BLAZE,blaze.time_data)
        elif model == 'orchidee':
            grid = orchidee.get_grid_emissions(year_adj,year_period*12,
                     orchidee.emis_ORCHIDEE,orchidee.grid_ORCHIDEE,
                     orchidee.landCover_ORCHIDEE,orchidee.time_data)
        elif model == 'inferno':
            grid = inferno.get_grid_emissions(year_adj,year_period*12,
                     inferno.emis_Inferno,inferno.grid_Inferno,
                     inferno.landmask_Inferno,inferno.landCover_Inferno)
    elif var == 'BA':
        if model == 'gfed':
            grid = gfed.get_grid_burnt_area(year_gfed,year_period*12,
                                          gfed.data_GFED, gfed.grid_GFED)
        elif model == 'jsbach':
            grid = jsbach.get_grid_burnt_area(year_adj,year_period*12,
                                         jsbach.BA_JSBACH, jsbach.grid_JSBACH)
        elif model == 'clm':
            grid = clm.get_grid_burnt_area(year_adj,year_period*12,
                                clm.BA_CLM,clm.grid_CLM,clm.time_data)
        elif model == 'ctem':
            grid = ctem.get_grid_burnt_area(year_ctem,year_period*12,
                         ctem.BA_CTEM,ctem.grid_CTEM,ctem.landCover_CTEM)
        elif model == 'blaze':
            grid = blaze.get_grid_burnt_area(year_adj,year_period*12,
                               blaze.BA_BLAZE,blaze.grid_BLAZE)
        elif model == 'orchidee':
            grid = orchidee.get_grid_burnt_area(year_adj,year_period*12,
                     orchidee.BA_ORCHIDEE,orchidee.grid_ORCHIDEE,
                     orchidee.landCover_ORCHIDEE)
        elif model == 'inferno':
            grid = inferno.get_grid_burnt_area(year_adj,year_period*12,
                     inferno.BA_Inferno,inferno.grid_Inferno,
                     inferno.landmask_Inferno,inferno.landCover_Inferno)
                     
    return grid


def plot_spatial_histogram(year, year_period, model, var='FC'):
    """
    Year is in absolute terms, e.g. 1997.
    
    Takes argument model, which is a string that can be one of the
    following: 'gfed', 'jsbach', 'clm', 'ctem', 'blaze', 'orchidee',
    'inferno'
    
    Also takes argument var, which is a string that can take either
    'FC', 'BA', or 'emis' as an input, which will return the
    spatial histogram for fuel consumption, burnt area, or emissions
    respectively.
    """
    grid = load_var_grid(year, year_period, model, var)
    
    if var == 'FC':
        title = 'Fuel Consumption'
        units = '($kg\, C\, m^{-2}\, burned$)'
    elif var == 'emis':
        title = 'Carbon Emissions'
        units = '($kg\, C$)'
    elif var == 'BA':
        title = 'Burnt Area'
        units = '($m^2$)'

    #flat_grid = np.ndarray.flatten(grid)
    flat_grid = grid[grid > 0]
    flat_grid = np.divide(flat_grid, year_period)
    plt.hist(flat_grid, 50, normed=1, facecolor='blue', alpha=0.75)
    
    plt.xlabel(title + ' ' + units)
    plt.ylabel('Frequency in \%')
    plt.title('Spatial Histogram of ' + title + ' for ' + model)

    plt.show()


def plot_multimodel_box(year, year_period, var='FC'):
    """
    Year is given in absolute terms, e.g. 1997.
    
    Also takes argument var, which is a string that can take either
    'FC', 'BA', or 'emis' as an input, which will return the
    boxplots for fuel consumption, burnt area, or emissions
    respectively.
    """
    model_list = ['gfed', 'jsbach', 'clm', 'ctem', 
                'blaze', 'orchidee', 'inferno']
    data=[]              
    for model in model_list:
        grid = load_var_grid(year, year_period, model, var)
        flat_grid = grid[grid > 0]
        flat_grid = np.divide(flat_grid, year_period)
        data.append(flat_grid)
      
    if var == 'FC':
        title = 'Fuel Consumption'
        units = '($kg\, C\, m^{-2}\, burned$)'
    elif var == 'emis':
        title = 'Carbon Emissions'
        units = '($kg\, C$)'
    elif var == 'BA':
        title = 'Burnt Area'
        units = '($m^2$)'
   
    plt.boxplot(data, showmeans=True, whis=[5,95], sym='', 
                labels = model_list)
    plt.ylabel(title +' '+ units)
    plt.show()



#plot_multimodel_box(1997,15,'FC')
#plot_spatial_histogram(1997, 15, 'orchidee', 'FC')
