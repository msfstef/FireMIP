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


#
# Regional Analysis Toolkit
#

def interp_regions(lons, lats):
    regions_GFED = gfed.grid_GFED["basis_regions"]
    regions_GFED = np.array(regions_GFED)
        
    lats_gfed = gfed.data_GFED["lat"]
    lats_gfed = np.array(lats_gfed)
    lons_gfed = gfed.data_GFED["lon"]
    lons_gfed = np.array(lons_gfed)
    
    regions = intrplt.griddata((lons_gfed.ravel(),lats_gfed.ravel()), 
                  regions_GFED.ravel(), (lons,lats), method='nearest')
    return regions


def plot_regions(model='gfed', plot=True):
    """
    Takes argument model, which is a string that can be one of the
    following: 'gfed', 'jsbach', 'clm', 'ctem', 'blaze', 'orchidee',
    'inferno'
    """
    no_interp = False
    if model=='gfed':
        lats = gfed.data_GFED["lat"][:,0]
        lons = gfed.data_GFED["lon"][0,:]
        no_interp = True
    elif model=='jsbach':
        lats = jsbach.grid_JSBACH["latitude"]
        lats = np.array(lats)
        lats = lats[::-1]
        lons = jsbach.grid_JSBACH["longitude"]
        lons = np.array(lons) - 180.
    elif model=='clm':
        lats = clm.grid_CLM["lat"]
        lons = clm.grid_CLM["lon"]
        lons = np.array(lons) - 180.
        
    elif model=='ctem':
        lats = ctem.grid_CTEM["lat"]
        lons = ctem.grid_CTEM["lon"]
        lons = np.array(lons) - 180.
    elif model=='blaze':
        lats = blaze.grid_BLAZE["lat"]
        lons = blaze.grid_BLAZE["lon"]
    elif model=='orchidee':
        lats = orchidee.grid_ORCHIDEE["latitude"]
        lats = np.array(lats)
        lats = lats[::-1]
        lons = orchidee.grid_ORCHIDEE["longitude"]
    elif model=='inferno':
        lats = inferno.grid_Inferno["latitude"]
        lons = inferno.grid_Inferno["longitude"]  
        lons = np.array(lons) - 180.
      
    lons, lats = np.meshgrid(lons, lats)
    
    if no_interp:
        map_data = gfed.grid_GFED["basis_regions"]
        map_data = map_data[::-1,:]
    else:
        map_data = interp_regions(lons, lats)
    
    if plot:
        fig=plt.figure()
        m = Basemap(llcrnrlon=-180,llcrnrlat=-90, 
            urcrnrlon=180,urcrnrlat=90)
        m.drawcoastlines()
        m.drawparallels(np.arange(-90.,91.,30.))
        m.drawmeridians(np.arange(-180.,181.,60.))
        m.drawmapboundary(fill_color='white')
        #cs = m.contourf(lons,lats, map_data, 100, cmap=plt.cm.rainbow, latlon=True)
        m.imshow(map_data, cmap=plt.cm.rainbow,  interpolation='none')
        #cb = m.colorbar(cs, "bottom", size="5%", pad="2%")
        #cb.set_label("Region number")
        plt.title("GFED Regions Interpolated for " + model)
        plt.show()


#plot_regions('jsbach')


#
# Spatial Global Analysis Toolkit
#

def load_var_grid(year, year_period, model, var='FC', for_map=False):
    """
    Year is in absolute terms, e.g. 1997.
    
    Takes argument model, which is a string that can be one of the
    following: 'gfed', 'jsbach', 'clm', 'ctem', 'blaze', 'orchidee',
    'inferno'
    
    Also takes argument var, which is a string that can take either
    'FC', 'BA', or 'emis' as an input, which will return the
    spatial histogram for fuel consumption, burnt area, or emissions
    respectively. Default is fuel consumption.
    
    for_map argument determines whether variables should be returned in
    units of var/m^2 or in absolute terms. Default is False.
    """
    year_adj = year-1700
    year_ctem = year-1861
    year_gfed = year-1997
    if var == 'FC':
        if model == 'gfed':
            grid = gfed.get_grid_fuel_consumption(year_gfed,year_period*12,
                                                     gfed.data_GFED)
            # Convert to standard format.
            grid = grid[::-1,:]
        elif model == 'jsbach':
            grid = jsbach.get_grid_fuel_consumption(year_adj,year_period*12,
                                     jsbach.emis_JSBACH,jsbach.BA_JSBACH)
            # Convert to standard format.
            grid = grid[::-1,:]
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
        elif model == 'clm':
            grid = clm.get_grid_fuel_consumption(year_adj,year_period*12,
                               clm.emis_CLM,clm.BA_CLM,clm.time_data)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
        elif model == 'ctem':
            grid = ctem.get_grid_fuel_consumption(year_ctem,year_period*12,
                         ctem.emis_CTEM,ctem.BA_CTEM,ctem.landCover_CTEM)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
        elif model == 'blaze':
            grid = blaze.get_grid_fuel_consumption(year_adj,year_period*12,
                               blaze.emis_BLAZE,blaze.BA_BLAZE,blaze.time_data)
        elif model == 'orchidee':
            grid = orchidee.get_grid_fuel_consumption(year_adj,year_period*12,
                     orchidee.emis_ORCHIDEE,orchidee.BA_ORCHIDEE,
                     orchidee.landCover_ORCHIDEE,orchidee.time_data)
            # Convert to standard format.
            grid = grid[::-1,:]
        elif model == 'inferno':
            grid = inferno.get_grid_fuel_consumption(year_adj,year_period*12,
                     inferno.emis_Inferno,inferno.BA_Inferno,
                     inferno.landmask_Inferno,inferno.landCover_Inferno)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            
    elif var == 'emis':
        if model == 'gfed':
            grid = gfed.get_grid_emissions(year_gfed,year_period*12,
                                       gfed.data_GFED, gfed.grid_GFED)
            # Convert to standard format.
            grid = grid[::-1,:]
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, gfed.grid_GFED["grid_cell_area"])
        elif model == 'jsbach':
            grid = jsbach.get_grid_emissions(year_adj,year_period*12,
                                         jsbach.emis_JSBACH, jsbach.grid_JSBACH)
            # Convert to standard format.
            grid = grid[::-1,:]
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, jsbach.grid_JSBACH["area"])
        elif model == 'clm':
            grid = clm.get_grid_emissions(year_adj,year_period*12,
                                clm.emis_CLM,clm.grid_CLM,clm.time_data)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, clm.grid_CLM["cell_area"])
        elif model == 'ctem':
            grid = ctem.get_grid_emissions(year_ctem,year_period*12,
                         ctem.emis_CTEM,ctem.grid_CTEM,ctem.landCover_CTEM)    
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, ctem.grid_CTEM["cell_area"])
        elif model == 'blaze':
            grid = blaze.get_grid_emissions(year_adj,year_period*12,
                               blaze.emis_BLAZE,blaze.grid_BLAZE,blaze.time_data)
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, blaze.grid_BLAZE["cell_area"])
        elif model == 'orchidee':
            grid = orchidee.get_grid_emissions(year_adj,year_period*12,
                     orchidee.emis_ORCHIDEE,orchidee.grid_ORCHIDEE,
                     orchidee.landCover_ORCHIDEE,orchidee.time_data)
            # Convert to standard format.
            grid = grid[::-1,:]
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, orchidee.grid_ORCHIDEE["cell_area"])
        elif model == 'inferno':
            grid = inferno.get_grid_emissions(year_adj,year_period*12,
                     inferno.emis_Inferno,inferno.grid_Inferno,
                     inferno.landmask_Inferno,inferno.landCover_Inferno)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, inferno.grid_Inferno["cell_area"])
                
    elif var == 'BA':
        if model == 'gfed':
            grid = gfed.get_grid_burnt_area(year_gfed,year_period*12,
                                          gfed.data_GFED, gfed.grid_GFED)
            # Convert to standard format.
            grid = grid[::-1,:]
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, gfed.grid_GFED["grid_cell_area"])
        elif model == 'jsbach':
            grid = jsbach.get_grid_burnt_area(year_adj,year_period*12,
                                         jsbach.BA_JSBACH, jsbach.grid_JSBACH)
            # Convert to standard format.
            grid = grid[::-1,:]
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, jsbach.grid_JSBACH["area"])
        elif model == 'clm':
            grid = clm.get_grid_burnt_area(year_adj,year_period*12,
                                clm.BA_CLM,clm.grid_CLM,clm.time_data)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, clm.grid_CLM["cell_area"])
        elif model == 'ctem':
            grid = ctem.get_grid_burnt_area(year_ctem,year_period*12,
                         ctem.BA_CTEM,ctem.grid_CTEM,ctem.landCover_CTEM)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, ctem.grid_CTEM["cell_area"])
        elif model == 'blaze':
            grid = blaze.get_grid_burnt_area(year_adj,year_period*12,
                               blaze.BA_BLAZE,blaze.grid_BLAZE)
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, blaze.grid_BLAZE["cell_area"])
        elif model == 'orchidee':
            grid = orchidee.get_grid_burnt_area(year_adj,year_period*12,
                     orchidee.BA_ORCHIDEE,orchidee.grid_ORCHIDEE,
                     orchidee.landCover_ORCHIDEE)
            # Convert to standard format.
            grid = grid[::-1,:]
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, orchidee.grid_ORCHIDEE["cell_area"])
        elif model == 'inferno':
            grid = inferno.get_grid_burnt_area(year_adj,year_period*12,
                     inferno.BA_Inferno,inferno.grid_Inferno,
                     inferno.landmask_Inferno,inferno.landCover_Inferno)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if for_map:
                grid = np.divide(grid, inferno.grid_Inferno["cell_area"])
                     
    return grid


def plot_map(year, year_period, model, var='FC'):
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
    grid = load_var_grid(year, year_period, model, var, True)
    
    # Convert zeroes to NaNs and take yearly mean.
    grid[grid==0]=np.nan
    grid = np.divide(grid, year_period)
    
    if var == 'FC':
        title = 'Fuel Consumption'
        units = '($kg\, C\, m^{-2}\, burned$)'
    elif var == 'emis':
        title = 'Carbon Emissions'
        units = '($kg\, C\, m^{-2}$)'
    elif var == 'BA':
        title = 'Burnt Area'
        units = '(Fraction Burned)'
        
    
    
    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90, 
    urcrnrlon=180,urcrnrlat=90)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs=m.imshow(grid, interpolation='none')
    cb=m.colorbar(cs, "bottom", size="5%", pad="2%")
    cb.set_label(title + ' ' + units)
    plt.title('Mean '+title+' for '+str(year)+
            '-'+str(year+year_period)+', '+model)
    plt.show()

plot_map(1997,1,'blaze','BA')

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
#plot_spatial_histogram(1997, 15, 'clm', 'FC')
