"""
This module is used for spatial analysis of the
models. For a temporal analysis, use the
model_comparison module.
"""

import numpy as np
import scipy.stats as stats
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


#
# Regional Analysis Toolkit
#

def create_box_regions(lons,lats):
    """
    Roughly based on GFED regions and
    GPP (gross primary production) data.
    See http://tinyurl.com/h6zvduw.
    """
    regions = np.empty(lons.shape)
    for index,lon in np.ndenumerate(lons):
        lat = lats[index]
        if -170.<=lon<-30.:
            if 45.<=lat<90.:
                regions[index]=1
            elif 20.<=lat<45.:
                regions[index]=2
            elif -25.<=lat<20.:
                regions[index]=3
            else:
                regions[index]=4
        elif -30.<=lon<60:
            if 45.<=lat<90.:
                regions[index]=5
            elif 15.<=lat<45.:
                regions[index]=6
            elif -15.<=lat<15.:
                regions[index]=7
            else:
                regions[index]=8
        else:
            if 45.<=lat<90.:
                regions[index]=9
            elif 20.<=lat<50.:
                regions[index]=10
            elif -10.<=lat<20.:
                regions[index]=11
            else:
                regions[index]=12
    return regions


def interp_regions(lons, lats, var='none'):
    regions_GFED = gfed.grid_GFED["basis_regions"]
    regions_GFED = np.array(regions_GFED)

    lats_gfed = np.arange(-89.875, 90.,0.25)
    lats_gfed = lats_gfed[::-1]
    lons_gfed = np.arange(-179.875, 180.,0.25)
    lons_gfed, lats_gfed = np.meshgrid(lons_gfed, lats_gfed)
    
    regions = intrplt.griddata((lons_gfed.ravel(),lats_gfed.ravel()), 
                  regions_GFED.ravel(), (lons,lats), method='nearest')
    return regions


def get_lons_lats(model):
    """
    Takes model string as argument and returns the
    desired model's longitude and latitude array
    for use in interpolation.
    """
    if model=='gfed':
        lats = np.arange(-89.875, 90.,0.25)
        lons = np.arange(-179.875, 180.,0.25)
    elif model=='jsbach':
        lats = jsbach.grid_JSBACH["latitude"]
        lons = jsbach.grid_JSBACH["longitude"]
        lats = np.array(lats)
        lats = lats[::-1]
        lons = np.array(lons) - np.max(lons)/2.
    elif model=='clm':
        lats = clm.grid_CLM["lat"]
        lons = clm.grid_CLM["lon"]
        lons = np.array(lons) - np.max(lons)/2.
        
    elif model=='ctem':
        lats = ctem.grid_CTEM["lat"]
        lons = ctem.grid_CTEM["lon"]
        lons = np.array(lons) - np.max(lons)/2.
    elif model=='blaze':
        lats = blaze.grid_BLAZE["lat"]
        lons = blaze.grid_BLAZE["lon"]
    elif model=='orchidee':
        lats = orchidee.grid_ORCHIDEE["latitude"]
        lats = np.array(lats)
        lats = lats[::-1]
        lons = orchidee.grid_ORCHIDEE["longitude"]
    elif model=='inferno':
        lats = inferno.grid_INFERNO["latitude"]
        lons = inferno.grid_INFERNO["longitude"]  
        lons = np.array(lons) - np.max(lons)/2.
    return lons, lats


def generate_regions(model='gfed', reg_type='boxes', plot=False):
    """
    Takes argument model, which is a string that can be one of the
    following: 'gfed', 'jsbach', 'clm', 'ctem', 'blaze', 'orchidee',
    'inferno'
        
    Argument reg_type can be set to either 'gfed' or 'boxes' to
    generate the interpolated GFED regions or boxed latlon regions
    respectively. The boxed latlon regions are the preferred method
    and are the default.
    
    Argument plot can be set to True to show regions on map.
    """
    no_interp = False
    if reg_type=='gfed' and model=='gfed':
        no_interp = True
    else:
        lons, lats = get_lons_lats(model)
      
    lons, lats = np.meshgrid(lons, lats)
    
    if no_interp:
        region_data = gfed.grid_GFED["basis_regions"]
        region_data = region_data[::-1,:]
    elif reg_type=='gfed':
        region_data = interp_regions(lons, lats)
    else:
        region_data = create_box_regions(lons,lats)
    if plot:
        fig=plt.figure()
        m = Basemap(llcrnrlon=-180,llcrnrlat=-90, 
            urcrnrlon=180,urcrnrlat=90)
        m.drawcoastlines()
        m.drawparallels(np.arange(-90.,91.,30.))
        m.drawmeridians(np.arange(-180.,181.,60.))
        m.drawmapboundary(fill_color='white')
        m.imshow(region_data, cmap=plt.cm.rainbow,  interpolation='none')
        plt.title("GFED Regions Interpolated for " + model)
        plt.show()
    else:
        return region_data


def get_regional_var_grid(year, year_period, region, model, 
                            var, reg_type, per_area=True):
    full_grid = load_var_grid(year,year_period,model,var,per_area)
    region_data = generate_regions(model, reg_type)
    
    # Set desired region to 1 and every other region to 0.
    region_data[region_data != region] = 0.
    region_data[region_data == region] = 1.
    
    regional_grid = np.multiply(full_grid, region_data)
    return regional_grid
    
    
#generate_regions('inferno', plot=True)


#
# Spatial Global Analysis Toolkit
#

def load_var_grid(year, year_period, model, var='FC', per_area=True):
    """
    Year is in absolute terms, e.g. 1997.
    
    Takes argument model, which is a string that can be one of the
    following: 'gfed', 'jsbach', 'clm', 'ctem', 'blaze', 'orchidee',
    'inferno'
    
    Also takes argument var, which is a string that can take either
    'FC', 'BA', or 'emis' as an input, which will return the
    spatial histogram for fuel consumption, burnt area, or emissions
    respectively. Default is fuel consumption.
    
    per_area argument determines whether variables should be returned in
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
                     inferno.emis_INFERNO,inferno.BA_INFERNO,
                     inferno.landmask_INFERNO,inferno.landCover_INFERNO)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            
    elif var == 'emis':
        if model == 'gfed':
            grid = gfed.get_grid_emissions(year_gfed,year_period*12,
                                       gfed.data_GFED, gfed.grid_GFED)
            # Convert to standard format.
            grid = grid[::-1,:]
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, gfed.grid_GFED["grid_cell_area"])
        elif model == 'jsbach':
            grid = jsbach.get_grid_emissions(year_adj,year_period*12,
                                         jsbach.emis_JSBACH, jsbach.grid_JSBACH)
            # Convert to standard format.
            grid = grid[::-1,:]
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, jsbach.grid_JSBACH["area"])
        elif model == 'clm':
            grid = clm.get_grid_emissions(year_adj,year_period*12,
                                clm.emis_CLM,clm.grid_CLM,clm.time_data)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, clm.grid_CLM["cell_area"])
        elif model == 'ctem':
            grid = ctem.get_grid_emissions(year_ctem,year_period*12,
                         ctem.emis_CTEM,ctem.grid_CTEM,ctem.landCover_CTEM)    
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, ctem.grid_CTEM["cell_area"])
        elif model == 'blaze':
            grid = blaze.get_grid_emissions(year_adj,year_period*12,
                               blaze.emis_BLAZE,blaze.grid_BLAZE,blaze.time_data)
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, blaze.grid_BLAZE["cell_area"])
        elif model == 'orchidee':
            grid = orchidee.get_grid_emissions(year_adj,year_period*12,
                     orchidee.emis_ORCHIDEE,orchidee.grid_ORCHIDEE,
                     orchidee.landCover_ORCHIDEE,orchidee.time_data)
            # Convert to standard format.
            grid = grid[::-1,:]
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, orchidee.grid_ORCHIDEE["cell_area"])
        elif model == 'inferno':
            grid = inferno.get_grid_emissions(year_adj,year_period*12,
                     inferno.emis_INFERNO,inferno.grid_INFERNO,
                     inferno.landmask_INFERNO,inferno.landCover_INFERNO)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, inferno.grid_INFERNO["cell_area"])
                
    elif var == 'BA':
        if model == 'gfed':
            grid = gfed.get_grid_burnt_area(year_gfed,year_period*12,
                                          gfed.data_GFED, gfed.grid_GFED)
            # Convert to standard format.
            grid = grid[::-1,:]
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, gfed.grid_GFED["grid_cell_area"])
        elif model == 'jsbach':
            grid = jsbach.get_grid_burnt_area(year_adj,year_period*12,
                                         jsbach.BA_JSBACH, jsbach.grid_JSBACH)
            # Convert to standard format.
            grid = grid[::-1,:]
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, jsbach.grid_JSBACH["area"])
        elif model == 'clm':
            grid = clm.get_grid_burnt_area(year_adj,year_period*12,
                                clm.BA_CLM,clm.grid_CLM,clm.time_data)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, clm.grid_CLM["cell_area"])
        elif model == 'ctem':
            grid = ctem.get_grid_burnt_area(year_ctem,year_period*12,
                         ctem.BA_CTEM,ctem.grid_CTEM,ctem.landCover_CTEM)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, ctem.grid_CTEM["cell_area"])
        elif model == 'blaze':
            grid = blaze.get_grid_burnt_area(year_adj,year_period*12,
                               blaze.BA_BLAZE,blaze.grid_BLAZE)
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, blaze.grid_BLAZE["cell_area"])
        elif model == 'orchidee':
            grid = orchidee.get_grid_burnt_area(year_adj,year_period*12,
                     orchidee.BA_ORCHIDEE,orchidee.grid_ORCHIDEE,
                     orchidee.landCover_ORCHIDEE)
            # Convert to standard format.
            grid = grid[::-1,:]
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, orchidee.grid_ORCHIDEE["cell_area"])
        elif model == 'inferno':
            grid = inferno.get_grid_burnt_area(year_adj,year_period*12,
                     inferno.BA_INFERNO,inferno.grid_INFERNO,
                     inferno.landmask_INFERNO,inferno.landCover_INFERNO)
            # Convert to standard format.
            grid = np.roll(grid, len(grid[0,:])/2,axis=1)
            # Convert to per m^2 units if output is for map.
            if per_area:
                grid = np.divide(grid, inferno.grid_INFERNO["cell_area"])
                     
    return grid


def plot_map(year, year_period, model, var='FC', 
         region=0, reg_type='boxes', binned=False, save=False):
    """
    Year is in absolute terms, e.g. 1997.
    
    Takes argument model, which is a string that can be one of the
    following: 'gfed', 'jsbach', 'clm', 'ctem', 'blaze', 'orchidee',
    'inferno'
    
    Also takes argument var, which is a string that can take either
    'FC', 'BA', or 'emis' as an input, which will return the
    spatial histogram for fuel consumption, burnt area, or emissions
    respectively.
    
    To plot regional data using boxed regions (for no interpolation
    errors) leave the reg_type argument to the default 'boxes' and
    set the region argument from 1 to 12, which corresponds to 
    Boreal North America (BONA), Temperate North America (TENA),
    Equatorial Central and South America (EQCSA), Southernmost
    America (SOMA), Northern Europe (NOEU), Mediterranean
    and Middle East (MEME), Equatorial Africa (EQAF), Southernmost 
    Africa (SOAF), Boreal Asia (BOAS), Central Asia (CEAS),
    Equatorial Asia (EQAS), Australia/Oceania (AUST).
    
    To plot regional data using the GFED regions, set reg_type argument
    to 'gfed' and then set value of region argument from 1 to 14, 
    which correspond to BONA, TENA, CEAM, NHSA. SHSA, EURO, MIDE, 
    NHAF, SHAF, BOAS, CEAS, SEAS, EQAS, AUST respectively.
    For details regarding the GFED regions, go to 
    http://www.globalfiredata.org/data.html.
    
    If binned is set to True, the map will be plotted with binned values
    and a discrete colorbar following the GFED bin convention.
    """
    if reg_type=='boxes':
        region_names = ['Global','BONA','TENA','EQCSA','SOMA','NOEU',
                        'MEME','EQAF','SOAF','BOAS','CEAS','EQAS','AUST']
    elif reg_type=='gfed':
        region_names = ['Global','BONA','TENA','CEAM','NHSA','SHSA',
                        'EURO','MIDE','NHAF','SHAF','BOAS','CEAS',
                        'SEAS','EQAS','AUST']
                
    grid = load_var_grid(year, year_period, model, var)
    if region != 0:
        grid = get_regional_var_grid(year, year_period, region, 
                                    model, var, reg_type)
    
    # Convert zeroes to NaNs and take yearly mean.
    grid[grid==0]=np.nan
    if var != 'FC':
        grid = np.divide(grid, year_period)
    
    if var == 'FC':
        title = 'Fuel Consumption'
        units = '($kg\, C\, m^{-2}\, burned$)'
        if binned:
            ticks=[0,.05,.1,.2,.5,1.,2.,5.,10.,'>20.0']
            bounds=[0,1,2,3,4,5,6,7,8,9]
            for index,value in np.ndenumerate(grid):
                bin = 'nan'
                if 0.<value<0.05:
                    bin = 0
                elif 0.05<=value<0.1:
                    bin = 1
                elif 0.1<=value<0.2:
                    bin = 2
                elif 0.2<=value<0.5:
                    bin = 3
                elif 0.5<=value<1.:
                    bin = 4
                elif 1.<=value<2.:
                    bin = 5
                elif 2.<=value<5.:
                    bin = 6
                elif 5.<=value<10.:
                    bin = 7
                elif 10.<=value:
                    bin = 8
                if bin != 'nan':    
                    grid[index] = bin
    elif var == 'emis':
        title = 'Carbon Emissions'
        units = '($kg\, C\, m^{-2}\, year^{-1}$)'
        if binned:
            ticks=[0.,.005,.01,.02,.05,.1,.2,.5,1.,'>5.0']
            bounds=[0,1,2,3,4,5,6,7,8,9]
            for index,value in np.ndenumerate(grid):
                bin = 'nan'
                if 0.<value<0.005:
                    bin = 0
                elif 0.005<=value<0.01:
                    bin = 1
                elif 0.01<=value<0.02:
                    bin = 2
                elif 0.02<=value<0.05:
                    bin = 3
                elif 0.05<=value<0.1:
                    bin = 4
                elif 0.1<=value<0.2:
                    bin = 5
                elif 0.2<=value<0.5:
                    bin = 6
                elif 0.5<=value<1.:
                    bin = 7
                elif 1.<=value:
                    bin = 8
                if bin != 'nan':    
                    grid[index] = bin
    elif var == 'BA':
        title = 'Burnt Area'
        # Due to recurring fires in <year, not normalised.
        units = '(Fraction Burned per Year)'
        if binned:
            ticks=[0,.002,.005,.01,.02,.05,.1,.2,.5,'>1.0']
            bounds=[0,1,2,3,4,5,6,7,8,9]
            for index,value in np.ndenumerate(grid):
                bin = 'nan'
                if 0.<value<0.002:
                    bin = 0
                elif 0.002<=value<0.005:
                    bin = 1
                elif 0.005<=value<0.01:
                    bin = 2
                elif 0.01<=value<0.02:
                    bin = 3
                elif 0.02<=value<0.05:
                    bin = 4
                elif 0.05<=value<0.1:
                    bin = 5
                elif 0.1<=value<0.2:
                    bin = 6
                elif 0.2<=value<0.5:
                    bin = 7
                elif 0.5<=value:
                    bin = 8
                if bin != 'nan':    
                    grid[index] = bin
    
    
    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90, 
    urcrnrlon=180,urcrnrlat=90)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs=m.imshow(grid, interpolation='none')
    if binned:
        cb=m.colorbar(cs, "bottom", boundaries=bounds)
        cb.ax.set_xticklabels(ticks)
    else:
        cb=m.colorbar(cs, "bottom")
    cb.set_label(title + ' ' + units)
    plt.title('Mean '+title+' for '+str(year)+'-'+
        str(year+year_period)+', '+region_names[region]+', '+model)
    if save:
        return fig
    else:
        plt.show()


def plot_spatial_histogram(year, year_period, model, 
                var='FC', region=0, reg_type='boxes', save=False):
    """
    Year is in absolute terms, e.g. 1997.
    
    Takes argument model, which is a string that can be one of the
    following: 'gfed', 'jsbach', 'clm', 'ctem', 'blaze', 'orchidee',
    'inferno'
    
    Also takes argument var, which is a string that can take either
    'FC', 'BA', or 'emis' as an input, which will return the
    spatial histogram for fuel consumption, burnt area, or emissions
    respectively.
    
    To plot regional data using boxed regions (for no interpolation
    errors) leave the reg_type argument to the default 'boxes' and
    set the region argument from 1 to 12, which corresponds to 
    Boreal North America (BONA), Temperate North America (TENA),
    Equatorial Central and South America (EQCSA), Southernmost
    America (SOMA), Northern Europe (NOEU), Mediterranean
    and Middle East (MEME), Equatorial Africa (EQAF), Southernmost 
    Africa (SOAF), Boreal Asia (BOAS), Central Asia (CEAS),
    Equatorial Asia (EQAS), Australia/Oceania (AUST).
    
    To plot regional data using the GFED regions, set reg_type argument
    to 'gfed' and then set value of region argument from 1 to 14, 
    which correspond to BONA, TENA, CEAM, NHSA. SHSA, EURO, MIDE, 
    NHAF, SHAF, BOAS, CEAS, SEAS, EQAS, AUST respectively.
    For details regarding the GFED regions, go to 
    http://www.globalfiredata.org/data.html.
    """
    if reg_type=='boxes':
        region_names = ['Global','BONA','TENA','EQCSA','SOMA','NOEU',
                        'MEME','EQAF','SOAF','BOAS','CEAS','EQAS','AUST']
    elif reg_type=='gfed':
        region_names = ['Global','BONA','TENA','CEAM','NHSA','SHSA',
                        'EURO','MIDE','NHAF','SHAF','BOAS','CEAS',
                        'SEAS','EQAS','AUST']
                
    grid = load_var_grid(year, year_period, model, var)
    if region != 0:
        grid = get_regional_var_grid(year, year_period, region, 
                                    model, var, reg_type)
    
    if var == 'FC':
        title = 'Fuel Consumption'
        units = '($kg\, C\, m^{-2}\, burned$)'
    elif var == 'emis':
        title = 'Carbon Emissions'
        units = '($kg\, C\, m^{-2} \, year^{-1}$)'
    elif var == 'BA':
        title = 'Burnt Area'
        # Due to recurring fires in <year, not normalised.
        units = '(Fraction Burned per Year)'
    
    # Not needed, removing 0 values flattens.
    #flat_grid = np.ndarray.flatten(grid)
    flat_grid = grid[grid > 0]
    if var!='FC':
        flat_grid = np.divide(flat_grid, year_period)
    fig = plt.figure()
    plt.hist(flat_grid, 50, normed=False,
                facecolor='blue', alpha=0.75)
    
    plt.xlabel(title + ' ' + units)
    plt.ylabel('No. of Grid Cells')
    plt.title('Spatial Histogram of '+title+' for '+str(year)+
              '-'+str(year+year_period)+', '
              +region_names[region]+', '+model)
    if save:
        return fig
    else:
        plt.show()


def plot_multimodel_box(year, year_period, var='FC', 
             region=0, reg_type='boxes', model='all', save=False):
    """
    Year is given in absolute terms, e.g. 1997.
    
    Also takes argument var, which is a string that can take either
    'FC', 'BA', or 'emis' as an input, which will return the
    boxplots for fuel consumption, burnt area, or emissions
    respectively.
    
    To plot regional data using boxed regions (for no interpolation
    errors) leave the reg_type argument to the default 'boxes' and
    set the region argument from 1 to 12, which corresponds to 
    Boreal North America (BONA), Temperate North America (TENA),
    Equatorial Central and South America (EQCSA), Southernmost
    America (SOMA), Northern Europe (NOEU), Mediterranean
    and Middle East (MEME), Equatorial Africa (EQAF), Southernmost 
    Africa (SOAF), Boreal Asia (BOAS), Central Asia (CEAS),
    Equatorial Asia (EQAS), Australia/Oceania (AUST).
    
    To plot regional data using the GFED regions, set reg_type argument
    to 'gfed' and then set value of region argument from 1 to 14, 
    which correspond to BONA, TENA, CEAM, NHSA. SHSA, EURO, MIDE, 
    NHAF, SHAF, BOAS, CEAS, SEAS, EQAS, AUST respectively.
    For details regarding the GFED regions, go to 
    http://www.globalfiredata.org/data.html.
    """
    if reg_type=='boxes':
        region_names = ['Global','BONA','TENA','EQCSA','SOMA','NOEU',
                        'MEME','EQAF','SOAF','BOAS','CEAS','EQAS','AUST']
    elif reg_type=='gfed':
        region_names = ['Global','BONA','TENA','CEAM','NHSA','SHSA',
                        'EURO','MIDE','NHAF','SHAF','BOAS','CEAS',
                        'SEAS','EQAS','AUST']
                        
    model_list = ['gfed', 'jsbach', 'clm', 'ctem', 
                'blaze', 'orchidee', 'inferno']
    data=[]              
    if model=='all':     
        x_labels = model_list
        title_end = region_names[region]
        for model_name in model_list:
            grid = load_var_grid(year, year_period, model_name, var)
            if region != 0:
                grid = get_regional_var_grid(year, year_period, 
                                       region, model_name, var, reg_type)
            flat_grid = grid[grid > 0]
            if var!='FC':
                flat_grid = np.divide(flat_grid, year_period)
            data.append(flat_grid)
    else:
        x_labels = region_names[1:]
        title_end = model
        for region in range(1,len(region_names)):
            grid = get_regional_var_grid(year, year_period, 
                                       region, model, var, reg_type)
            flat_grid = grid[grid > 0]
            if var!='FC':
                flat_grid = np.divide(flat_grid, year_period)
            data.append(flat_grid)
      
    if var == 'FC':
        title = 'Fuel Consumption'
        units = '($kg\, C\, m^{-2}\, burned$)'
    elif var == 'emis':
        title = 'Carbon Emissions'
        units = '($kg\, C\, m^{-2} \, year^{-1}$)'
    elif var == 'BA':
        title = 'Burnt Area'
        # Due to recurring fires in <year, not normalised.
        units = '(Fraction Burned per Year)'
   
    fig = plt.figure()
    plt.boxplot(data, showmeans=True, whis=[5,95], sym='', 
                labels = x_labels)
    plt.ylabel(title +' '+ units)
    plt.title('Box plot of '+title+' for '+str(year)+
            '-'+str(year+year_period)+', '+title_end)
    if save:
        return fig
    else:
        plt.show()


#plot_map(1997,15,'inferno', 'FC', binned=True)
#plot_multimodel_box(1997,15,'FC',model='ctem')
#plot_spatial_histogram(1997, 15, 'ctem', var='FC')


#
# Spatial Comparison Maps and Correlations
#


def interp_GFED_func(lons, lats, var, year, year_period, method):
    """
    Interplates the GFED map for a given variable for a given
    year and year period using the given interploation method
    to the desired grid, determined by the lons and lats arguments.
    
    Used in interp_GFED_grid, see its associated docstring for
    more information on the use of the function.
    """
    GFED_data = load_var_grid(year,year_period,'gfed',var)
    
    
    lats_gfed = np.arange(-89.875, 90.,0.25)
    lons_gfed = np.arange(-179.875, 180.,0.25)
    lons_gfed, lats_gfed = np.meshgrid(lons_gfed, lats_gfed)
    
    new_grid = intrplt.griddata((lons_gfed.ravel(),lats_gfed.ravel()), 
                  GFED_data.ravel(), (lons,lats), method=method)
    return new_grid


def interp_GFED_grid(year, year_period, model, var='FC',
                    method='nearest', plot=False):
    """
    Returns the GFED grid for a given variable for a given year
    and year period, interpolated to the desired model's grid.
    
    Year is in absolute terms, e.g. 1997.
    
    Takes argument model, which is a string that can be one of the
    following: 'gfed', 'jsbach', 'clm', 'ctem', 'blaze', 'orchidee',
    'inferno'
    
    Argument var can take values 'FC', 'emis', 'BA', to get the
    grid for fuel consumption, burnt area, or emissions respectively.
    
    Argument method determines the interpolation method, and can
    take values 'nearest', 'linear', or 'cubic' for the associated
    interpolation methods (self explanatory. Default is set to 
    nearest neighbour, as it is the fastest and simplest, without
    too much loss in accuracy.
    
    Argument plot can be set to True to show map of interpolated
    data, useful for checks.
    """
    if model=='gfed':
        lats = np.arange(-89.875, 90.,0.25)
        lons = np.arange(-179.875, 180.,0.25)
    elif model=='jsbach':
        lats = jsbach.grid_JSBACH["latitude"]
        lats = np.array(lats)
        lats = lats[::-1]
        lons = jsbach.grid_JSBACH["longitude"]
        lons = np.array(lons) - np.max(lons)/2.
    elif model=='clm':
        lats = clm.grid_CLM["lat"]
        lons = clm.grid_CLM["lon"]
        lons = np.array(lons) - np.max(lons)/2.
        
    elif model=='ctem':
        lats = ctem.grid_CTEM["lat"]
        lons = ctem.grid_CTEM["lon"]
        lons = np.array(lons) - np.max(lons)/2.
    elif model=='blaze':
        lats = blaze.grid_BLAZE["lat"]
        lons = blaze.grid_BLAZE["lon"]
    elif model=='orchidee':
        lats = orchidee.grid_ORCHIDEE["latitude"]
        lats = np.array(lats)
        lats = lats[::-1]
        lons = orchidee.grid_ORCHIDEE["longitude"]
    elif model=='inferno':
        lats = inferno.grid_INFERNO["latitude"]
        lons = inferno.grid_INFERNO["longitude"]  
        lons = np.array(lons) - np.max(lons)/2.
      
    lons, lats = np.meshgrid(lons, lats)
    
    new_grid = interp_GFED_func(lons,lats,var,year,
                                year_period,method)

    if plot:
        fig=plt.figure()
        m = Basemap(llcrnrlon=-180,llcrnrlat=-90, 
            urcrnrlon=180,urcrnrlat=90)
        m.drawcoastlines()
        m.drawparallels(np.arange(-90.,91.,30.))
        m.drawmeridians(np.arange(-180.,181.,60.))
        m.drawmapboundary(fill_color='white')
        m.imshow(new_grid,  interpolation='none')
        plt.title("GFED Data Interpolated for " + model)
        plt.show()
    else:
        return new_grid

def calc_spatial_correlation(year, year_period, model, var):
    """
    Calculates the Pearson correlation between a model and
    GFED for a given variable for a given year and year period.
    
    Function is used in get_spacial_correlations to get table
    of correlations for all models. See its associated docstring
    for more information.
    """
    GFED_grid = interp_GFED_grid(year,year_period,model,var)
    model_grid = load_var_grid(year,year_period,model,var)
    
    GFED_flat = GFED_grid.flatten()
    model_flat = model_grid.flatten()
    
    pearson = stats.pearsonr(GFED_flat, model_flat)
    return pearson
    
    
def get_spatial_correlations(year, year_period, var):
    """
    Calculates and returns a table of Pearson correlations
    of the given variable var (which takes values 'FC', 'emis',
    and 'BA' for fuel consumption, emissions, and burnt area
    respectively) for a given year and year period of all the
    models' maps for that period with GFED's map.
    
    The correlation is calculated by interpolating the GFED
    variable map to each individual model's grid and then
    flattening both grids and calculating their Pearson
    correlation. Should give a measure of the spatial
    accuracy of the models using GFED as the performance
    metric/observations.
    """
    model_list = ['gfed', 'jsbach', 'clm', 'ctem', 
                'blaze', 'orchidee', 'inferno']
    pearson_list = []
    for model in model_list[1:]:
        pearson = calc_spatial_correlation(year,year_period,model,var)
        pearson_list.append(pearson)
    
    print('Table of Correlations for '+var+' for '+
            str(year)+'-'+str(year+year_period))
    print('===========================================')
    for i in range(len(pearson_list)):
        print(model_list[i+1]+': '+str(pearson_list[i]))
    return pearson_list


def plot_diff_map(year, year_period, model, var, 
                    binned=True, save=False):
    """
    Plots a map of the relative difference of the given model
    with GFED for a given variable.
    
    Year is in absolute terms, e.g. 1997.
    
    Takes argument model, which is a string that can be one of the
    following: 'gfed', 'jsbach', 'clm', 'ctem', 'blaze', 'orchidee',
    'inferno'
    
    Also takes argument var, which is a string that can take either
    'FC', 'BA', or 'emis' as an input, which will return the
    difference map for fuel consumption, burnt area, or emissions
    respectively.
    
    Argument binned is set by default to true. If set to false,
    it will arrange the colorbar and colormap automatically.
    It is suggested that bins are used for standardised maps
    and easier comparison.
    """
    GFED_grid = interp_GFED_grid(year,year_period,model,var,
                                 method='linear')
    GFED_grid = np.array(GFED_grid)
    model_grid = load_var_grid(year,year_period,model,var)
    model_grid = np.array(model_grid)
    
    diff_grid = np.divide(model_grid-GFED_grid,GFED_grid)
    diff_grid = np.multiply(diff_grid, 100)
    
    if binned:
            ticks=[-100,-50,0,50,100,150,200,300,'>500']
            bounds=[-120,-60,0,20,40,60,80,100,120]
            for index,value in np.ndenumerate(diff_grid):
                bin = 'nan'
                if -100<=value<-50:
                    bin = bounds[0]
                elif -50<=value<0:
                    bin = bounds[1]
                elif 0<=value<50:
                    bin = bounds[2]
                elif 50<=value<100:
                    bin = bounds[3]
                elif 100<=value<150:
                    bin = bounds[4]
                elif 150<=value<200:
                    bin = bounds[5]
                elif 200<=value<300:
                    bin = bounds[6]
                elif 300<=value:
                    bin = bounds[7]
                if bin != 'nan':    
                    diff_grid[index] = bin
    
    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90, 
        urcrnrlon=180,urcrnrlat=90)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs=m.imshow(diff_grid,  interpolation='none',
                cmap=plt.cm.coolwarm)
    plt.title("Relative Difference in "+var+
                ", GFED vs "+ model)
    if binned:
        cb=m.colorbar(cs, "bottom", boundaries=bounds)
        cb.ax.set_xticklabels(ticks)
    else:
        cb=m.colorbar(cs, "bottom")
    cb.set_label('Relative Difference (%)')
    if save:
        return fig
    else:
        plt.show()


def interp_std_func(model, var, year, year_period, 
                    method, ref_grid):
    """
    Interplates the GFED map for a given variable for a given
    year and year period using the given interploation method
    to the desired grid, determined by the lons and lats arguments.
    
    Used in interp_GFED_grid, see its associated docstring for
    more information on the use of the function.
    """
    model_data = load_var_grid(year,year_period,model,var)
    lons, lats = get_lons_lats(model)
      
    lons, lats = np.meshgrid(lons, lats)
    
    if ref_grid=='ctem':
        lats_ref = ctem.grid_CTEM["lat"]
        lons_ref = ctem.grid_CTEM["lon"]
        lons_ref = np.array(lons_ref) - np.max(lons_ref)/2.
    
    if ref_grid=='gfed':
        lats_ref = np.arange(-89.875, 90.,0.25)
        lons_ref = np.arange(-179.875, 180.,0.25)
    
    lons_ref, lats_ref = np.meshgrid(lons_ref, lats_ref)
    
    new_grid = intrplt.griddata((lons.ravel(),lats.ravel()), 
                  model_data.ravel(), (lons_ref,lats_ref), method=method)
    return new_grid
    
    
def plot_std_map(year, year_period, var='FC', method='nearest',
                 ref_grid='gfed', binned='True'):
    model_list = ['jsbach', 'clm', 'blaze', 
                    'orchidee', 'inferno','ctem']
    grid_list = []
    if ref_grid=='gfed':
        print 'Using GFED resolution.'
    if ref_grid=='ctem':
        print 'Using CTEM resolution.'
        model_list=model_list[:-1]
        model_grid = load_var_grid(year,year_period,'ctem',var)
        model_grid = np.array(model_grid)
        grid_list.append(model_grid)
    
    for model in model_list:
        model_grid=interp_std_func(model, var, year, year_period, 
                                      method, ref_grid)
        model_grid = np.array(model_grid)
        grid_list.append(model_grid)
        
    grid_list = np.array(grid_list)
    std_grid = np.std(grid_list, axis=0)
    std_grid[std_grid==0]=np.nan

    
    
    if var == 'FC':
        title = 'Fuel Consumption'
        units = '($kg\, C\, m^{-2}\, burned$)'
        if binned:
            ticks=[0,0.1,0.5,1,2,3,4,5,10,'>20.0']
            bounds=[0,1,2,3,4,5,6,7,8,9]
            for index,value in np.ndenumerate(std_grid):
                bin = 'nan'
                if 0.<value<0.1:
                    bin = 0
                elif 0.1<=value<0.5:
                    bin = 1
                elif 0.5<=value<1:
                    bin = 2
                elif 1<=value<2:
                    bin = 3
                elif 2<=value<3:
                    bin = 4
                elif 3<=value<4:
                    bin = 5
                elif 4<=value<5:
                    bin = 6
                elif 5<=value<10:
                    bin = 7
                elif 10<=value:
                    bin = 8
                if bin != 'nan':    
                    std_grid[index] = bin
    elif var == 'emis':
        title = 'Carbon Emissions'
        units = '($kg\, C\, m^{-2} \, year^{-1}$)'
    elif var == 'BA':
        title = 'Burnt Area'
        # Due to recurring fires in <year, not normalised.
        units = '(Fraction Burned per Year)'
    
    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90, 
        urcrnrlon=180,urcrnrlat=90)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs=m.imshow(std_grid,  interpolation='none')
    plt.title("Multimodel Standard Deviations for "+title+
                ", "+str(year)+'-'+str(year+year_period))
    if binned:
        cb=m.colorbar(cs, "bottom", boundaries=bounds)
        cb.ax.set_xticklabels(ticks)
    else:
        cb=m.colorbar(cs, "bottom")
    cb.set_label('Standard Deviation '+units)
    plt.show()



#get_spatial_correlations(1997,15,'FC')
#plot_diff_map(1997,15,'clm','FC')
#plot_std_map(1997,15,ref_grid='ctem')

