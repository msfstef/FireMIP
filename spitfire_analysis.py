from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm


emis_SPITFIRE = Dataset('../../model_data/LPJ-GUESS-SPITFIRE_SF1_fFirepft.nc', 
                'r', format = 'NETCDF4')
BA_SPITFIRE = Dataset('../../model_data/LPJ-GUESS-SPITFIRE_SF1_burntArea.nc', 
                'r', format = 'NETCDF4')            
grid_SPITFIRE = Dataset('../../model_data/HalfDegree-gridarea-8975.nc', 
                    'r', format = 'NETCDF4')                   


#
# Burnt Area Analysis
#

def get_grid_burnt_area(year, month_period, BA_data, grid_data):
    BA = BA_data["burntArea"][year*12:year*12+month_period]
    BA = np.divide(BA, 100.)
    
    burnt_area_data = np.multiply(BA, grid_data["cell_area"])
    burnt_area_data = np.sum(burnt_area_data, axis=0)
    return burnt_area_data
    
def get_global_BA_yearly(year, BA_data, grid_data):
    BA_grid = get_grid_burnt_area(year, 12, BA_data, grid_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area
    
    
def plot_global_BA_yearly(no_years, BA_data, grid_data):
    years = int(len(BA_data["time"])/12)
    x_data = range(years)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_BA_yearly(x,BA_data,grid_data))
    # Convert to millions of km^2.
    y_data = np.array(y_data)
    y_data = np.divide(y_data, 1e12)
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='LPJ SPITFIRE Results')
    plt.ylabel('Burnt Area (millions of $km^2$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

#plot_global_BA_yearly(16,BA_BLAZE, grid_BLAZE)
#print get_global_BA_yearly(300, BA_BLAZE, grid_BLAZE)

#
# Carbon Emissions Analysis
#


def get_grid_emissions(year_start, month_period, emis_data, grid_data):
    time = int(year_start*12)
    
    # Original output is in 'per month' units, so using
    # exact conversion unit to get correct results.
    sec_per_month = 1/0.000000388024691
    
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    
    emis = emis_data["fFirepft"][time:time+month_period]
     
    actual_emis_data = np.sum(emis, axis=1) 
    emis_rate_data = np.multiply(actual_emis_data, 
                                grid_data["cell_area"])
    emis_per_month = np.multiply(emis_rate_data,sec_per_month)
    emissions = np.sum(emis_per_month, axis = 0)
    return emissions

def get_global_emissions_yearly(year, emis_data, grid_data):
    emissions_grid = get_grid_emissions(year, 12, emis_data, grid_data)
    emissions = np.sum(emissions_grid)
    return emissions
   

def plot_global_emissions_yearly(no_years, emis_data, grid_data):
    years = int(len(emis_data["time"])/12)
    x_data = range(years-1)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_emissions_yearly(x, emis_data, grid_data)/(10**12))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='LPJ SPITFIRE Results')
    plt.ylabel('Carbon Emitted ($Pg/year$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

def plot_map_period(year_start, year_period, emis_data, grid_data):
    month_period = int(year_period*12)
    map_data = get_grid_emissions(year_start,month_period, 
                emis_data, grid_data)
    map_data = np.divide(map_data,grid_data["cell_area"])
    map_data = np.divide(map_data, year_period)
    map_data[map_data==0]=np.nan
    
    lats = grid_data["latitude"]
    lons = grid_data["longitude"]
    lons, lats = np.meshgrid(lons, lats)

    binned=True
    if binned:
            ticks=[0.,.005,.01,.02,.05,.1,.2,.5,1.,'>5.0']
            bounds=[0,1,2,3,4,5,6,7,8,9]
            for index,value in np.ndenumerate(map_data):
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
                    map_data[index] = bin
    
    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90, 
    urcrnrlon=180,urcrnrlat=90)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs=m.imshow(map_data, interpolation='none')
    if binned:
        cb=m.colorbar(cs, "bottom", boundaries=bounds)
        cb.ax.set_xticklabels(ticks)
    else:
        cb=m.colorbar(cs, "bottom")
    cb.set_label("kg C per m^2")
    plt.title("Total Emissions 1997-2012, SPITFIRE")
    plt.show()
    
 
#plot_global_emissions_yearly(20,emis_SPITFIRE, grid_SPITFIRE, landCover_SPITFIRE)
#print get_global_emissions_yearly(312,emis_SPITFIRE,grid_SPITFIRE, landCover_SPITFIRE)
#plot_map_period(297,15,emis_SPITFIRE,grid_SPITFIRE, landCover_SPITFIRE)


#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year_start, month_period, emis_data,
                 BA_data, landCover_data, monthly=False):
    time = int(year_start*12)
    
    # Original output is in 'per month' units, so using
    # exact conversion unit to get correct results.
    sec_per_month = 1/0.000000388024691
    
    # Ignore division by zero warning. Returns NaN.
    np.seterr(divide='ignore')
    
    BA = BA_data["burntArea"][time:time+month_period]
    fractional_BA = np.divide(BA, 100)
    fractional_BA = np.sum(fractional_BA, axis=1)
    if not monthly:
        fractional_BA = np.sum(fractional_BA, axis=0)
    inverse_BA = 1./np.array(fractional_BA)
    # Remove infinities due to division by 0.
    inverse_BA[inverse_BA == np.inf] = 0
    
    emis = emis_data["fFirepft"][time:time+month_period]
    emis = np.multiply(emis, sec_per_month)
    if not monthly:
        emis = np.sum(emis, axis=0)
        
    fuel_consumption = np.multiply(emis,inverse_BA)
    if monthly:
        fuel_consumption = np.nansum(fuel_consumption, axis = 0)
    return fuel_consumption

    
def get_global_mean_FC_yearly(year, emis_data, BA_data, grid_data):
    time = int(year*12)   
    FC_data = get_grid_fuel_consumption(year,12,emis_data, BA_data)
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["cell_area"])
    return global_mean_FC 
    
    
def get_global_mean_FC_yearly_rough(year, emis_data, BA_data, grid_data):
    total_emis = get_global_emissions_yearly(year,emis_data,grid_data)
    total_BA = get_global_BA_yearly(year,BA_data,grid_data)
    global_mean_FC = total_emis/total_BA
    return global_mean_FC

    
def plot_global_mean_FC_yearly(no_years, emis_data, BA_data, grid_data):
    years = int(len(time_data["time"])/12)
    x_data = range(years-1)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_mean_FC_yearly(x, emis_data, BA_data, grid_data))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='LPJ BLAZE Results')
    plt.ylabel('Fuel Consumption ($kg\, C/m^2 \,burned$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
    
#print get_global_FC_yearly(300, emis_CLM, BA_CLM, landCover_CLM)
#plot_FC_map_period(140,1 , emis_CTEM,BA_CTEM,landCover_CTEM)
#plot_global_mean_FC_yearly(20, emis_CTEM,BA_CTEM,grid_CTEM, landCover_CTEM)
#print get_global_mean_FC_yearly(150,emis_CTEM,BA_CTEM,grid_CTEM, landCover_CTEM)






