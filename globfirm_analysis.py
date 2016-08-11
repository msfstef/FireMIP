from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm

emis_GLOBFIRM = Dataset('../../model_data/LPJ-GUESS-globfirm_SF1_Cfire.nc','r', format = 'NETCDF4')
BA_GLOBFIRM = Dataset('../../model_data/LPJ-GUESS-globfirm_SF1_burntArea.nc', 'r', format = 'NETCDF4')
grid_GLOBFIRM = Dataset('../../model_data/HalfDegree-gridarea-8975.nc',
                        'r', format = 'NETCDF4')

#
# Burnt Area Analysis
#

def get_grid_burnt_area(year, year_period, BA_data, grid_data,
                         keep_time=False):
    BA = BA_data["burntArea."][year:year+year_period]
    BA = np.divide(BA, 100.)
    
    burnt_area_data = np.multiply(BA, grid_data["cell_area"])
    burnt_area_data = np.array(burnt_area_data)
    burnt_area_data[burnt_area_data<0.]=0
    if keep_time:
        return burnt_area_data
    burnt_area_data = np.sum(burnt_area_data, axis=0)
    return burnt_area_data
    
def get_global_BA_yearly(year, BA_data, grid_data):
    BA_grid = get_grid_burnt_area(year, 1, BA_data, grid_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area
    
    
def plot_global_BA_yearly(no_years, BA_data, grid_data):
    years = len(BA_data["time"])
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
               label='LPJ BLAZE Results')
    plt.ylabel('Burnt Area (millions of $km^2$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

def plot_map_period(year_start, year_period, BA_data, grid_data):
    map_data = get_grid_burnt_area(year_start,year_period, BA_data, grid_data)
    map_data = np.divide(map_data, year_period)
    map_data = np.divide(map_data, grid_data["cell_area"])
    map_data[map_data==0]=np.nan
    
    lats = grid_data["latitude"]
    lons = grid_data["longitude"]
    lons, lats = np.meshgrid(lons, lats)

    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-90, 
    urcrnrlon=180,urcrnrlat=90)
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs=m.imshow(map_data, interpolation='none')
    if False:
        cb=m.colorbar(cs, "bottom", boundaries=bounds)
        cb.ax.set_xticklabels(ticks)
    else:
        cb=m.colorbar(cs, "bottom")
    cb.set_label("Fraction Burned")
    plt.title("Total Burnt Area 1997-2012, GlobFIRM")
    plt.show()


#plot_map_period(297,16,BA_GLOBFIRM, grid_GLOBFIRM)
#plot_global_BA_yearly(16,BA_GLOBFIRM, grid_GLOBFIRM)
#print get_global_BA_yearly(300, BA_GLOBFIRM, grid_GLOBFIRM)


#
# Carbon Emissions Analysis
#
                   
def get_grid_emissions(year_start, year_period, emis_data, grid_data,
                        keep_time=False):
    # Using 1 year equal to 365.25 days.
    sec_per_year = 31557600.
    
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    emis = np.array(emis_data["fFire."][year_start:year_start+year_period])
    emis = np.multiply(emis, grid_data["cell_area"])
    emis = np.multiply(emis, sec_per_year)
    emis[emis<0]=0.
    
    if keep_time:
        return emissions
    emissions = np.sum(emissions, axis = 0)
    return emissions

def get_global_emissions_yearly(year, emis_data, grid_data):
    emissions_grid = get_grid_emissions(year, 1, emis_data, grid_data)
    emissions = np.sum(emissions_grid)
    return emissions


def plot_global_emissions_yearly(no_years, emis_data, grid_data):
    years = len(emis_data["time"])
    print years
    x_data = range(years)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_emissions_yearly(x, emis_data, grid_data)/(10**12))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='LPJ GlobFIRM Results')
    plt.ylabel('Carbon Emitted ($Pg/year$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

def plot_map_period(year_start, year_period, emis_data, grid_data):
    map_data = get_grid_emissions(year_start,year_period, emis_data, grid_data)
    map_data = np.divide(map_data, year_period)
    map_data[map_data==0]=np.nan
    
    lats = grid_data["latitude"]
    lons = grid_data["longitude"]
    lons, lats = np.meshgrid(lons, lats)

    binned=True
    if binned:
            map_data=np.divide(map_data,grid_data["cell_area"])
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
    plt.title("Total Emissions 1997-2012, GlobFIRM")
    plt.show()

#plot_map_period(297,16,emis_GLOBFIRM, grid_GLOBFIRM)    
#plot_global_emissions_yearly(20, emis_GLOBFIRM, grid_GLOBFIRM)
#print get_global_emissions_yearly(60,emis_GLOBFIRM,grid_GLOBFIRM)
