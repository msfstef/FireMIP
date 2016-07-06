from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm
                    
emis_Inferno = Dataset('../../model_data/Inferno_S1_fFirepft.nc4', 
                'r', format = 'NETCDF4')
BA_Inferno = Dataset('../../model_data/Inferno_S1_burntArea.nc4', 
                'r', format = 'NETCDF4')              
grid_Inferno = Dataset('../../model_data/Inferno_grid.nc', 
                    'r', format = 'NETCDF4')
landmask_Inferno = Dataset('../../model_data/CRU-NCEP-LandMask.nc',
                    'r', format = 'NETCDF4')
landCover_Inferno = Dataset('../../model_data/Inferno_S1_LandCoverFrac.nc4',
                    'r', format = 'NETCDF4')

#a=landmask_Inferno["lsm"]                    
#a=landCover_Inferno["LandCoverFrac"][300*12:301*12]
#a=BA_Inferno["burntArea"][300*12:301*12]
#a = np.array(a)
#a[a>10**10]=0
#print len(a[a>0.7])
#print len(a[a>0.1])
#print len(a[a>=0])
#print len(a[a<0])
#print a.size

                    

#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year_start, month_period, emis_data, BA_data, landmask, landCover_data):
    time = int(year_start*12)
    sec_per_month = []
    for i in range(month_period):
        if time+i == 0:
            sec_per_month.append(emis_data["time"][0])
        else:
            sec_per_month.append(emis_data["time"][time+i]-emis_data["time"][time+i-1])
    sec_per_month= np.array(sec_per_month)
    
    # Ignore division by zero warning. Returns NaN.
    np.seterr(divide='ignore')
    
    # Remove fill values.
    landmask = np.array(landmask["lsm"])
    landmask[landmask > 1.01] = 0

    actual_landCover_data = np.multiply(landCover_data["LandCoverFrac"][time:time+month_period, :9],
                                         landmask)
    actual_emission_data= np.multiply(actual_landCover_data, emis_data["fFirepft"][time:time+month_period])
    
    # Possibly convert to decimals and multiply by landCover (?)
    BA = BA_data["burntArea"][time:time+month_period]
    BA[BA<0]=0
    actual_BA = np.divide(BA, 100)
    actual_BA = np.multiply(actual_BA, actual_landCover_data)
    inverse_BA = 1./actual_BA
    
    # Remove infinities due to division by 0.
    inverse_BA[inverse_BA == np.inf] = 0
    
    FC_data_pft = np.multiply(actual_emission_data, inverse_BA)
    FC_per_month_pft = np.multiply(FC_data_pft, 
                sec_per_month[:, np.newaxis, np.newaxis, np.newaxis])
    
    FC_per_month = np.sum(FC_per_month_pft, axis = 1)
    fuel_consumption = np.sum(FC_per_month, axis = 0)
    return np.sum(np.sum(actual_BA, axis = 0), axis = 0)
    
    
def get_global_total_FC_yearly(year, emis_data, BA_data, landmask, landCover_data):
    FC_grid = get_grid_fuel_consumption(year, 12, emis_data, BA_data, landmask, landCover_data)
    fuel_consumption = np.sum(FC_grid)
    return fuel_consumption
    
    
def get_global_mean_FC_yearly(year, emis_data, BA_data, grid_data, landCover_data, time_data):
    time = int(year*12)   
    FC_data = get_grid_fuel_consumption(year,12,emis_data, BA_data, landCover_data, time_data) 
    
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["cell_area"])
    return global_mean_FC 
    
    
def plot_global_mean_FC_yearly(no_years, emis_data, BA_data, grid_data, landCover_data, time_data):
    years = int(len(time_data["time"])/12)
    x_data = range(years-1)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_mean_FC_yearly(x,emis_data,BA_data,grid_data,landCover_data,time_data))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='CLM Results')
    plt.ylabel('Fuel Consumption ($kg\, C / m^2 \, burned$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
    
    
def plot_FC_map_period(year_start, year_period, emis_data, BA_data, landmask, landCover_data):
    month_period = int(year_period*12)
    map_data = get_grid_fuel_consumption(year_start,month_period,emis_data, BA_data, landCover_data, time_data)    
    
    lats = landCover_data["latitude"]
    lons = landCover_data["longitude"]
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
    cb.set_label("Carbon emitted per area burned ($kg\, C / m^2 \, burned$)")
    plt.title("Fuel Consumption, CLM")
    plt.show()
    
#print get_global_FC_yearly(300, emis_CLM, BA_CLM, landCover_CLM, time_data)
#plot_FC_map_period(300,6 , emis_Inferno, BA_Inferno, landmask_Inferno, landCover_Inferno)
#plot_global_FC_yearly(30, emis_CLM, BA_CLM, landCover_CLM, time_data)
#plot_global_mean_FC_yearly(200, emis_Inferno, BA_Inferno, grid_Inferno, landmask_Inferno, landCover_Inferno)
#print get_global_mean_FC_yearly(300,emis_Inferno,BA_Inferno, grid_Inferno, landmask_Inferno, landCover_Inferno)


#
# Carbon Emissions Analysis
#
                    
def get_grid_emissions(year_start, month_period, 
                    emis_data, grid_data, landmask, landCover_data):
    time = int(year_start*12)
    sec_per_month = []
    for i in range(month_period):
        if time+i == 0:
            sec_per_month.append(emis_data["time"][0])
        else:
            sec_per_month.append(emis_data["time"][time+i]-emis_data["time"][time+i-1])
    sec_per_month= np.array(sec_per_month)
    
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    # Remove fill values.    
    landmask = np.array(landmask["lsm"])
    landmask[landmask > 1.01] = 0
    actual_grid_area = np.multiply(grid_data["cell_area"],landmask)
    
    landCover = landCover_data["LandCoverFrac"][time:time + month_period,:9]
    complete_area_data = np.multiply(landCover, actual_grid_area)
    
    emis_rate_data = np.multiply(emis_data["fFirepft"][time:time + month_period], complete_area_data)
    emis_per_month = np.multiply(emis_rate_data, 
                sec_per_month[:, np.newaxis, np.newaxis, np.newaxis])
    emissions_pft = np.sum(emis_per_month, axis = 0)
    emissions = np.sum(emissions_pft, axis = 0)
    return emissions
    
def get_global_emissions_yearly(year, emis_data, grid_data, landmask, landCover_data):
    emissions_grid = get_grid_emissions(year, 12, emis_data, grid_data, landmask, landCover_data)
    emissions = np.sum(emissions_grid)
    return emissions

def get_global_emissions_monthly(month, emis_data, grid_data, landmask, landCover_data):
    emissions_grid = get_grid_emissions(month/12., 1, emis_data, grid_data, landmask, landCover_data)
    emissions = np.sum(emissions_grid)
    return emissions
   
   
    
def plot_global_emissions_monthly(no_months, emis_data, grid_data):
    x_data = range(len(emis_data["time"])-1)
    y_data = []
    for x in x_data[-no_months:]:
        print("%.2f" % ((x-x_data[-no_months])/
                float(len(x_data[-no_months:]))))
        y_data.append(get_global_emissions_monthly(x, emis_data, grid_data)/(10**12))
    plt.plot(x_data[-no_months:], y_data, color='r', linewidth=2.0,
               label='Inferno Results')
    plt.ylabel('Carbon Emitted ($Pg/month$)')
    plt.xlabel('Month')
    plt.legend()
    plt.show()
    print "success"
   

def plot_global_emissions_yearly(no_years, emis_data, grid_data, landmask, landCover_data):
    years = int(len(emis_data["time"])/12 +1)
    x_data = range(years-1)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        # Progress bar
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_emissions_yearly(x, emis_data, grid_data, landmask, landCover_data)/(10**12))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='Inferno Results')
    plt.ylabel('Carbon Emitted ($Pg/year$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
    print "success"

def plot_map_period(year_start, year_period, emis_data, grid_data, landmask, landCover_data):
    month_period = int(year_period*12)
    map_data = get_grid_emissions(year_start,month_period,emis_data,grid_data,landmask, landCover_data)
    # Convert to billions of kg.
    map_data = np.multiply(map_data, (1./(10**9)))
    
    
    lats = emis_data["latitude"]
    lons = emis_data["longitude"]
    
    lons, lats = np.meshgrid(lons, lats)
    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-60, urcrnrlon=180,urcrnrlat=80,projection='mill')
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='aqua')
    cs = m.contourf(lons,lats, map_data, 100, cmap=plt.cm.YlOrRd, latlon=True)
    cb = m.colorbar(cs, "bottom", size="5%", pad="2%")
    cb.set_label("Billions of kg of Carbon")
    plt.title("Total Emissions, Inferno")
    plt.show()
        

#print emis_Inferno["longitude"][-1]
#print grid_Inferno["longitude"][-1]

#plot_map_period(178,1,emis_Inferno, grid_Inferno, landmask_Inferno, landCover_Inferno)
#print get_grid_emissions(300,12,data_JSBACH,grid_data_JSBACH)

#plot_global_emissions_yearly(200,emis_Inferno, grid_Inferno, landmask_Inferno, landCover_Inferno)
#print get_global_emissions_yearly(178, emis_Inferno, grid_Inferno, landmask_Inferno, landCover_Inferno)
