from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm


emis_ORCHIDEE = Dataset('../../model_data/ORCHIDEE_SF1_fFirepft.nc', 
                'r', format = 'NETCDF4')             
BA_ORCHIDEE = Dataset('../../model_data/ORCHIDEE_SF1_burntArea.nc', 
                'r', format = 'NETCDF4')                              
grid_ORCHIDEE = Dataset('../../model_data/HalfDegree-gridarea-8975-inverted.nc', 
                    'r', format = 'NETCDF4')
landCover_ORCHIDEE = Dataset('../../model_data/ORCHIDEE_SF1_landCoverFrac.nc', 
                'r', format = 'NETCDF4')                                  
                    
# ORCHIDEE has no useful separate time data, must use JSBACH data which matches it.             
time_data = Dataset('../../model_data/JSBACH_SF1_fFirepft.nc', 
                'r', format = 'NETCDF4')



#
# Burnt Area Analysis
#

def get_grid_burnt_area(year, month_period, BA_data, grid_data, landCover_data):
    # Remove meaningless values (fill and small negative), and make decimal.
    landCover = landCover_data["landCoverFrac"][year:int(year+month_period/12)]
    landCover[landCover>10**10] = 0
    landCover[landCover<0] = 0
    landCover = np.divide(landCover,100)

    BA = BA_data["burntArea"][year*12:year*12+month_period]
    BA = np.multiply(landCover, BA)
    BA[BA<0.] = 0
    BA[BA>1.] = 0
    
    burnt_area_data = np.multiply(BA, grid_data["cell_area"])
    burnt_area_data = np.sum(burnt_area_data, axis=0)
    burnt_area_data = np.sum(burnt_area_data, axis=0)
    return burnt_area_data
    
def get_global_BA_yearly(year, BA_data, grid_data, landCover_data):
    BA_grid = get_grid_burnt_area(year, 12, BA_data, grid_data, landCover_data)
    burnt_area = np.sum(BA_grid)
    return burnt_area
    
    
def plot_global_BA_yearly(no_years, BA_data, grid_data, landCover_data):
    years = int(len(BA_data["time"])/12)
    x_data = range(years)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_BA_yearly(x,BA_data,grid_data,landCover_data))
    # Convert to millions of km^2.
    y_data = np.array(y_data)
    y_data = np.divide(y_data, 1e12)
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='ORCHIDEE Results')
    plt.ylabel('Burnt Area (millions of $km^2$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

def plot_burnt_area_map_period(year_start, year_period, BA_data, grid_data, landCover_data):
    BA_list=[]                        
    for i in range(year_period):
        BA_list.append(get_grid_burnt_area(year_start+i,12,BA_data,grid_data,landCover_data))
    map_data = np.sum(BA_list, axis=0)
    # Change units.
    map_data = np.divide(map_data, 1e9)
    # Replace 0 by nan, take mean.
    map_data[map_data==0]=np.nan
    map_data = np.divide(map_data,year_period)
    
    
    lats = BA_data["latitude"]
    lons = BA_data["longitude"]
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
    plt.title("Mean Burnt Area, ORCHIDEE")
    plt.show()
 
    
#plot_global_BA_yearly(16,BA_ORCHIDEE, grid_ORCHIDEE, landCover_ORCHIDEE)
#plot_burnt_area_map_period(297,16,BA_ORCHIDEE, grid_ORCHIDEE, landCover_ORCHIDEE)    
#print get_global_BA_yearly(300, BA_ORCHIDEE, grid_ORCHIDEE, landCover_ORCHIDEE)


#
# Carbon Emissions Analysis
#

                
def get_grid_emissions(year_start, month_period, emis_data, grid_data, landCover_data, time_data):
    time = int(year_start*12)
    days_per_month = []
    if (time+12) == len(time_data["time"]):
        for i in range(month_period):
            if time+i+1 < len(time_data["time"]):
                days_per_month.append(time_data["time"][time+i+1]-time_data["time"][time+i])
            else:
                days_per_month.append(31)
    else:
        for i in range(month_period):
            days_per_month.append(time_data["time"][time+i+1]-time_data["time"][time+i])
    sec_per_month = np.multiply(days_per_month, 86400)
    
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    # Remove meaningless values (fill and small negative), and make decimal.
    landCover = landCover_data["landCoverFrac"][year_start:int(year_start+month_period/12.)]
    landCover[landCover>10**10] = 0
    landCover[landCover<0] = 0
    landCover = np.divide(landCover,100)
       
    complete_area_data = np.multiply(landCover, grid_data["cell_area"])
    emis_rate_data = np.multiply(emis_data["fFirepft"][time:time+month_period], complete_area_data)
    emis_per_month = np.multiply(emis_rate_data, 
                sec_per_month[:, np.newaxis, np.newaxis, np.newaxis])
    emissions_pft = np.sum(emis_per_month, axis = 0)
    emissions = np.sum(emissions_pft, axis = 0)
    return emissions

def get_global_emissions_yearly(year, emis_data, grid_data, landCover_data, time_data):
    emissions_grid = get_grid_emissions(year, 12, emis_data, grid_data, landCover_data, time_data)
    emissions = np.sum(emissions_grid)
    return emissions
 

def plot_global_emissions_yearly(no_years, emis_data, grid_data, landCover_data, time_data):
    years = int(len(time_data["time"])/12)
    x_data = range(years-1)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_emissions_yearly(x, emis_data, grid_data, landCover_data, time_data)/(10**12))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='ORCHIDEE Results')
    plt.ylabel('Carbon Emitted ($Pg/year$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()

def plot_map_period(year_start, year_period, emis_data, grid_data, landCover_data, time_data):
    yearly_data = []
    for i in range(year_period):
        yearly_data.append(get_grid_emissions(year_start+i,12, emis_data, grid_data, landCover_data, time_data))
    map_data = np.mean(yearly_data, axis=0)
    # Convert to billions of kg.
    map_data = np.multiply(map_data, (1./(10**9)))
    # Replace 0 by nan.
    map_data[map_data==0]=np.nan
    
    lats = grid_data["latitude"]
    lons = grid_data["longitude"]
    emis_data = map_data
    lons, lats = np.meshgrid(lons, lats)

    fig=plt.figure()
    m = Basemap(llcrnrlon=-180,llcrnrlat=-60, urcrnrlon=180,urcrnrlat=80,projection='mill')
    m.drawcoastlines()
    m.drawparallels(np.arange(-90.,91.,30.))
    m.drawmeridians(np.arange(-180.,181.,60.))
    m.drawmapboundary(fill_color='white')
    cs = m.contourf(lons,lats, emis_data, 100, cmap=plt.cm.YlOrRd, latlon=True)
    cb = m.colorbar(cs, "bottom", size="5%", pad="2%")
    cb.set_label("Billions of kg of Carbon")
    plt.title("Mean Emissions, ORCHIDEE")
    plt.show()
    
#print get_global_emissions_yearly(312,emis_ORCHIDEE,grid_ORCHIDEE, landCover_ORCHIDEE, time_data)
#plot_global_emissions_yearly(20,emis_ORCHIDEE,grid_ORCHIDEE, landCover_ORCHIDEE, time_data)
#plot_map_period(297,16,emis_ORCHIDEE,grid_ORCHIDEE, landCover_ORCHIDEE, time_data)


#
# Fuel Consumption Analysis
#

def get_grid_fuel_consumption(year_start, month_period, emis_data, BA_data, landCover_data, time_data):
    time = int(year_start*12)
    days_per_month = []
    if (time+12) == len(time_data["time"]):
        for i in range(month_period):
            if time+i+1 < len(time_data["time"]):
                days_per_month.append(time_data["time"][time+i+1]-
                                        time_data["time"][time+i])
            else:
                days_per_month.append(31)
    else:
        for i in range(month_period):
            days_per_month.append(time_data["time"][time+i+1]-time_data["time"][time+i])
    sec_per_month = np.multiply(days_per_month, 86400)
    
    # Ignore division by zero warning. Returns inf or NaN.
    np.seterr(divide='ignore')
    # Ignore overflow warning.
    np.seterr(over='ignore')
    
    # Remove meaningless values (fill and small negative), and make decimal.
    landCover = landCover_data["landCoverFrac"][year_start:int(year_start+month_period/12.)]
    landCover[landCover>10**10] = 0
    landCover[landCover<0] = 0
    landCover = np.divide(landCover,100)

    BA = BA_data["burntArea"][time:time+month_period]
    BA = np.multiply(landCover, BA)
    BA[BA<0.] = 0
    BA[BA>10**10] = 0
    # Assume fractional, add up pft dependency.
    BA = np.sum(BA, axis=1)
    inverse_BA = 1./BA

    
    actual_emission_data= np.multiply(landCover, emis_data["fFirepft"][time:time+month_period])
    # Add up pft dependency.
    actual_emission_data = np.sum(actual_emission_data, axis=1)
    
    FC_data = np.multiply(actual_emission_data, inverse_BA)  
    FC_per_month = np.multiply(FC_data, 
                sec_per_month[:, np.newaxis, np.newaxis])
    fuel_consumption = np.sum(FC_per_month, axis = 0)
    return fuel_consumption

    
def get_global_mean_FC_yearly(year, emis_data, BA_data, grid_data, landCover_data, time_data):
    time = int(year*12)   
    FC_data = get_grid_fuel_consumption(year,12,emis_data, BA_data, landCover_data, time_data) 
    
    #Taking the weighted average of fuel consumption.
    global_mean_FC = np.average(FC_data, weights=grid_data["cell_area"])
    return global_mean_FC 
    
    
def get_global_mean_FC_yearly_rough(year, emis_data, BA_data, grid_data, landCover_data, time_data):
    total_emis = get_global_emissions_yearly(year,emis_data,grid_data,landCover_data,time_data)
    total_BA = get_global_BA_yearly(year,BA_data,grid_data,landCover_data)
    global_mean_FC = total_emis/total_BA
    return global_mean_FC  
    
    
def plot_global_mean_FC_yearly(no_years, emis_data, BA_data, grid_data, landCover_data, time_data):
    years = int(len(time_data["time"])/12)
    x_data = range(years-1)
    x_data_plot = [x+1700 for x in x_data]
    y_data = []
    for x in x_data[-no_years:]:
        print("%.2f" % ((x-x_data[-no_years])/
                float(len(x_data[-no_years:]))))
        y_data.append(get_global_mean_FC_yearly(x,emis_data,BA_data,grid_data,landCover_data, time_data))
    plt.plot(x_data_plot[-no_years:], y_data, color='r', linewidth=2.0,
               label='ORCHIDEE Results')
    plt.ylabel('Fuel Consumption ($kg\, C / m^2 \, burned$)')
    plt.xlabel('Year')
    plt.legend()
    plt.show()
    
    
def plot_FC_map_period(year_start, year_period, emis_data, BA_data, landCover_data, time_data):
    yearly_data = []
    for i in range(year_period):
        yearly_data.append(get_grid_fuel_consumption(year_start,12,
                                emis_data, BA_data, landCover_data, time_data))                             
    map_data = np.mean(yearly_data, axis=0)
    # Replace 0 by nan.
    map_data[map_data==0]=np.nan
    
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
    plt.title("Mean Fuel Consumption, ORCHIDEE")
    plt.show()


#print get_global_mean_FC_yearly_rough(300, emis_ORCHIDEE, BA_ORCHIDEE, grid_ORCHIDEE, landCover_ORCHIDEE,time_data)
    
#plot_FC_map_period(297,16 , emis_ORCHIDEE, BA_ORCHIDEE, landCover_ORCHIDEE, time_data)
#plot_global_mean_FC_yearly(20, emis_ORCHIDEE,BA_ORCHIDEE,grid_ORCHIDEE, landCover_ORCHIDEE, time_data)
#print get_global_mean_FC_yearly(300,emis_ORCHIDEE, BA_ORCHIDEE, grid_ORCHIDEE, landCover_ORCHIDEE, time_data)
