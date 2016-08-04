import numpy as np
import scipy.stats as stats
from scipy import spatial
import matplotlib.pyplot as plt


import spatial_comparison as spt


#FIELD OBSERVATIONS DATA
observ_latlon = [[18.35,-95.05],[36,-79.1],[-33.93,115.46],
[-2.29,-60.09],[34.1,-117.47],[40,-2],[-7.98,-38.32],
[-15.84,-47.95],[-5.35,-49.15],[35.21,-83.48],[-4.3,-49.03],
[-2.34,-60.09],[-12.4,132.5],[-12.3,133],[60.49,-150.98],
[-5.3,-49.15],[-4.3,-49.03],[-2.61,-60.17],[-25.15,31.14],
[-9.2,-60.5],[-12.35,30.21],[-16.6,27.15],[8.56,-67.25],
[62.69,-141.77],[-15.85,-60.52],[36.6,-118.81],[-15.84,-47.95],
[60.45,89.25],[2.54,-61.28],[61.61,-149.04],[19.3,-105.3],
[-2.5,-48.12],[-9.17,-63.18],[19.3,-105.3],[-2.61,-60.17],
[-10.16,-60.81],[-9.11,-63.16],[-14.52,24.49],[63.38,-158.25],
[61.67,-117.65],[-9.52,-56.06],[-0.52,117.01],[17.59,81.55],
[64.4,-145.74],[55.85,-107.67],[-15,23],[17.65,81.75],[-43.22,146.54],
[54.93,-114.17],[42.4,-124.1],[-2.75,114.32],[38.9,-120.67],
[53.92,-105.7],[-12.38,133.55],[46.73,-117.18],[65,-146],
[-35.77,148.03],[65.03,-147.85],[-9.97,-56.34],[-12.43,131.49],
[12.22,-2.7],[19.5,-99.5],[-2.52,113.79],[-33.68,116.25],
[68.58,-149.72],[-34.2,116.34],[-12.53,-54.88],[-33.91,116.16],
[-12.38,133.55],[34.73,-120.57],[40,-98],[-7.9,-72.44],[34.01,-80.72]]

# Units are in tonnes of C per hectare burned.
observ_FC = [380,2.3,43,77.4,45,1.1,64.2,7.2,51.6,92.7,102.9,
111.4,4.5,5.1,32.9,44,55,106.6,3.5,4,4.2,4.5,5.5,48.8,165.5,
212,8.2,16.9,24,51.4,91,43.1,20.7,23.6,82.2,109.9,190.5,2.2,
36.7,42.7,191.5,120,4,35,42,2.9,7.7,299,42.5,38,109,111.4,26.6,
1.4,11,40.2,47,56.7,72.8,2.4,4.2,17.4,331.7,28,40,58,154.6,53,
4.8,116,2.1,226.2,6.3]
# Convert to kg C per m^2 burned.
observ_FC = np.array(observ_FC)*0.1

observ_year = [1980,1983,1983,1984,1986,1989,1989,1990,1990,1990,
1990,1990,1991,1991,1991,1991,1991,1991,1992,1992,1992,1992,1992,
1992,1992,1992,1993,1993,1993,1993,1993,1994,1995,1995,1995,1995,
1995,1996,1997,1997,1997,1998,1999,1999,1999,2000,2000,2000,2001,
2002,2002,2002,2003,2004,2004,2004,2004,2004,2004,2005,2005,2006,
2006,2007,2007,2007,2007,2008,2009,2009,2010,2010,2011]


def find_indices(model):
    lons_list, lats_list = spt.get_lons_lats(model)
    lats, lons = np.meshgrid(lats_list,lons_list)
    latlon_data = zip(lats.ravel(), lons.ravel())
    model_grid = spatial.KDTree(latlon_data)
    query = model_grid.query(observ_latlon)
    
    indices = query[1]
    model_shape = (len(lats_list), len(lons_list))
    lat_ind=np.mod(indices,len(lats_list))
    lon_ind= np.divide(indices,len(lats_list))
    model_indices = [(lat,lon) for lat,lon 
                        in zip(lat_ind,lon_ind)]
    return model_indices, model_shape

def get_observ_grid(model):
    indices, shape = find_indices(model)
    grid = np.zeros(shape,dtype=object)
    for i in range(len(indices)):
        index = indices[i]
        value = observ_FC[i]
        if grid[index] == 0:
            grid[index] = [value]
        else:
            grid[index].append(value)
    for index,FC_list in np.ndenumerate(grid):
        grid[index] = np.mean(FC_list)
    grid = grid.astype('float64')
    return grid
    

def compare_points(model, region=0):
    observ_grid = get_observ_grid(model)
    observ_grid = spt.get_regional_var_grid(1970,42,
                                    region,model,'FC',
                                    grid = observ_grid)
    model_grid = spt.get_regional_var_grid(1970,42,
                                    region,model,'FC')
                                    
    diff_grid = model_grid-observ_grid
    # Avoid division by 0.
    observ_grid[observ_grid==0]=np.inf
    diff_grid = np.divide(diff_grid,observ_grid)
    diff_grid = np.multiply(diff_grid, 100)
    return diff_grid


def calc_mean_dev_points(model, region=0):
    diff_grid = compare_points(model, region)
    diff_vals = diff_grid[diff_grid!=0]
    mean_dev = np.mean(diff_vals)
    no_obs = len(diff_vals)
    stderr_dev = np.std(diff_vals)/np.sqrt(no_obs)
    return (mean_dev, stderr_dev)


def calc_mean_dev_total(model, region=0):
    observ_grid = get_observ_grid(model)
    observ_grid = spt.get_regional_var_grid(1980,32,
                                    region,model,'FC',
                                    grid = observ_grid)
    model_grid = spt.get_regional_var_grid(1980,32,
                                    region,model,'FC')
    
    observ_vals = observ_grid[observ_grid>0]
    model_vals = model_grid[model_grid>0]
    no_obs = len(observ_vals)
    no_model = len(model_vals)
    observ_mean = np.mean(observ_vals)
    observ_stderr = np.std(observ_vals)/np.sqrt(no_obs)
    model_mean = np.mean(model_vals)
    model_stderr = np.std(model_vals)/np.sqrt(no_model)
    
    mean_dev = ((model_mean-observ_mean)/observ_mean)*100
    stderr_dev = ((((model_stderr+observ_stderr)/
              (model_mean+observ_mean)) +
              (observ_stderr/observ_mean))*mean_dev)
    return (mean_dev, stderr_dev)


def plot_bar_chart(model='all'):
    model_list = ['gfed', 'jsbach', 'clm', 'ctem', 
                'blaze', 'orchidee', 'inferno']
    colour_list = ['k','b','g','r','c','m','y']
    region_names = ['Global','BONA','TENA','EQCSA','SOMA','NOEU',
                   'MEME','EQAF','SOAF','BOAS','CEAS','EQAS','AUST']
    #model_list=['jsbach','ctem','clm']
    fig, ax = plt.subplots()
    error_config = {'ecolor': '0.3'}
    no_of_bars = len(region_names)
    ind = np.arange(no_of_bars)*2
    bar_width = 0.20
        
    for i in range(len(model_list)):
        model_name = model_list[i]
        colour = colour_list[i]
        means = []
        stderrs = []
        for j in range(13):
            dev = calc_mean_dev_points(model_name,region=j)
            means.append(dev[0])
            stderrs.append(dev[1])
        
        ax.bar(ind+i*bar_width,means,bar_width,
             color=colour,yerr=stderrs, 
             error_kw=error_config,label=model_name)
        
    ax.set_ylabel('Relative Deviation (%)')
    ax.set_xticks(ind+bar_width)
    ax.set_xticklabels(region_names)
    ax.legend()
    
    plt.show()
    
     
#plot_bar_chart()
#print calc_mean_dev_points('clm',12)
#print calc_mean_dev_total('clm',12)

