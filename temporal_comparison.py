import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.colors as clb

import spatial_comparison as spt


def plot_time_series(year, year_period, var, region=0, model='all', save=False):
    model_list = ['gfed', 'jsbach', 'clm', 'ctem', 
                'blaze', 'orchidee', 'inferno', 'spitfire']       
    if model != 'all':
        model_list = [model]
    model_names = [label.upper() for label in model_list] 
    region_names = ['Global','BONA','TENA','EQCSA','SOMA','NOEU',
                    'MEME','EQAF','SOAF','BOAS','CEAS','EQAS','AUST']
    years = np.arange(year,year+year_period,1)
    data_list=[]
    
    if var != 'FC':
        for model in model_list:
            grid = spt.get_regional_var_grid(year,year_period,region,
                                            model,var,per_area=False,
                                            keep_time=True)
            yearly_data = np.sum(np.split(grid.sum(axis=1).sum(axis=1),
                                    year_period),axis=1)
            data_list.append(yearly_data)
            print(model.upper() + ' finished.')
    elif var == 'FC':
        for model in model_list:
            emis_grid = spt.get_regional_var_grid(year,year_period,region,
                                            model,'emis',per_area=False,
                                            keep_time=True)
            emis_data = np.sum(np.split(emis_grid.sum(axis=1).sum(axis=1),
                                    year_period),axis=1)
                                    
            BA_grid = spt.get_regional_var_grid(year,year_period,region,
                                            model,'BA',per_area=False,
                                            keep_time=True)
            BA_data = np.sum(np.split(BA_grid.sum(axis=1).sum(axis=1),
                                    year_period),axis=1)
            
            yearly_data = np.divide(emis_data, BA_data)
            data_list.append(yearly_data)
            print(model.upper() + ' finished.')    
    
    if var == 'FC':
        title = 'Total Fuel Consumption'
        units = '($kg\, C\, m^{-2}\, burned \, year^{-1}$)'
        unit_conv = 1.
    elif var == 'emis':
        title = 'Total Carbon Emissions'
        units = '($Pg\, C\, year^{-1}$)'
        unit_conv=1e12
    elif var == 'BA':
        title = 'Total Burnt Area'
        units = '(millions of $km^2$)'
        unit_conv=1e12
    
    # Script to save raw yearly data to CSV file.
    #data = [years]+data_list
    #f = open("yearly_data.csv", "w")
    #f.write("Year,JSBACH,CLM,CTEM,BLAZE,ORCHIDEE,INFERNO,SPITFIRE\n")
    #np.savetxt(f, np.array(data).T, delimiter=',')
    
    
    fig = plt.figure(figsize=(12,8))
    
    if len(model_list)>1:
        multimodel_list = np.array(data_list[1:])
        multimodel_mean = np.mean(multimodel_list, axis=0)
        multimodel_std = np.std(multimodel_list, axis=0)
        
        title_end = region_names[region]
        colour_map=iter(plt.cm.Dark2(
                np.linspace(0,1,len(multimodel_list))))
        for i in range(len(multimodel_list)):
            c = next(colour_map)
            plt.plot(years,np.divide(multimodel_list[i],unit_conv),
                 c=c, label=model_names[i+1], linewidth=1.5)
        
        plt.plot(years, np.divide(data_list[0],unit_conv), color='k', 
                        linewidth=1.5, linestyle='--', label='GFED')
        
        plt.plot(years, np.divide(multimodel_mean,unit_conv), 'k-.',
                         label='Multimodel Mean')
        plt.fill_between(years, 
                np.divide(multimodel_mean-multimodel_std,unit_conv), 
                np.divide(multimodel_mean+multimodel_std,unit_conv), 
                facecolor='grey', alpha=0.2)
    else:
        title_end = model_names[0]
        plt.plot(years, np.divide(data_list[0],unit_conv),'r',
                    label=model_names[0], linewidth=1.5)
                    
    
    
    plt.xlabel('Year')
    plt.xlim([np.min(years),np.max(years)])  
    plt.ylabel(title +' '+ units)
    plt.ylim([np.min(np.array(data_list)/unit_conv)*0.80, 
            np.max(np.array(data_list)/unit_conv)*1.20])
    plt.title('Yearly Mean of '+title+' for '+str(year)+
            '-'+str(year+year_period-1)+', '+title_end)
    plt.legend(ncol=3)
    if save:
        return fig
    else:
        plt.show()

#plot_time_series(1997,16,'BA', 11)
