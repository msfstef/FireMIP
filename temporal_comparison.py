import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as clb
from matplotlib import gridspec

import spatial_comparison as spt


def plot_time_series(year, year_period, var, region=0, model='all', 
                    means=False, corr_gfed=False, corr_multimodel=False, 
                    all_regions=False, save=False):
    """
    Plots a time series from the given year and year period (in absolute
    terms, e.g. 1997).
    
    The argument var determines the variable which is plotted, and can take
    values 'emis', 'BA', and 'FC', for emissions, burnt area, and fuel
    consumption respectively.
    
    The argument region takes values from 0 to 12, where 0 is global and
    set as default, and the rest are regions specified in the spatial module
    spatial_comparison.
    
    By default, the time series is plotted for all models, but a model
    can be specified to plot only that model's time series.
    
    The arguments means, corr_gfed, and corr_multimodel determine whether
    tables of means, temporal correlations with GFED, and temporal correlations
    with the multimodel mean will be saved respectively. They are set to
    False by default to not bloat the folder in which they are saved.
    
    The argument all_regions can be set to True to create a plot consisting
    of multiple subplots for all individual regions.
    
    The argument save is used by the generate_figures module and is set
    to False by default.
    """
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
                                            keep_time=True, 
                                            all_regions=all_regions)
            if all_regions:
                yearly_data = np.sum(np.split(grid.sum(axis=2).sum(axis=2),
                                    year_period,axis=1),axis=2)
            else:
                yearly_data = np.sum(np.split(grid.sum(axis=1).sum(axis=1),
                                    year_period),axis=1)
            data_list.append(yearly_data)
            print(model.upper() + ' finished.')
    elif var == 'FC':
        for model in model_list:
            emis_grid = spt.get_regional_var_grid(year,year_period,region,
                                            model,'emis',per_area=False,
                                            keep_time=True, 
                                            all_regions=all_regions)
            if all_regions:
                emis_data = np.sum(np.split(emis_grid.sum(axis=2).sum(axis=2),
                                        year_period,axis=1),axis=2)
            else:
                emis_data = np.sum(np.split(emis_grid.sum(axis=1).sum(axis=1),
                                        year_period),axis=1)
                                        
            BA_grid = spt.get_regional_var_grid(year,year_period,region,
                                            model,'BA',per_area=False,
                                            keep_time=True, 
                                            all_regions=all_regions)
            if all_regions:
                BA_data = np.sum(np.split(BA_grid.sum(axis=2).sum(axis=2),
                                     year_period,axis=1),axis=2)
            else:
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
        units = '(millions of $km^2 \, year^{-1}$)'
        unit_conv=1e12
    
    
    if len(model_list)>1 and not all_regions:
        fig = plt.figure(figsize=(12,8))
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
                
        save_table(data_list,year,year_period,var,region,model_list,
                    means, corr_gfed, corr_multimodel)
        
        plt.xlabel('Year')
        plt.xlim([np.min(years),np.max(years)])  
        plt.ylabel(title +' '+ units)
        plt.ylim([np.min(np.array(data_list)/unit_conv)*0.80, 
                np.max(np.array(data_list)/unit_conv)*1.20])
        plt.title('Yearly Mean of '+title+' for '+str(year)+
            '-'+str(year+year_period-1)+', '+title_end)
        plt.legend(ncol=3)
                    
    elif len(model_list)==1 and not all_regions:
        fig = plt.figure(figsize=(12,8))
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
    
    else:
        fig = plt.figure(figsize=(16,12.5))
        ax_list = []
        gs = gridspec.GridSpec(6,3)
        data_list = np.array(data_list)
        for region in range(13):
            if region==0:
                ax = fig.add_subplot(gs[0:2,:])
            elif 0<region<=4:
                ax = fig.add_subplot(gs[region+1,0])
                
            elif 4<region<=8:
                ax = fig.add_subplot(gs[region-4+1,1])
            elif 8<region:
                ax = fig.add_subplot(gs[region-8+1,2])
            
            data=data_list[:,:,region]
            
            multimodel_list = np.array(data[1:])
            multimodel_mean = np.mean(multimodel_list, axis=0)
            multimodel_std = np.std(multimodel_list, axis=0)

            colour_map=iter(plt.cm.Dark2(
                np.linspace(0,1,len(multimodel_list))))
            for j in range(len(multimodel_list)):
                c = next(colour_map)
                ax.plot(years,np.divide(multimodel_list[j],unit_conv),
                     c=c, label=model_names[j+1], linewidth=1.5)
        
            ax.plot(years, np.divide(data[0],unit_conv), color='k', 
                        linewidth=1.5, linestyle='--', label='GFED')
        
            ax.plot(years, np.divide(multimodel_mean,unit_conv), 'k-.',
                         label='Multimodel Mean')
            ax.fill_between(years, 
                np.divide(multimodel_mean-multimodel_std,unit_conv), 
                np.divide(multimodel_mean+multimodel_std,unit_conv), 
                facecolor='grey', alpha=0.2)
            
           
            ax.set_xlim([np.min(years),np.max(years)])  

            if region==0:
                ax.set_ylim([np.min(np.array(data)/unit_conv)*0.50, 
                   np.max(np.array(data)/unit_conv)*1.50])
                ax.set_ylabel(title +' '+ units)
                ax.set_title('Yearly Mean of '+title+' for '+str(year)+
                 '-'+str(year+year_period-1)+', Global + All Regions')
                ax.legend(ncol=3)
            elif np.mod(region,4)!=0:
                ax.set_xticklabels([])
                ax.set_ylim([np.min(np.array(data)/unit_conv)*0.80, 
                   np.max(np.array(data)/unit_conv)*1.20])
                start, end = ax.get_ylim()
                ax.yaxis.set_ticks(np.linspace(start, end,4))
                ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
                ax.set_title(region_names[region])
            else:
                ax.set_xlabel('Year')
                ax.set_ylim([np.min(np.array(data)/unit_conv)*0.80, 
                   np.max(np.array(data)/unit_conv)*1.20])
                start, end = ax.get_ylim()
                ax.yaxis.set_ticks(np.linspace(start, end,4))
                ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
                ax.set_title(region_names[region])
            ax_list.append(ax)
        plt.tight_layout()
            
    if save:
        return fig
    else:
        plt.show()


def save_table(data, year, year_period, var, region, model_list, 
                means, corr_gfed, corr_multimodel):
    """
    Used by the plot_time_series function to create and save tables
    of means and various correlations. See docstring of said function
    for more details.
    """
    region_names = ['Global','BONA','TENA','EQCSA','SOMA','NOEU',
                    'MEME','EQAF','SOAF','BOAS','CEAS','EQAS','AUST']
    model_names = [label.upper() for label in model_list] 
    model_names = np.array(model_names)
    data = np.array(data)
    
    if var == 'FC':
        title = 'Total Fuel Consumption'
        units = '(kg C per m^2 burned per year)'
        unit_conv = 1.
    elif var == 'emis':
        title = 'Total Carbon Emissions'
        units = '(Pg C per year)'
        unit_conv=1e12
    elif var == 'BA':
        title = 'Total Burnt Area'
        units = '(millions of km^2 per year)'
        unit_conv=1e12
    
    if means:
        means = np.mean(data, axis=1)/unit_conv
        std = np.std(data, axis=1)/unit_conv
        table_data = np.array([model_names, means, std])
        f = open("./figures/temporal_comparison/means_table_"+var+
                +"_"+str(year)+"-"+str(year+year_period-1)+"_"+
                region_names[region]+".csv", "w")
        f.write("Table of Mean Yearly "+title+" for "+str(year)+"-"+
                str(year+year_period-1)+" ~ "+region_names[region]+"\n")
        f.write("Model Name, Mean "+units+", Standard Deviation "+units+"\n")
        np.savetxt(f, table_data.T, delimiter=',', fmt='%s')
        print 'Means table finished!'
    
    if corr_gfed:
        corr_tuple = [stats.pearsonr(data[0],model_data) 
                        for model_data in data[1:]]
        corrval = np.array([corr[0] for corr in corr_tuple])
        pval = np.array([corr[1] for corr in corr_tuple])
        table_data = np.array([model_names[1:], corrval, pval])
        f = open("./figures/temporal_comparison/"+
                "GFED_correlations_table_"+var+
                +"_"+str(year)+"-"+str(year+year_period-1)+"_"+
                region_names[region]+".csv", "w")
        f.write("Table of Temporal Pearson Correlations of "+title+
                " with GFED for "+str(year)+"-"+str(year+year_period-1)+" ~ "
                +region_names[region]+"\n")
        f.write("Model Name,Pearson's r,p-value\n")
        np.savetxt(f, table_data.T, delimiter=',', fmt='%s')
        print 'GFED Correlations table finished!'
    
    if corr_multimodel:
        multimodel = np.mean(data[1:], axis=0)
        corr_tuple = [stats.pearsonr(multimodel,model_data) 
                        for model_data in data]
        corrval = np.array([corr[0] for corr in corr_tuple])
        pval = np.array([corr[1] for corr in corr_tuple])
        table_data = np.array([model_names, corrval, pval])
        f = open("./figures/temporal_comparison/"+
                "multimodel_correlations_table_"+var+
                +"_"+str(year)+"-"+str(year+year_period-1)+"_"+
                region_names[region]+".csv", "w")
        f.write("Table of Temporal Pearson Correlations of "+title+
                " with the Multimodel Mean for "+str(year)+"-"+
                str(year+year_period-1)+" ~ "+region_names[region]+"\n")
        f.write("Model Name,Pearson's r,p-value\n")
        np.savetxt(f, table_data.T, delimiter=',', fmt='%s')
        print 'Multimodel Correlations table finished!'


#plot_time_series(1997,16,'emis',all_regions=True)
            
