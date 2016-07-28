"""
This module is used for temporal analysis of the
models. For a spatial analysis, use the
spatial_comparison module.
"""

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt

import jsbach_analysis as jsbach
import clm_analysis as clm
import ctem_analysis as ctem
import blaze_analysis as blaze
import orchidee_analysis as orchidee
import inferno_analysis as inferno
import spitfire_analysis as spitfire
import globfirm_analysis as globfirm

import gfed_analysis as gfed

def plot_present_emissions(table=True, pearson=False):
    """
    Plots a time series of global carbon emissions for all
    models from 1997 to 2012, along with a multimodel mean
    and a shaded region that shows the 1 sigma range from
    the mean. The GFED 4.1s results are also plotted for
    reference.
    
    If argument pearson is set to True, tables with the
    models' temporal correlations with GFED and the multimodel
    mean will be printed to the console.
    """
    jsbach_data,clm_data,ctem_data = [],[],[]
    blaze_data,orchidee_data,inferno_data = [],[],[]
    gfed_data=[]
    years = range(1997,2013,1)
    for i in range(16):
        # Progress Bar
        print("%.2f" % (float(i)/15.))
        jsbach_data.append(jsbach.get_global_emissions_yearly(297+i,
            jsbach.emis_JSBACH, jsbach.grid_JSBACH))
        clm_data.append(clm.get_global_emissions_yearly(297+i,
            clm.emis_CLM, clm.grid_CLM, clm.time_data))
        ctem_data.append(ctem.get_global_emissions_yearly(136+i,
            ctem.emis_CTEM, ctem.grid_CTEM, ctem.landCover_CTEM))
        blaze_data.append(blaze.get_global_emissions_yearly(297+i,
            blaze.emis_BLAZE, blaze.grid_BLAZE, blaze.time_data))
        orchidee_data.append(orchidee.get_global_emissions_yearly(297+i,
            orchidee.emis_ORCHIDEE, orchidee.grid_ORCHIDEE,
            orchidee.landCover_ORCHIDEE, orchidee.time_data))
        inferno_data.append(inferno.get_global_emissions_yearly(297+i,
            inferno.emis_INFERNO, inferno.grid_INFERNO, 
            inferno.landmask_INFERNO, inferno.landCover_INFERNO))
        
        gfed_data.append(gfed.get_global_emissions_yearly(i,
                        gfed.data_GFED,gfed.grid_GFED))
    
    multimodel_list = np.array([jsbach_data,clm_data,ctem_data,blaze_data,orchidee_data])
    multimodel_mean = np.mean(multimodel_list, axis=0)
    multimodel_std = np.std(multimodel_list, axis=0)
    
    if table:
        print "Table of Mean Emissions"
        print "Model Name || (mean,stdev in Pg C)"
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print('JSBACH: '+str(np.mean(jsbach_data)/1e12)+
                ','+str(np.std(jsbach_data)/1e12))
        print('CLM: '+str(np.mean(clm_data)/1e12)+
                ','+str(np.std(clm_data)/1e12))
        print('CTEM: '+str(np.mean(ctem_data)/1e12)+
                ','+str(np.std(ctem_data)/1e12))
        print('LPJ_GUESS_BLAZE: '+str(np.mean(blaze_data)/1e12)+
                ','+str(np.std(blaze_data)/1e12))
        print('ORCHIDEE: '+str(np.mean(orchidee_data)/1e12)+
                ','+str(np.std(orchidee_data)/1e12))
        print('INFERNO: '+str(np.mean(inferno_data)/1e12)+
                ','+str(np.std(inferno_data)/1e12))
        print('Multimodel: '+str(np.mean(multimodel_mean)/1e12)+
                ','+str(np.std(multimodel_mean)/1e12))     

    
    if pearson:
        jsbach_r = stats.pearsonr(gfed_data, jsbach_data)
        clm_r = stats.pearsonr(gfed_data, clm_data)
        ctem_r = stats.pearsonr(gfed_data, ctem_data)
        blaze_r = stats.pearsonr(gfed_data, blaze_data)
        orchidee_r = stats.pearsonr(gfed_data, orchidee_data)
        inferno_r = stats.pearsonr(gfed_data, inferno_data)
        
        jsbach_r2 = stats.pearsonr(multimodel_mean, jsbach_data)
        clm_r2 = stats.pearsonr(multimodel_mean, clm_data)
        ctem_r2 = stats.pearsonr(multimodel_mean, ctem_data)
        blaze_r2 = stats.pearsonr(multimodel_mean, blaze_data)
        orchidee_r2 = stats.pearsonr(multimodel_mean, orchidee_data)
        inferno_r2 = stats.pearsonr(multimodel_mean, inferno_data)
        
        print "Table of Emissions correlations with GFED"
        print "Model Name || (Pearson's r, p-value)"
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print 'JSBACH:', jsbach_r
        print 'CLM:', clm_r
        print 'CTEM:', ctem_r
        print 'LPJ_GUESS_BLAZE:', blaze_r
        print 'ORCHIDEE:', orchidee_r
        print 'INFERNO:', inferno_r
        print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        print "Table of Emissions correlations with Multimodel Mean"
        print "Model Name || (Pearson's r, p-value)"
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print 'JSBACH:', jsbach_r2
        print 'CLM:', clm_r2
        print 'CTEM:', ctem_r2
        print 'LPJ_GUESS_BLAZE:', blaze_r2
        print 'ORCHIDEE:', orchidee_r2
        print 'INFERNO:', inferno_r2
    
    plt.plot(years, jsbach_data, 'r', label='JSBACH')
    plt.plot(years, clm_data, 'b', label='CLM')
    plt.plot(years, ctem_data, 'g', label='CTEM')
    plt.plot(years, blaze_data, 'm', label='LPJ-GUESS-BLAZE')
    plt.plot(years, orchidee_data, 'y', label='ORCHIDEE')
    plt.plot(years, inferno_data, 'c', label='INFERNO')
    
    plt.plot(years, gfed_data, color='k', linewidth=1.5, linestyle='--', label='GFED')
    
    plt.plot(years, multimodel_mean, 'k-.', label='Multimodel Mean')
    plt.fill_between(years, multimodel_mean-multimodel_std, 
                    multimodel_mean+multimodel_std, facecolor='grey', alpha=0.2)
    
    plt.xlabel('Year')
    plt.ylabel('Total Carbon Emissions ($Pg\, C/year$)')
    
    plt.legend()
    plt.show()


def plot_present_burnt_area(table=True, pearson=False):
    """
    Plots a time series of global total burnt area for all
    models from 1997 to 2012, along with a multimodel mean
    and a shaded region that shows the 1 sigma range from
    the mean. The GFED 4.1s results are also plotted for
    reference.
    
    If argument pearson is set to True, tables with the
    models' temporal correlations with GFED and the multimodel
    mean will be printed to the console.
    """
    jsbach_data,clm_data,ctem_data = [],[],[]
    blaze_data,orchidee_data,inferno_data = [],[],[]
    gfed_data=[]
    years = range(1997,2013,1)
    for i in range(16):
        # Progress Bar
        print("%.2f" % (float(i)/15.))
        jsbach_data.append(jsbach.get_global_BA_yearly(297+i,
            jsbach.BA_JSBACH, jsbach.grid_JSBACH))
        clm_data.append(clm.get_global_BA_yearly(297+i,
            clm.BA_CLM, clm.grid_CLM, clm.time_data))
        ctem_data.append(ctem.get_global_BA_yearly(136+i,
            ctem.BA_CTEM, ctem.grid_CTEM, ctem.landCover_CTEM))
        blaze_data.append(blaze.get_global_BA_yearly(297+i,
            blaze.BA_BLAZE, blaze.grid_BLAZE))
        orchidee_data.append(orchidee.get_global_BA_yearly(297+i,
            orchidee.BA_ORCHIDEE, orchidee.grid_ORCHIDEE,
            orchidee.landCover_ORCHIDEE))
        inferno_data.append(inferno.get_global_BA_yearly(297+i,
            inferno.BA_INFERNO, inferno.grid_INFERNO, 
            inferno.landmask_INFERNO, inferno.landCover_INFERNO))
        
        gfed_data.append(gfed.get_global_BA_yearly(i,
                        gfed.data_GFED,gfed.grid_GFED))
    
    multimodel_list = np.array([jsbach_data,clm_data,ctem_data,blaze_data,orchidee_data])
    multimodel_mean = np.mean(multimodel_list, axis=0)
    multimodel_std = np.std(multimodel_list, axis=0)
    
    if table:
        print "Table of Mean Burnt Area"
        print "Model Name || (mean,stdev in millions of km^2)"
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print('JSBACH: '+str(np.mean(jsbach_data)/1e12)+
                ','+str(np.std(jsbach_data)/1e12))
        print('CLM: '+str(np.mean(clm_data)/1e12)+
                ','+str(np.std(clm_data)/1e12))
        print('CTEM: '+str(np.mean(ctem_data)/1e12)+
                ','+str(np.std(ctem_data)/1e12))
        print('LPJ_GUESS_BLAZE: '+str(np.mean(blaze_data)/1e12)+
                ','+str(np.std(blaze_data)/1e12))
        print('ORCHIDEE: '+str(np.mean(orchidee_data)/1e12)+
                ','+str(np.std(orchidee_data)/1e12))
        print('INFERNO: '+str(np.mean(inferno_data)/1e12)+
                ','+str(np.std(inferno_data)/1e12))
        print('Multimodel: '+str(np.mean(multimodel_mean)/1e12)+
                ','+str(np.std(multimodel_mean)/1e12))   
    
    if pearson:
        jsbach_r = stats.pearsonr(gfed_data, jsbach_data)
        clm_r = stats.pearsonr(gfed_data, clm_data)
        ctem_r = stats.pearsonr(gfed_data, ctem_data)
        blaze_r = stats.pearsonr(gfed_data, blaze_data)
        orchidee_r = stats.pearsonr(gfed_data, orchidee_data)
        inferno_r = stats.pearsonr(gfed_data, inferno_data)
        
        jsbach_r2 = stats.pearsonr(multimodel_mean, jsbach_data)
        clm_r2 = stats.pearsonr(multimodel_mean, clm_data)
        ctem_r2 = stats.pearsonr(multimodel_mean, ctem_data)
        blaze_r2 = stats.pearsonr(multimodel_mean, blaze_data)
        orchidee_r2 = stats.pearsonr(multimodel_mean, orchidee_data)
        inferno_r2 = stats.pearsonr(multimodel_mean, inferno_data)
        
        print "Table of BA correlations with GFED"
        print "Model Name || (Pearson's r, p-value)"
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print 'JSBACH:', jsbach_r
        print 'CLM:', clm_r
        print 'CTEM:', ctem_r
        print 'LPJ_GUESS_BLAZE:', blaze_r
        print 'ORCHIDEE:', orchidee_r
        print 'INFERNO:', inferno_r
        print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        print "Table of BA correlations with Multimodel Mean"
        print "Model Name || (Pearson's r, p-value)"
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print 'JSBACH:', jsbach_r2
        print 'CLM:', clm_r2
        print 'CTEM:', ctem_r2
        print 'LPJ_GUESS_BLAZE:', blaze_r2
        print 'ORCHIDEE:', orchidee_r2
        print 'INFERNO:', inferno_r2
    
    plt.plot(years, np.divide(jsbach_data, 1e12), 'r', label='JSBACH')
    plt.plot(years, np.divide(clm_data, 1e12), 'b', label='CLM')
    plt.plot(years, np.divide(ctem_data, 1e12), 'g', label='CTEM')
    plt.plot(years, np.divide(blaze_data, 1e12), 'm', label='LPJ-GUESS-BLAZE')
    plt.plot(years, np.divide(orchidee_data, 1e12), 'y', label='ORCHIDEE')
    plt.plot(years, np.divide(inferno_data, 1e12), 'c', label='INFERNO')
    
    plt.plot(years, np.divide(gfed_data,1e12), color='k', 
                    linewidth=1.5, linestyle='--', label='GFED')
    
    plt.plot(years, np.divide(multimodel_mean,1e12), 'k-.', label='Multimodel Mean')
    plt.fill_between(years, np.divide(multimodel_mean-multimodel_std,1e12), 
                    np.divide(multimodel_mean+multimodel_std,1e12), 
                    facecolor='grey', alpha=0.2)
    
    plt.xlabel('Year')
    plt.ylabel('Total Burnt Area (millions of $km^2$)')
    
    plt.legend()
    plt.show()
    
    
def plot_present_fuel_consumption(table=True, pearson = False):
    """
    Plots a time series of mean global fuel consumption (FC),
    calculated by dividing total global emissions by total
    global burnt area over the desired time period, for all
    models from 1997 to 2012, along with a multimodel mean
    and a shaded region that shows the 1 sigma range from
    the mean. The GFED 4.1s results are also plotted for
    reference.
    
    If argument pearson is set to True, tables with the
    models' temporal correlations with GFED and the multimodel
    mean will be printed to the console.
    """
    jsbach_data,clm_data,ctem_data = [],[],[]
    blaze_data,orchidee_data, inferno_data = [],[],[]
    gfed_data=[]
    years = range(1997,2013,1)
    for i in range(16):
        # Progress Bar
        print("%.2f" % (float(i)/15.))
        jsbach_data.append(jsbach.get_global_mean_FC_yearly_rough(297+i,
            jsbach.emis_JSBACH, jsbach.BA_JSBACH, jsbach.grid_JSBACH))
        clm_data.append(clm.get_global_mean_FC_yearly_rough(297+i,
            clm.emis_CLM, clm.BA_CLM, clm.grid_CLM, clm.time_data))
        ctem_data.append(ctem.get_global_mean_FC_yearly_rough(136+i,
            ctem.emis_CTEM, ctem.BA_CTEM, ctem.grid_CTEM, ctem.landCover_CTEM))
        blaze_data.append(blaze.get_global_mean_FC_yearly_rough(297+i,
            blaze.emis_BLAZE, blaze.BA_BLAZE, blaze.grid_BLAZE, blaze.time_data))
        orchidee_data.append(orchidee.get_global_mean_FC_yearly_rough(297+i,
            orchidee.emis_ORCHIDEE, orchidee.BA_ORCHIDEE, orchidee.grid_ORCHIDEE,
            orchidee.landCover_ORCHIDEE, orchidee.time_data))
        inferno_data.append(inferno.get_global_mean_FC_yearly_rough(297+i,
            inferno.emis_INFERNO, inferno.BA_INFERNO, inferno.grid_INFERNO,
            inferno.landmask_INFERNO, inferno.landCover_INFERNO))
            
        gfed_data.append(gfed.get_global_mean_FC_yearly_rough(i,
                        gfed.data_GFED,gfed.grid_GFED))


    multimodel_list = np.array([jsbach_data,clm_data,ctem_data,
                            blaze_data,orchidee_data, inferno_data])
    multimodel_mean = np.mean(multimodel_list, axis=0)
    multimodel_std = np.std(multimodel_list, axis=0)
    
    if table:
        print "Table of Temporal Mean of Mean Fuel Consumption"
        print "Model Name || (mean,stdev in kg C per m^2 burned)"
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print('JSBACH: '+str(np.mean(jsbach_data))+
                ','+str(np.std(jsbach_data)))
        print('CLM: '+str(np.mean(clm_data))+
                ','+str(np.std(clm_data)))
        print('CTEM: '+str(np.mean(ctem_data))+
                ','+str(np.std(ctem_data)))
        print('LPJ_GUESS_BLAZE: '+str(np.mean(blaze_data))+
                ','+str(np.std(blaze_data)))
        print('ORCHIDEE: '+str(np.mean(orchidee_data))+
                ','+str(np.std(orchidee_data)))
        print('INFERNO: '+str(np.mean(inferno_data))+
                ','+str(np.std(inferno_data)))
        print('Multimodel: '+str(np.mean(multimodel_mean))+
                ','+str(np.std(multimodel_mean))) 
    
    if pearson:
        jsbach_r = stats.pearsonr(gfed_data, jsbach_data)
        clm_r = stats.pearsonr(gfed_data, clm_data)
        ctem_r = stats.pearsonr(gfed_data, ctem_data)
        blaze_r = stats.pearsonr(gfed_data, blaze_data)
        orchidee_r = stats.pearsonr(gfed_data, orchidee_data)
        inferno_r = stats.pearsonr(gfed_data, inferno_data)
        
        jsbach_r2 = stats.pearsonr(multimodel_mean, jsbach_data)
        clm_r2 = stats.pearsonr(multimodel_mean, clm_data)
        ctem_r2 = stats.pearsonr(multimodel_mean, ctem_data)
        blaze_r2 = stats.pearsonr(multimodel_mean, blaze_data)
        orchidee_r2 = stats.pearsonr(multimodel_mean, orchidee_data)
        inferno_r2 = stats.pearsonr(multimodel_mean, inferno_data)
        
        print "Table of FC correlations with GFED"
        print "Model Name || (Pearson's r, p-value)"
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print 'JSBACH:', jsbach_r
        print 'CLM:', clm_r
        print 'CTEM:', ctem_r
        print 'LPJ_GUESS_BLAZE:', blaze_r
        print 'ORCHIDEE:', orchidee_r
        print 'INFERNO:', inferno_r
        print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        print "Table of FC correlations with Multimodel Mean"
        print "Model Name || (Pearson's r, p-value)"
        print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        print 'JSBACH:', jsbach_r2
        print 'CLM:', clm_r2
        print 'CTEM:', ctem_r2
        print 'LPJ_GUESS_BLAZE:', blaze_r2
        print 'ORCHIDEE:', orchidee_r2
        print 'INFERNO:', inferno_r2
    
    plt.plot(years, jsbach_data, 'r', label='JSBACH')
    plt.plot(years, clm_data, 'b', label='CLM')
    plt.plot(years, ctem_data, 'g', label='CTEM')
    plt.plot(years, blaze_data, 'm', label='LPJ-GUESS-BLAZE')
    plt.plot(years, orchidee_data, 'y', label='ORCHIDEE')
    plt.plot(years, inferno_data, 'c', label='INFERNO')
    
    plt.plot(years, gfed_data, color='k', linewidth=1.5, linestyle='--', label='GFED')
    
    plt.plot(years, multimodel_mean, 'k-.', label='Multimodel Mean')
    plt.fill_between(years, multimodel_mean-multimodel_std, 
                    multimodel_mean+multimodel_std, facecolor='grey', alpha=0.2)
    
    plt.xlabel('Year')
    plt.ylabel('Mean Fuel Consumption ($kg\, C \, m^{-2} \, burned \, year^{-1}$)')
    
    plt.legend()
    plt.show()


#plot_present_burnt_area(True)
#plot_present_fuel_consumption_rough(True)
#plot_present_fuel_consumption()
#plot_present_emissions(True)


#
# Past Transient Analysis
#    

def movingaverage(data, period=20):
    """
    Simple moving average function, takes arguments
    data, self explanatory, and period in years over
    which the averaging is done. The default is 20
    years, as it is the period of the climatological
    cycle used in past transient simulations.
    """
    weights = np.repeat(1.0, period)/period
    ma = np.convolve(data, weights, 'valid')
    return ma


def plot_past_emissions():
    """
    Plots a time series of global carbon emissions for all
    models from 1700 to 1996, along with a multimodel mean
    and a shaded region that shows the 1 sigma range from
    the mean.
    """
    jsbach_data,clm_data,ctem_data = [],[],[]
    blaze_data,orchidee_data, inferno_data = [],[],[]
    years = range(1700,1997,1)
    for i in range(len(years)):
        # Progress Bar
        print("%.2f" % (float(i)/float(len(years))))
        jsbach_data.append(jsbach.get_global_emissions_yearly(i,
            jsbach.emis_JSBACH, jsbach.grid_JSBACH))
        clm_data.append(clm.get_global_emissions_yearly(i,
            clm.emis_CLM, clm.grid_CLM, clm.time_data))    
        blaze_data.append(blaze.get_global_emissions_yearly(i,
            blaze.emis_BLAZE, blaze.grid_BLAZE, blaze.time_data))
        orchidee_data.append(orchidee.get_global_emissions_yearly(i,
            orchidee.emis_ORCHIDEE, orchidee.grid_ORCHIDEE,
            orchidee.landCover_ORCHIDEE, orchidee.time_data))
        inferno_data.append(inferno.get_global_emissions_yearly(i,
            inferno.emis_INFERNO, inferno.grid_INFERNO, 
            inferno.landmask_INFERNO, inferno.landCover_INFERNO))
            
    for i in range(len(years[161:])):
        # Progress Bar
        print("%.2f" % (float(i)/float(len(years[161:]))))
        ctem_data.append(ctem.get_global_emissions_yearly(i,
                ctem.emis_CTEM, ctem.grid_CTEM, ctem.landCover_CTEM))

    
    jsbach_data_ma = movingaverage(jsbach_data)
    clm_data_ma = movingaverage(clm_data)
    empty = np.empty(161)*np.nan
    ctem_data_ma = np.concatenate((empty,movingaverage(ctem_data)))
    blaze_data_ma = movingaverage(blaze_data)    
    orchidee_data_ma = movingaverage(orchidee_data)
    inferno_data_ma = movingaverage(inferno_data)    
    
    multimodel_list = np.array([jsbach_data_ma,clm_data_ma,ctem_data_ma,
                            blaze_data_ma,orchidee_data_ma])                       
    multimodel_mean = np.nanmean(multimodel_list, axis=0)
    multimodel_std = np.nanstd(multimodel_list, axis=0)
    
    plt.plot(years[19:], jsbach_data_ma, 'r', label='JSBACH')
    plt.plot(years[19:], clm_data_ma, 'b', label='CLM')
    plt.plot(years[19:], ctem_data_ma, 'g', label='CTEM')
    plt.plot(years[19:], blaze_data_ma, 'm', label='LPJ-GUESS-BLAZE')
    plt.plot(years[19:], orchidee_data_ma, 'y', label='ORCHIDEE')
    plt.plot(years[19:], inferno_data_ma, 'c', label='INFERNO')
    
    plt.plot(years[19:], multimodel_mean, 'k--', label='Multimodel Mean')
    plt.fill_between(years[19:], multimodel_mean-multimodel_std, 
                    multimodel_mean+multimodel_std, facecolor='grey', alpha=0.2)
    
    plt.xlabel('Year')
    plt.ylabel('Total Carbon Emissions ($Pg\, C/year$)')
    
    plt.legend()
    plt.show()


def plot_past_burnt_area():
    """
    Plots a time series of global total burnt area for all
    models from 1700 to 1996, along with a multimodel mean
    and a shaded region that shows the 1 sigma range from
    the mean.
    """
    jsbach_data,clm_data,ctem_data = [],[],[]
    blaze_data,orchidee_data, inferno_data = [],[],[]
    years = range(1700,1997,1)
    for i in range(len(years)):
        # Progress Bar
        print("%.2f" % (float(i)/float(len(years))))
        jsbach_data.append(jsbach.get_global_BA_yearly(i,
            jsbach.BA_JSBACH, jsbach.grid_JSBACH))
        clm_data.append(clm.get_global_BA_yearly(i,
            clm.BA_CLM, clm.grid_CLM, clm.time_data))
        blaze_data.append(blaze.get_global_BA_yearly(i,
            blaze.BA_BLAZE, blaze.grid_BLAZE))
        orchidee_data.append(orchidee.get_global_BA_yearly(i,
            orchidee.BA_ORCHIDEE, orchidee.grid_ORCHIDEE,
            orchidee.landCover_ORCHIDEE))
        inferno_data.append(inferno.get_global_BA_yearly(i,
            inferno.BA_INFERNO, inferno.grid_INFERNO, 
            inferno.landmask_INFERNO, inferno.landCover_INFERNO))
            
    for i in range(len(years[161:])):
        # Progress Bar
        print("%.2f" % (float(i)/float(len(years[161:]))))
        ctem_data.append(ctem.get_global_BA_yearly(i,
            ctem.BA_CTEM, ctem.grid_CTEM, ctem.landCover_CTEM))

    
    jsbach_data_ma = movingaverage(jsbach_data)
    clm_data_ma = movingaverage(clm_data)
    empty = np.empty(161)*np.nan
    ctem_data_ma = np.concatenate((empty,movingaverage(ctem_data)))
    blaze_data_ma = movingaverage(blaze_data)    
    orchidee_data_ma = movingaverage(orchidee_data)
    inferno_data_ma = movingaverage(inferno_data)    
    
    multimodel_list = np.array([jsbach_data_ma,clm_data_ma,ctem_data_ma,
                            blaze_data_ma,orchidee_data_ma])                       
    multimodel_mean = np.nanmean(multimodel_list, axis=0)
    multimodel_std = np.nanstd(multimodel_list, axis=0)
    
    plt.plot(years[19:], np.divide(jsbach_data_ma,1e12), 'r', label='JSBACH')
    plt.plot(years[19:], np.divide(clm_data_ma,1e12), 'b', label='CLM')
    plt.plot(years[19:], np.divide(ctem_data_ma,1e12), 'g', label='CTEM')
    plt.plot(years[19:], np.divide(blaze_data_ma,1e12), 'm', label='LPJ-GUESS-BLAZE')
    plt.plot(years[19:], np.divide(orchidee_data_ma,1e12), 'y', label='ORCHIDEE')
    plt.plot(years[19:], np.divide(inferno_data_ma,1e12), 'c', label='INFERNO')
    
    plt.plot(years[19:], np.divide(multimodel_mean, 1e12), 'k--', label='Multimodel Mean')
    plt.fill_between(years[19:], np.divide(multimodel_mean-multimodel_std,1e12), 
                    np.divide(multimodel_mean+multimodel_std,1e12), 
                    facecolor='grey', alpha=0.2)
    
    plt.xlabel('Year')
    plt.ylabel('Total Burnt Area (millions of $km^2$)')
    
    plt.legend()
    plt.show()


#plot_past_emissions()
#plot_past_burnt_area()
