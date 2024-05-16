# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 12:19:01 2024

@author: nivfe
"""
import pandas as pd
import geopandas as gpd
import json
import numpy as np
import matplotlib.pyplot as plt
import ast
import csv

def export_tuples(data,file_name):
    '''
    
    '''
    with open(file_name,mode='w',newline='') as file:
        writer=csv.writer(file)
        writer.writerow([['N Repairs','Exceedance rate']])
        for row in data:
            x_str=json.dumps(row[0])
            y_str=json.dumps(row[1])
            writer.writerow([x_str,y_str])
            

def createRR_df(RRFileName):
    #create repair rate dataFrame containing N_PGV,N_PGD and RR for 
    RR_df=pd.read_csv(RRFileName,header=0)
    #load geo Dataframe
    
   #Group RR by global site Id for all the pipes 
    RR_agg=RR_df.groupby(['Global_site_id'])[['event_id','RR','N_PGV','N_PGD']].agg(list).reset_index()
    
    
    RR_agg['event_id']=RR_agg.apply(lambda x:np.array(x['event_id']),axis=1 )
    RR_agg['RR']=RR_agg.apply(lambda x:np.array(x['RR']),axis=1)
    RR_agg['N_PGV']=RR_agg.apply(lambda x:np.array(x['N_PGV']),axis=1)
    RR_agg['N_PGD']=RR_agg.apply(lambda x:np.array(x['N_PGD']),axis=1)
    return RR_agg

def compute_TotalsumFailures(RR_agg,Field_Name):
    TN=np.vstack(np.array(RR_agg[Field_Name].transform(lambda x:np.array(x)))).sum(axis=0)
    return TN
def annualExceedanceX(Vector,X,simmulationTime):
    return [np.sum(np.heaviside(-1*np.ones_like(Vector) *threshold + Vector,0.0))/simmulationTime for threshold in X]
def get_exceedanceTotalFailurestuple(RR_agg,FieldName,simmulationTime):
    TN=compute_TotalsumFailures(RR_agg,FieldName)
    TotalN=np.logspace(-3, np.log10(np.max(np.array(TN))),num=100)
    exceedanceN=annualExceedanceX(TN,TotalN,simmulationTime)
    return(TotalN,exceedanceN)
def compute_exceedance_rate(T, data,simmulationTime):
    """
    Compute the X (RR,N) for a given exceedance rate corresponding to a return period T.

    Parameters:
    - T: Return period (in years) for which to compute the exceedance rate.
    - RR: Array of Repair Rates from seismic data.

    Returns:
    - exceedance_rr: The RR value that corresponds to the exceedance rate for the return period T.
    - exceedance_rate: The calculated exceedance rate for the exceedance_rr.
    """
    # Total number of years in the catalog
    total_years = simmulationTime  # As specified
    
    RR=np.array(data)
    
    # Calculate the exceedance rate for 1/T
    desired_exceedance_rate = 1 / T

    # Sort the RR data in descending order to work from highest to lowest
    sorted_RR = np.sort(RR)[::-1]

    # Convert exceedances to exceedance rates
    exceedance_rates = np.arange(1, len(sorted_RR) + 1) / total_years

    # Find the RR value that corresponds to the closest exceedance rate for 1/T (or directly matches)
    # This is the value where the exceedance rate just drops below the desired rate
    try:
        idx = np.where(exceedance_rates >= desired_exceedance_rate)[0][0]
        exceedance_rr = sorted_RR[idx]
    except:
        exceedance_rr=0

    return exceedance_rr

def get_exceedanceTuple2(T,RR_agg,FieldName,simmulationTime):
    TN=compute_TotalsumFailures(RR_agg,FieldName)
    data_exceedance=np.array([compute_exceedance_rate(Ti,TN,simmulationTime) for Ti in T])
    exceedance=1/np.array(T)
    return (data_exceedance,exceedance)



def create_plots(data, title, x_label, y_label,Field, legends ,log_scale, figsize=(5, 4)):
    """
    Generate a plots based on the provided data.

    Parameters:
    - data: A list of tuples, each containing (x, y) data for plotting.
    - plot_titles: A list of strings for the titles of each plot.
    - x_labels: A list of strings for the x-axis labels of each plot.
    - y_labels: A list of strings for the y-axis labels of each plot.
    - log_scale: A boolean indicating whether to use a logarithmic scale. Default is False.
    - figsize: A tuple indicating the figure size. Default is (10, 10).
    """

    plt.style.use('ggplot')
    plt.figure(figsize=figsize)

    
    for i in range(len(data)):
        x,y=data[i]
        AnnualExpected_N=np.dot(x,1-np.exp(np.array(y)))
        legends[i]=f'{legends[0]} Annual N= {AnnualExpected_N:.3f}'
        plt.plot(x, y, marker='o', linestyle='-')
        
    if log_scale:
        plt.xscale('log')
        #plt.yscale('log')
        plt.xlim(0.1,None)
    else:
        plt.xlim(0,None)
            
    plt.yscale('log')    
    plt.title(title, fontsize=28, fontweight='bold')
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel(x_label, fontsize=24, fontweight='bold')
    plt.ylabel(y_label, fontsize=24, fontweight='bold')
    plt.legend(legends,loc='upper right',fontsize=20)
    plt.grid(True, color='black', linestyle='-', linewidth=1.0)
    plt.grid(True, color='k', linestyle='-', linewidth=0.2,which='minor')
    plt.ylim(0.001,1)
    # plt.text(0.95, 0.95, f'AnnualExpected Failures = {AnnualExpected_N}',
    #          transform=plt.gca().transAxes, ha='right', va='top',
    #          bbox=dict(facecolor='white', alpha=0.8, edgecolor='black'))
    plt.show()
    
    plt.savefig(f'Aggregated_Loss_exceedance2_{Field}.png', dpi=100)
    plt.tight_layout()
    plt.show()
    
    
    
def main():
    RR_FileNames=['RR_results/RR_pipe.csv'
                   # ,'RR_results/RR_pipe_corr_model1.csv'
                   # ,'RR_results/RR_pipe_corr_model4.csv'
                   # ,'RR_results/RR_pipe_corr_model5.csv'
                  ]
    RR_dfs=[pd.read_csv(RR_FileNames[0])
            # , pd.read_csv(RR_FileNames[1])
            # , pd.read_csv(RR_FileNames[2])
            # , pd.read_csv(RR_FileNames[3])
            ]
    simmulationTime=3000.0
    T=[1,2.6,3,5.0,7.5,10.0,25.0,50.0,75.0,100.0,250.0,500.0]
    
    IMLabel='PGD'
    print(f'creating plot_{IMLabel}')
    title=f'Loss Exceedance {IMLabel}'
    y_label=r'Annual Exceedance rate $\lambda$'
    x_label=f'Expected Number of Repairs {IMLabel}'
    FieldName=f'N_{IMLabel}'
    labels=['No Correlation'
            # ,'Goda_Hong_PGV'
            # ,'Wang_Takada_2005'
            ]
    
    
    RR_dfs[0][FieldName]=[eval(x)for x in RR_dfs[0][FieldName]]
    # RR_dfs[1][FieldName]=[eval(x)for x in RR_dfs[1][FieldName]]
    # RR_dfs[2][FieldName]=[eval(x)for x in RR_dfs[2][FieldName]]
    #RR_dfs[3][FieldName]=[eval(x)for x in RR_dfs[3][FieldName]]
           
    data=[get_exceedanceTotalFailurestuple(x,FieldName,simmulationTime) for x in RR_dfs]
    #data=[get_exceedanceTotalFailurestuple(x,FieldName,simmulationTime) for x in RR_dfs[:-1]]
    #data=[get_exceedanceTuple2(T,x,FieldName,simmulationTime) for x in RR_dfs[:-1]] 
    create_plots(data, title, x_label, y_label,FieldName,labels, log_scale=True, figsize=(10, 8))
    
    
    # IMLabel='PGD'
    # print(f'creating plot_{IMLabel}')
    # title=f'Loss Exceedance {IMLabel}'
    # y_label=r'Annual Exceedance rate $\lambda$'
    # x_label=f'Expected Number of Repairs {IMLabel}'
    # FieldName=f'N_{IMLabel}'
    # labels=['No Correlation'
    #         ,'Boore_2003'
    #         ,'Goda_Atk_2009'
    #         ,'J_Baker_2009'
    #         ]
    
    # RR_dfs[0][FieldName]=[eval(x)for x in RR_dfs[0][FieldName]]
    # RR_dfs[1][FieldName]=[eval(x)for x in RR_dfs[1][FieldName]]
    # RR_dfs[2][FieldName]=[eval(x)for x in RR_dfs[2][FieldName]]
    # RR_dfs[3][FieldName]=[eval(x)for x in RR_dfs[3][FieldName]]
    
    # data=[get_exceedanceTotalFailurestuple(x,FieldName,simmulationTime) for x in RR_dfs]
    # #data=[get_exceedanceTuple2(T,x,FieldName,simmulationTime) for x in RR_dfs]
    # create_plots(data, title, x_label, y_label,FieldName,labels, log_scale=True, figsize=(10, 8))
   
    
if __name__:
    main()
