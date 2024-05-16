# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 04:55:19 2024

@author: nivfe
"""

import pandas as pd
import numpy as np
import argparse

def annualExceedanceX(Vector,X,simmulationTime,Index):
    return [np.sum(np.heaviside(-1*np.ones_like(Vector) *threshold + Vector,0.0))/simmulationTime for threshold in X]
def compute_poe_vectorized(rate, t):
    """
    Compute the Probability of Exceedance (POE) for a given time period.

    Parameters:
        rate (float or numpy array): Exceedance rate(s) for the IM.
        t (float): Time period.

    Returns:
        float or numpy array: Probability of Exceedance.
    """
    poe = 1 - np.exp(-rate * t)
    return poe

# Function to get IM for a given return period in a vectorized way
def logarithmic_interpolation(x_list, y_list, x_new):
    """
    Performs logarithmic interpolation to find the new y value at x_new.
    
    Parameters:
    - x_list: List of x values.
    - y_list: List of y values.
    - x_new: The new x value to interpolate the y value for.
    
    Returns:
    - y_new: The interpolated y value at x_new.
    """
    # Ensure inputs are numpy arrays for vectorized operations
    x_array = np.array(x_list)
    y_array = np.array(y_list)

    # Sort the x_array and y_array based on x_array to ensure increasing order
    sorted_indices = np.argsort(x_array)
    x_sorted = x_array[sorted_indices]
    y_sorted = y_array[sorted_indices]
    x_sorted[x_sorted<=0]=1/30000

    # Use logarithmic interpolation only if it makes sense (x > 0 and y > 0)
    if np.any(x_sorted <= 0) or np.any(y_sorted <= 0):
        raise ValueError("Logarithmic interpolation requires all x and y values to be > 0.")
    
    # Convert x and y to log-space since we're doing logarithmic interpolation
    log_x = np.log(x_sorted)
    log_y = np.log(y_sorted)
    

    
    # Interpolate in log-space
    log_x_new = np.log(x_new)
    log_y_new = np.interp(log_x_new, log_x, log_y)
    
    # Convert the interpolated value back from log-space
    y_new = np.exp(log_y_new)
    
    return y_new

def main():
    # #Here I need to extend for parsing the inputs
    parser = argparse.ArgumentParser(description='Process simulation parameters')
    parser.add_argument('--IMFileName', type=str, help='Path to file containing RR and N for each event of PGV')
    parser.add_argument('--SimID',type=str,help='job ID from OQ you wish to analyze')
    parser.add_argument('--simmulationTimeSpan', type=float, help='Simulation time span')
    parser.add_argument('--exportName',type=str,help='Path or Filename of the desired GeoJson with no format')
    args = parser.parse_args()
    jobID=args.SimID 
    rupturesFileName=f"ruptures_{jobID}.csv"
    eventsFileName=f"events_{jobID}.csv"
    events_df=pd.read_csv(eventsFileName,header=1)
    ruptures_df=pd.read_csv(rupturesFileName,header=1)
    IMFileName=args.IMFileName
    IM_df=pd.read_csv(IMFileName,header=0)
    meshFileName=f'sitemesh_{jobID}.csv'
    mesh_df=pd.read_csv(meshFileName,header=1)
    simmulationTime=args.simmulationTimeSpan
    model=args.exportName
    
    print('Merging Data Frames...')
    IM_merged_df=pd.merge(IM_df,mesh_df,on='site_id')
    del IM_df
    
    #Group by location and aggregate 'IM' to a list 
    print('organizing data frame in aggregated lists...')
    IM_agg_df=IM_merged_df.groupby(['site_id','lon','lat'])[['event_id','gmv_PGA','gmv_PGV','gmv_PGD']].agg(list).reset_index()
    
    merged_df = pd.merge(events_df, ruptures_df, on='rup_id', how='left')
    merged_df = merged_df[['event_id','mag','multiplicity',
                           'centroid_lon','centroid_lat']]
    
    
    
    #Create Exceedance curve for each IM 
    IM_agg_df['PGD_exceedance']=IM_agg_df.apply(lambda row:
                                                annualExceedanceX(np.array(row['gmv_PGD']),
                                                                  np.logspace(-6, np.log10(np.max(np.array(row['gmv_PGD']))),num=100),simmulationTime,
                                                                  row['site_id']),axis=1)
    IM_agg_df['PGV_exceedance']=IM_agg_df.apply(lambda row:
                                                annualExceedanceX(np.array(row['gmv_PGV']),
                                                                  np.logspace(-6, np.log10(np.max(np.array(row['gmv_PGV']))),num=100),simmulationTime,
                                                                  row['site_id']),axis=1)
    IM_agg_df['PGA_exceedance']=IM_agg_df.apply(lambda row:
                                                annualExceedanceX(np.array(row['gmv_PGA']),
                                                                  np.logspace(-6, np.log10(np.max(np.array(row['gmv_PGA']))),num=100),simmulationTime,
                                                                  row['site_id']),axis=1)
    #IM_agg_df[['site_id','exceedance_rate_PGD','exceedance_rate_PGV','exceedance_rate_PGA']].to_csv('Excedance_66.csv',index=False)
    
    T=np.array([1.0,10.0,50.0,100.0,500.0])
    keys=['site_id','lon','lat']
    for index in T:
        IM_agg_df[f'PGD_POE_{index}_yr']=IM_agg_df.apply(
            lambda row: logarithmic_interpolation(compute_poe_vectorized(np.array(row['PGD_exceedance']),1.0),
                                                  np.logspace(-6, np.log10(np.max(np.array(row['gmv_PGD']))),num=100),
                                                  1/index),axis=1)
        IM_agg_df[f'PGV_POE_{index}_yr']=IM_agg_df.apply(
            lambda row: logarithmic_interpolation(compute_poe_vectorized(np.array(row['PGV_exceedance']),1.0),
                                                  np.logspace(-6, np.log10(np.max(np.array(row['gmv_PGV']))),num=100),
                                                  1/index),axis=1)
        IM_agg_df[f'PGA_POE_{index}_yr']=IM_agg_df.apply(
            lambda row: logarithmic_interpolation(compute_poe_vectorized(np.array(row['PGA_exceedance']),1.0),
                                                  np.logspace(-6, np.log10(np.max(np.array(row['gmv_PGA']))),num=100),
                                                  1/index),axis=1)
    IM_agg_df[np.append(keys,[f'PGD_POE_{index}_yr' for index in T])].to_csv(f'RR_Results/PGD_16_T_{model}_yr.csv',index=False)
    IM_agg_df[np.append(keys,[f'PGV_POE_{index}_yr' for index in T])].to_csv(f'RR_Results/PGV_16_T_{model}_yr.csv',index=False)
    IM_agg_df[np.append(keys,[f'PGA_POE_{index}_yr' for index in T])].to_csv(f'RR_Results/PGA_16_T_{model}_yr.csv',index=False)
     
    
    
    print('HI')
if __name__ == '__main__':
    main()   
    '''
    python3.8 IM_hazard_corr.py --IMFileName gmf_corr/complete_gmf-data_16.csv --exportName NoCorr --simmulationTime 30000.0 --SimID 16
    python3.8 IM_hazard_corr.py --IMFileName gmf_corr/complete_gmf-data_corr_16_PGA_Boore_2003_PGV_Goda_Hong_PGV.csv --exportName model1 --simmulationTime 30000.0 --SimID 16
    python3.8 IM_hazard_corr.py --IMFileName gmf_corr/complete_gmf-data_corr_16_PGA_Boore_2003_PGV_Wang_Takada_2005.csv --exportName model2 --simmulationTime 30000.0 --SimID 16
    python3.8 IM_hazard_corr.py --IMFileName gmf_corr/complete_gmf-data_corr_16_PGA_Goda_Atk_2009_PGV_Goda_Hong_PGV.csv --exportName model3 --simmulationTime 30000.0 --SimID 16
    python3.8 IM_hazard_corr.py --IMFileName gmf_corr/complete_gmf-data_corr_16_PGA_Goda_Atk_2009_PGV_Wang_Takada_2005.csv --exportName model4 --simmulationTime 30000.0 --SimID 16
    python3.8 IM_hazard_corr.py --IMFileName gmf_corr/complete_gmf-data_corr_16_PGA_J_Baker_2009_PGV_Goda_Hong_PGV.csv --exportName model5 --simmulationTime 30000.0 --SimID 16
    python3.8 IM_hazard_corr.py --IMFileName gmf_corr/complete_gmf-data_corr_16_PGA_J_Baker_2009_PGV_Wang_Takada_2005.csv --exportName model6 --simmulationTime 30000.0 --SimID 16
    '''