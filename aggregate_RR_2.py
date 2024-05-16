# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 06:49:28 2024

@author: nivfe
"""
import geopandas as gpd
import pandas as pd
import numpy as np
import hazard as hz
import argparse
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import gc
import ast

def get_exceedance(index_list,df):
    exceedance_array=[df['annual_lambda'][index] for index in index_list]
    return exceedance_array

def get_delta_excedance(df):
    sorted_=np.sort(df["exceedance_rate"])
    index=np.argsort(df["exceedance_rate"])
    delta=sorted_-np.append(0,sorted_[1:])
    return 

#def function to create the individual Epf for each element
    
def expectedProbabilityofFailure(data,IMLabel):
    '''
    

    Parameters
    ----------
    data : data frame 
        Data Frame containing the fields object
         and array for number of repair for each event N

    Returns
    -------
    Array containing the probability of failure given an .

    '''
    #Extract the data contiaining info about the number of failures for each pipe
    N=np.vstack(np.array(data[f'N_{IMLabel}'].transform(lambda x:np.array(x))))
    
    pf=1-np.exp(-N)
    Epf=np.average(pf,axis=1)
    return Epf

def failure_Rate(data,IMLabel,simmulationTime):
    '''
    

    Parameters
    ----------
    data : data frame 
        Data Frame containing the fields object
         and array for number of repair for each event N

    Returns
    -------
    Array containing the probability of failure given an .

    '''
    #Extract the data contiaining info about the number of failures for each pipe
    N=np.vstack(np.array(data[f'N_{IMLabel}'].transform(lambda x:np.array(x))))
    
    pf=1-np.exp(-N)
    fRate=np.sum(pf,axis=1)*1/simmulationTime
    return fRate


def computeExcedanceVector(vector,data,simmulationTime):
    '''
    

    Parameters
    ----------
    vector : TYPE
        DESCRIPTION.
    data : TYPE
        DESCRIPTION.
    simmulationTime : TYPE
        DESCRIPTION.

    Returns
    -------
    excedanceVector : TYPE
        DESCRIPTION.

    '''
    excedanceVector=[np.sum(np.heaviside(-1*np.ones_like(data) *threshold + data, 0.0))/simmulationTime for threshold in vector]
    
    return excedanceVector
def compute_X_for_Excedance(T,data,simmulationTime):
    '''
    returns the value corresponding to the given return period T

    Parameters
    ----------
    T : float
        return period.
    data : list, array
        contains vector for which the correspondance excedance willl be computed.
    
    simmulationTime : 
        simmulation time to compute excedance vector.

    Returns
    -------
    expectedP : float
        expecte value X  for the required excedance.

    '''
    vector=np.linspace(0,max(np.array(data)))
    #compute_excedance Vector for the the input data
    excedanceVector=computeExcedanceVector(vector, data,simmulationTime)
    
    expectedX=log_interpolation(1/T,vector,excedanceVector)
    
    
    return expectedX

def compute_exceedance_rate(T, data,simmulationTime):
    """
    Compute the Repair Rate (RR) for a given exceedance rate corresponding to a return period T.

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

def log_interpolation(x,X,Y):
    '''
    


    Parameters
    ----------
    yx : TYPE
        Probability of exceedance to be found.
    X : TYPE
        Exceedance Vector.
    Y : TYPE
        Vector.

    Returns
    -------
    xx : TYPE
        DESCRIPTION.

    '''
    log_X=np.log(X)#exceedanceVector
    log_X[log_X==-np.inf]=-10000
    log_fp=np.log(Y)#Vector
    log_fp[log_fp==-np.inf]=-10000
    xx=np.exp(np.interp(np.log(x),log_X,log_fp))
    return xx
def computeNfailures(N_PGV,N_PGD):
    '''
    

    Parameters
    ----------
    N_PGV : TYPE
        DESCRIPTION.
    N_PGD : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    return np.array(N_PGV)*0.6+0.4*np.array(N_PGD)

def parsePipeCosttoDict(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    prices_dict = {}
    
    for diameter in root.findall("./Diameters/Diameter"):
        diameter_value = diameter.get("value")
        materials_dict = {}
        for material in diameter.find("Materials").findall("Material"):
            material_name = material.get("name")
            unit_price = material.get("unitPrice")
            materials_dict[material_name] = float(unit_price)  # Convert string to float for price
        prices_dict[diameter_value] = materials_dict
        
    return prices_dict

def get_price(prices_dict, diameter, material):
    """
    Get the price for a given diameter and material from the dictionary.
    
    :param prices_dict: Dictionary containing the prices information
    :param diameter: Diameter of the pipe
    :param material: Material of the pipe
    :return: Price or None if not found
    """
    diameter_str = str(diameter)  # Ensure diameter is a string, matching dictionary keys
    material_price = prices_dict.get(diameter_str, {}).get(material)
    return material_price

def expectedLoss(xml_file,data):
    '''
    

    Parameters
    ----------
    xml_file : TYPE
        DESCRIPTION.
    data : TYPE
        DESCRIPTION.

    Returns
    -------
    expectedLoss : TYPE
        Returns expected loss of pipe.

    '''
    
    return expectedLoss

def main():
    # #Here I need to extend for parsing the inputs
    parser = argparse.ArgumentParser(description='Process simulation parameters')
    parser.add_argument('--RRFileName', type=str, help='Path to file containing RR and N for each event of PGV')
    
    #parser.add_argument('--GeoJSONPAth',type=str,help='Path to file containing GeoJson')
    # parser.add_argument('--RRPGVFileName', type=str, help='Path to file containing RR and N for each event of PGD')
    # #parser.add_argument('--meshFileName', type=str, help='Path to the mesh file')
    # #parser.add_argument('--IMtype',type=str,help='IM type that you wish to analyze')
    # #parser.add_argument('--simmulationTimeSpan', type=float, help='Simulation time span')
    # #parser.add_argument('--timePeriod', type=float, help='Time period for hazard curves')
    parser.add_argument('--exportName',type=str,help='Path or Filename of the desired GeoJson with no format')
    # #parser.add_argument('--sig_eps',type=str,help='Path or Filename of the inter event residuals sig_eps')
    args = parser.parse_args()
    #---------------------------------------------
    print('Loading Data')
    #RRFileName='RR_pipe_corr_3_16.csv'
    RRFileName=args.RRFileName
    RR_agg=pd.read_csv(RRFileName,header=0)
    #outputFile=args.csvName
    
    geoJSONFile='extended_data.geojson'
    #events_exceedance='events_annual_rate_66.csv'
    
    #events_df=pd.read_csv(events_exceedance,header=0).set_index('event_id')
    
    geo_df = gpd.read_file(geoJSONFile)
    geo_df['Global_site_id']=geo_df.apply(
        lambda row: row['GLOBALID']+'_'+str(row['site_id']),axis=1)
    '''
    RR_agg=RR_df.groupby(['Global_site_id'])[['event_id',
                                                      'RR',
                                                      'N_PGV','N_PGD']].agg(list).reset_index()
   '''
    geo_df=geo_df.merge(RR_agg[['Global_site_id']],on='Global_site_id')
    #geo_df['exceedance_rate']=geo_df.apply(lambda row: get_exceedance(np.array(row['event_id']),events_df),axis=1)
    geo_df['RR_PGV']=[eval(x)for x in RR_agg['RR_PGV']]
    geo_df['RR_PGD']=[eval(x)for x in RR_agg['RR_PGD']]
    geo_df['N_PGV']=[eval(x)for x in RR_agg['N_PGV']]
    geo_df['N_PGD']=[eval(x)for x in RR_agg['N_PGD']]
    
    #RR_agg.apply(lambda x:np.array(x['RR']),axis=1)
   
    del RR_agg
    model=args.exportName
    gc.collect()
    T=np.array([1.0,5.0,10.0,50.0,100.0,250.0,500.0])
    simmulationTime=30000.0
    
    print('Epf')
    #geo_df['N_tot']=geo_df.apply(lambda row:computeNfailures(row['N_PGV'],row['N_PGD']),axis=1)
    #compute pf given there is a seismic event
    geo_df['Epf_PGD']=expectedProbabilityofFailure(geo_df,'PGD')

    geo_df['Epf_PGV']=expectedProbabilityofFailure(geo_df,'PGV')
    #geo_df['Epf_tot']=expectedProbabilityofFailure(geo_df,'tot')
    
    geo_df['failRate_PGD']=failure_Rate(geo_df,'PGD',simmulationTime)
    geo_df['failRate_PGV']=failure_Rate(geo_df,'PGV',simmulationTime)
    #get Failure Rate  np.sum((1-np.exp(-N))*1/simmulationTime)
    
    
    keys=['Global_site_id'
          ,'MATERIAL'
          ,'DIAMETER'
          ,'Length'
          ,'Epf_PGD'
          ,'Epf_PGV'
          ,'failRate_PGD'
          ,'failRate_PGV'
          ]
    for T_i in T:
        geo_df[f'Pf_PGD_T_{T_i}']=1-np.exp(-1*np.array(geo_df['failRate_PGD'])*T_i)
        keys=np.append(keys,f'Pf_PGD_T_{T_i}')
        geo_df[f'Pf_PGV_T_{T_i}']=1-np.exp(-1*np.array(geo_df['failRate_PGV'])*T_i)
        keys=np.append(keys,f'Pf_PGV_T_{T_i}')
        geo_df[f'Pf_T_{T_i}']=1-np.exp(-1*(np.array(geo_df['failRate_PGD'])+np.array(geo_df['failRate_PGV']))*T_i)
        keys=np.append(keys,f'Pf_T_{T_i}')

    keys=np.append(keys,'geometry')
    
    gpd.GeoDataFrame(geo_df[keys],crs="EPSG:4326").to_file(f'Results_2_{model}', driver='GeoJSON')
    
    
    
    # print('RR')
    # label='RR_PGV'
    # keys=['Global_site_id']
    # for index in T:
    #     # geo_df[f'{label}_{index}_yr']=geo_df.apply(
    #     #     lambda row:(compute_X_for_Excedance(
    #     #         index,row[f'{label}'],
    #     #         simmulationTime)),axis=1)
    #     # keys=np.append(keys,f'{label}_{index}_yr')
    #     geo_df[f'{label}_{index}_yr']=geo_df.apply(
    #         lambda row:(compute_exceedance_rate(
    #             index,row[f'{label}'],
    #             simmulationTime)),axis=1)
    #     keys=np.append(keys,f'{label}_{index}_yr')
        
    # keys=np.append(keys,'geometry')
    
    # gpd.GeoDataFrame(geo_df[keys],crs="EPSG:4326").to_file(f'Results_2_{label}_{model}', driver='GeoJSON')
   
    # label='RR_PGD'
    # keys=['Global_site_id']
    # for index in T:
    #     # geo_df[f'{label}_{index}_yr']=geo_df.apply(
    #     #     lambda row:(compute_X_for_Excedance(
    #     #         index,row[f'{label}'],
    #     #         simmulationTime)),axis=1)
    #     # keys=np.append(keys,f'{label}_{index}_yr')
    #     geo_df[f'{label}_{index}_yr']=geo_df.apply(
    #         lambda row:(compute_exceedance_rate(
    #             index,row[f'{label}'],
    #             simmulationTime)),axis=1)
    #     keys=np.append(keys,f'{label}_{index}_yr')
        
    # keys=np.append(keys,'geometry')
    
    # gpd.GeoDataFrame(geo_df[keys],crs="EPSG:4326").to_file(f'Results_2_{label}_{model}', driver='GeoJSON')
    # #
   
    # print('N')
    
    # label='N_PGV'
    # keys=['Global_site_id']
    # for index in T:
    #     # geo_df[f'{label}_{index}_yr']=geo_df.apply(
    #     #     lambda row:(compute_X_for_Excedance(
    #     #         index,row[f'{label}'],
    #     #         simmulationTime)),axis=1)
    #     # keys=np.append(keys,f'{label}_{index}_yr')
    #     geo_df[f'{label}_{index}_yr']=geo_df.apply(
    #         lambda row:(compute_exceedance_rate(
    #             index,row[f'{label}'],
    #             simmulationTime)),axis=1)
    #     keys=np.append(keys,f'{label}_{index}_yr')
        
    # keys=np.append(keys,'geometry')
    
    # gpd.GeoDataFrame(geo_df[keys],crs="EPSG:4326").to_file(f'Results_2_{label}_{model}', driver='GeoJSON')
    
    # label='N_PGD'
    # keys=['Global_site_id']
    # for index in T:
    #     # geo_df[f'{label}_{index}_yr']=geo_df.apply(
    #     #     lambda row:(compute_X_for_Excedance(
    #     #         index,row[f'{label}'],
    #     #         simmulationTime)),axis=1)
    #     # keys=np.append(keys,f'{label}_{index}_yr')
    #     geo_df[f'{label}_{index}_yr']=geo_df.apply(
    #         lambda row:(compute_exceedance_rate(
    #             index,row[f'{label}'],
    #             simmulationTime)),axis=1)
    #     keys=np.append(keys,f'{label}_{index}_yr')
        
    # keys=np.append(keys,'geometry')
    
    # gpd.GeoDataFrame(geo_df[keys],crs="EPSG:4326").to_file(f'Results_2_{label}_{model}', driver='GeoJSON')
    
    
    
   
    
    
  
    
    # # Expected Failure probability for different T
    # print('pf fot differnt T')
    # IMT='PGD'
    # label='pf_{IMT}'
    
    # for index in T:
    #     geo_df[f'{label}_{index}_yr']=geo_df.apply(
    #         lambda row:1-np.exp(-compute_X_for_Excedance(
    #             index,row[f'N_{IMT}'],
    #             simmulationTime)),axis=1)
    #     eys=np.append(keys,f'{label}_{index}_yr')
        
    # keys=np.append(keys,'geometry')
    
    # gpd.GeoDataFrame(geo_df[keys],crs="EPSG:4326").to_file(f'Results_{label}_{model}', driver='GeoJSON')
    
    # IMT='PGV'
    # label='pf_{IMT}'
    
    # for index in T:
    #     geo_df[f'{label}_{index}_yr']=geo_df.apply(
    #         lambda row:1-np.exp(-compute_X_for_Excedance(
    #             index,row[f'N_{IMT}'],
    #             simmulationTime)),axis=1)
    #     eys=np.append(keys,f'{label}_{index}_yr')
        
    # keys=np.append(keys,'geometry')
    
    # gpd.GeoDataFrame(geo_df[keys],crs="EPSG:4326").to_file(f'Results_{label}_{model}', driver='GeoJSON')
    
    
    
    
    
    #get expected replacement cost
    '''
    xml_file="Expected cost.xml"
    
    costDict=parsePipeCosttoDict(xml_file)
    
    geo_df['unit_Price']=geo_df.apply(
        lambda row:get_price(costDict, row['DIAMETER'], row['MATERIAL']),axis=1)
   
    #Failure Rate for pipe
    
    geo_df['failure_rate_PGV']=geo_df.apply(
        lambda row:np.dot(1-np.exp(row['N_tot']),row['exceedance_rate'] ),axis=1)
    
    geo_df['failure_rate_PGD']=geo_df.apply(
        lambda row:np.dot(1-np.exp(row['N_tot']),row['exceedance_rate'] ),axis=1)
    #Expected Repair Rate for different T
    '''

    
    
    
    '''
    #expected loss for different T
    label='EL_USD'
    
    for index in T:
        geo_df[f'EL_USD_{index}_yr']=geo_df.apply(
            lambda row:row['unit_Price']*row['Length']*row[f'pF_{index}_yr'],axis=1)
    keys=['Global_site_id',f'pF_{T[0]}_yr',f'pF_{T[1]}_yr',
          f'pF_{T[2]}_yr',f'pF_{T[3]}_yr',f'pF_{T[4]}_yr','geometry']
    
    geo_df[keys].to_file('Example_{label}_{model}_T', driver='GeoJSON')
    
    
    
    #Compute the pnf

    '''
    
    print('END')
    
if __name__ == '__main__':
    
    '''
    *python3.8 aggregate_RR_2.py --RRFileName RR_results/RR_pipe.csv --exportName NoCorr
    *python3.8 aggregate_RR_2.py --RRFileName RR_results/RR_pipe_corr_model1.csv --exportName model1
    python3.8 aggregate_RR.py --RRFileName RR_results/RR_pipe_corr_model2.csv --exportName model2
    python3.8 aggregate_RR.py --RRFileName RR_results/RR_pipe_corr_model3.csv --exportName model3   
    *python3.8 aggregate_RR_2.py --RRFileName RR_results/RR_pipe_corr_model4.csv --exportName model4
    *python3.8 aggregate_RR_2.py --RRFileName RR_results/RR_pipe_corr_model5.csv --exportName model5
    python3.8 aggregate_RR.py --RRFileName RR_results/RR_pipe_corr_model6.csv --exportName model6
    '''
    
    main() 
    
    
    #--RRPGVFileName OBJECT_RR_PGV.csv --RRPGDFileName OBJECT_RR_PGD.csv --exportName Example.csv
    #--RRPGVFileName RR_PGV_pipe_corr.csv --RRPGDFileName RR_PGD_pipe_corr.csv --exportName Risk_Network_corr
    #--RRPGVFileName RR_PGV_pipe.csv --RRPGDFileName RR_PGD_pipe.csv --exportName Risk_Network
