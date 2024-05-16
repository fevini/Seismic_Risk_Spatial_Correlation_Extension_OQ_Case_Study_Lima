# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 21:03:27 2024

@author: nivfe
"""
import geopandas as gpd
import pandas as pd
import numpy as np
import argparse
import ast
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import json

def RR_Fragility_PGV(PGV, K):
    '''


    Parameters
    ----------
    PGV : Double
        PGV IM .
    K : Double
        Constant, of pipe.

    Returns
    -------
    Double
        Repair Rate associated to Transient wave n/km.

    '''
    PGV = PGV.T
    return [0.002416*np.multiply(K, PGV_i) for PGV_i in PGV]


def PGV_from_RR(RR, K):
    '''


    Parameters
    ----------
    RR : TYPE
        DESCRIPTION.
    K : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    return RR/(0.002416*K)


def RR_Fragility_PGD(PGD, K):
    '''


    Parameters
    ----------
    PGD : TYPE
        PGD IM.
    K : TYPE
        Constant associated to pipe.

    Returns
    -------
    TYPE
        Repair rate asssociated to ground deformation n/km.

    '''
    PGD = PGD.T
    return [2.5829*np.multiply(K, np.power(PGD_i, 0.319)) for PGD_i in PGD]


def RR_PGD(RR, K):
    '''


    Parameters
    ----------
    RR : TYPE
        DESCRIPTION.
    K : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    return np.emath.logn(0.319, RR/(K*2.5829))


def RRfromN(N, L):
    '''


    Parameters
    ----------
    N : Array
        N repairs per km.
    L : float
        lenght of pipe in m.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    return N/L


def extract_site_id(data):

    # Use a set to store unique events
    unique_events = set()

    # Iterate through the dataframe to extract events
    for event_list in data['site_id']:
        unique_events.update(event_list)

    # Convert the set to a sorted list
    all_events = sorted(list(unique_events))
    return all_events


def get_exceedance(index_list, df):
    '''
    only to get annual lambda for each event  with the rupture data frame 
    with attribute 'annual rate'.
    Parameters
    ----------
    index_list : TYPE
        list containing the name of the events.
    df : TYPE
        rupture data frame with  annual exceedance by event. .

    Returns
    -------
    exceedance_array : TYPE
        DESCRIPTION.

    '''
    exceedance_array = [df['annual_lambda'][index] for index in index_list]
    return exceedance_array


# def function to create the individual Epf for each element
def plotExcedance(X, Y):
    '''


    Parameters
    ----------
    X : TYPE
        DESCRIPTION.
    Y : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    plt.yscale('log')
    # plt.style.use('seaborn-whitegrid')
    plt.plot(X, Y)


def expectedProbabilityofFailure(data, IMLabel):
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
    # Extract the data contiaining info about the number of failures for each pipe
    N = np.vstack(
        np.array(data[f'N_{IMLabel}'].transform(lambda x: np.array(x))))

    pf = 1-np.exp(-N)
    Epf = np.average(pf, axis=1)
    return Epf


def computeAnnualOcurrence(vector, simmulationTimeSpan):
    '''


    Parameters
    ----------
    vector : TYPE
        DESCRIPTION.
    simmulationTimeSpan : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    '''
    return np.sum(vector, axis=0)/simmulationTimeSpan


def computeExcedanceVector(vector, data, simmulationTime):
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
    excedanceVector = [np.sum(np.heaviside(-1*np.ones_like(data) *
                              threshold + data, 0.0))/simmulationTime for threshold in vector]

    return excedanceVector


def computeDeltaofVector(Vector):
    # Convert list to a numpy array if it's not already

    # Step 1: Get the sorted indices
    sorted_indices = np.argsort(Vector)

    # Step 2: Sort the vector
    sorted_vector = Vector[sorted_indices]

    # Step 3: Calculate deltas in the sorted order
    # Use np.diff to calculate the difference between consecutive elements in the sorted vector
    sorted_deltas = np.diff(sorted_vector)

    # Step 4: Create a placeholder for deltas with an appropriate size
    # Initialize with np.nan or any other placeholder value since the last element won't have a delta
    deltas_with_original_order = np.full_like(Vector, np.nan, dtype=np.float64)

    # Step 5: Place the sorted deltas back into the positions corresponding to the first element
    # of each pair in their original order. We use the sorted indices to map back.
    # Note that the last delta corresponds to the second last element in sorted order,
    # so we use sorted_indices[:-1] to correctly place deltas.
    deltas_with_original_order[sorted_indices[:-1]] = sorted_deltas

    return deltas_with_original_order


def compute_X_for_Excedance(T, data, simmulationTime):
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
    vector = np.linspace(0, max(np.array(data)))
    # compute_excedance Vector for the the input data
    excedanceVector = computeExcedanceVector(vector, data, simmulationTime)

    expectedX = log_interpolation(1/T, vector, excedanceVector)

    return expectedX


def log_interpolation(yx, X, Y):
    '''



    Parameters
    ----------
    yx : TYPE
        DESCRIPTION.
    X : TYPE
        DESCRIPTION.
    Y : TYPE
        DESCRIPTION.

    Returns
    -------
    xx : TYPE
        DESCRIPTION.

    '''
    xx = np.exp(np.interp(np.log(yx), np.log(X), np.log(Y)))
    return xx


def computeNfailures(N_PGV, N_PGD):
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
            # Convert string to float for price
            materials_dict[material_name] = float(unit_price)
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
    diameter_str = str(
        diameter)  # Ensure diameter is a string, matching dictionary keys
    material_price = prices_dict.get(diameter_str, {}).get(material)
    return material_price


def expectedLoss(xml_file, data):
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

# -------------------------------------------


# -------------------------------------------


def main():
    # Here I need to extend for parsing the inputs
    parser = argparse.ArgumentParser(
        description='Process simulation parameters')
    parser.add_argument('--IMFileName', type=str, help='Path to the IM file')
    #parser.add_argument('--PGDFileName', type=str, help='Path to the PGV file')
    #parser.add_argument('--meshFileName', type=str, help='Path to the mesh file')
    #parser.add_argument('--IMtype',type=str,help='IM type that you wish to analyze')
    parser.add_argument('--simmulationTimeSpan', type=float,
                        help='Simulation time span')
    #parser.add_argument('--timePeriod', type=float, help='Time period for hazard curves')
    parser.add_argument('--csvName', type=str,
                        help='Path or Filename of the desired Hazard curve')
    parser.add_argument('--corrFlag',type=str,
                        help='Flag stating if file is spatialy correlated')
    #parser.add_argument('--sig_eps',type=str,help='Path or Filename of the inter event residuals sig_eps')
    # # Parse the command-line arguments
    args = parser.parse_args()
# ----------------------------------------º-----
    # IMFileName_PGV='gmf-data_corr_66.csv'
    IMFileName_PGV = args.IMFileName
    IM_df = pd.read_csv(IMFileName_PGV, header=0)
    meshFileName = 'sitemesh_66.csv'
    mesh_df = pd.read_csv(meshFileName, header=1)
    # simmulationTime=3000.0

    IMLabel_PGV = 'PGV'
    IMLabel_PGD = 'PGD'
    if args.corrFlag=='True':
        IMCol_PGV = f'gmv_{IMLabel_PGV}_corr'
        IMCol_PGD = f'gmv_{IMLabel_PGD}_corr'
    else:
        IMCol_PGV = f'gmv_{IMLabel_PGV}'
        IMCol_PGD = f'gmv_{IMLabel_PGD}'
    
    geoJSONFile = 'extended_data.geojson'

    # outputName='OBJECT_RR.csv'
    outputName = args.csvName

    geo_df = gpd.read_file(geoJSONFile)

    print('Merging Data Frames...')
    IM_merged_df = pd.merge(IM_df, mesh_df, on='site_id')
    IM_merged_df.set_index(['site_id', 'event_id'], inplace=True)

    PGV_pivot_event_df = IM_merged_df.pivot_table(
        index='site_id', columns='event_id', values=[IMCol_PGV])
    PGV_pivot_event_df = PGV_pivot_event_df.fillna(0)
    PGD_pivot_event_df = IM_merged_df.pivot_table(
        index='site_id', columns='event_id', values=[IMCol_PGD])
    PGD_pivot_event_df = PGD_pivot_event_df.fillna(0)
    print('aggregating events')
    All_events_list = PGV_pivot_event_df.columns.get_level_values(
        'event_id').unique().to_numpy()


    
    #use vectorized operations for the computation of RR and N at each opject.
    PGV_event = PGV_pivot_event_df.xs(5, axis=1, level='event_id').reset_index()
   
    K_df = geo_df[['GLOBALID', 'k1', 'k2', 'site_id', 'Length','geometry']].merge(
        PGV_event['site_id'], on='site_id')
    K_df['Global_site_id'] = K_df.apply(lambda x:
                                        x['GLOBALID']+'_'+str(x['site_id']), axis=1)
    
    index = np.array(K_df['site_id'])
    PGV = PGV_pivot_event_df[IMCol_PGV].values[index]
    PGD = PGD_pivot_event_df[IMCol_PGD].values[index]
    # rows objects columns events
    RR_PGV = np.vstack(RR_Fragility_PGV(PGV, K_df['k1'].values)).T
    RR_PGD = np.vstack(RR_Fragility_PGD(PGD, K_df['k2'].values)).T
    N_PGV = np.vstack([np.multiply(K_df['Length'].values,
                      RR_PGV_i)/1000.0 for RR_PGV_i in RR_PGV.T]).T
    N_PGD = np.vstack([np.multiply(K_df['Length'].values,
                      RR_PGD_i)/1000.0 for RR_PGD_i in RR_PGD.T]).T
    K_df['event_id']=K_df.apply(lambda x: All_events_list.tolist(), axis=1)
    K_df['RR_PGV']=RR_PGV.tolist()

    K_df['RR_PGD']=RR_PGD.tolist()

    K_df['N_PGD']=N_PGD.tolist()

    K_df['N_PGV']=N_PGV.tolist()
    
    
    keys=['Global_site_id','site_id','event_id','Length','RR_PGV','RR_PGD','N_PGV','N_PGD']
    K_df[keys].to_csv(outputName,index=False)
    


if __name__ == '__main__': 
    '''
    example
    python3.8 risk_2.py --IMFileName gmf-data_66_complete.csv --simmulationTimeSpan 3000.0 --csvName RR_results/RR_pipe.csv --corrFlag False
    
    python3.8 risk_2.py --IMFileName gmf_corr/gmf-data_16_complete.csv --simmulationTimeSpan 30000.0 --csvName RR_results/RR_pipe.csv --corrFlag False
    python3.8 risk_2.py --IMFileName gmf_corr/complete_gmf-data_corr_16_PGA_Boore_2003_PGV_Goda_Hong_PGV.csv  --simmulationTimeSpan 30000.0 --csvName RR_results/RR_pipe_corr_model1.csv --corrFlag False
    python3.8 risk_2.py --IMFileName gmf_corr/complete_gmf-data_corr_16_PGA_Boore_2003_PGV_Wang_Takada_2005.csv  --simmulationTimeSpan 30000.0 --csvName RR_results/RR_pipe_corr_model2.csv --corrFlag False
    python3.8 risk_2.py --IMFileName gmf_corr/complete_gmf-data_corr_16_PGA_Goda_Atk_2009_PGV_Goda_Hong_PGV.csv  --simmulationTimeSpan 30000.0 --csvName RR_results/RR_pipe_corr_model3.csv --corrFlag False
    python3.8 risk_2.py --IMFileName gmf_corr/complete_gmf-data_corr_16_PGA_Goda_Atk_2009_PGV_Wang_Takada_2005.csv  --simmulationTimeSpan 30000.0 --csvName RR_results/RR_pipe_corr_model4.csv --corrFlag False
    python3.8 risk_2.py --IMFileName gmf_corr/complete_gmf-data_corr_16_PGA_J_Baker_2009_PGV_Goda_Hong_PGV.csv  --simmulationTimeSpan 30000.0 --csvName RR_results/RR_pipe_corr_model5.csv --corrFlag False
    python3.8 risk_2.py --IMFileName gmf_corr/complete_gmf-data_corr_16_PGA_J_Baker_2009_PGV_Wang_Takada_2005.csv  --simmulationTimeSpan 30000.0 --csvName RR_results/RR_pipe_corr_model6.csv --corrFlag False
    '''
    main() 

'''
gmf_corr/complete_gmf-data_corr_16_PGA_Boore_2003_PGV_Goda_Hong_PGV.csv
gmf_corr/complete_gmf-data_corr_16_PGA_Boore_2003_PGV_Wang_Takada_2005.csv

gmf_corr/complete_gmf-data_corr_16_PGA_Goda_Atk_2009_PGV_Goda_Hong_PGV.csv
gmf_corr/complete_gmf-data_corr_16_PGA_Goda_Atk_2009_PGV_Wang_Takada_2005.csv

gmf_corr/complete_gmf-data_corr_16_PGA_J_Baker_2009_PGV_Goda_Hong_PGV.csv
gmf_corr/complete_gmf-data_corr_16_PGA_J_Baker_2009_PGV_Wang_Takada_2005.csv

'''    
