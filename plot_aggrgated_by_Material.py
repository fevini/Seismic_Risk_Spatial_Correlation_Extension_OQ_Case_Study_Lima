# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 12:03:33 2024

@author: nivfe
"""
#rutina para unir los resultados de el ana+àlisis de riesgo de las tuberias para mostrar 
import geopandas as gpd
import pandas as pd 
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import argparse
def replace_material_labels(material):
    replacements = {'HD': 'DI', 'FV': 'GF', 'CAN': 'AR', 'FOFO': 'CI'}
    if material == 'HD':
       return 'DI'
    elif material == 'FV':
       return 'GF'
    elif material == 'CAN':
       return 'AR'
    elif material == 'FOFO':
       return 'CI'
    else:
       return material

def plot_attribute_distribution(dataframe, field_name,title, label, figsize,log_scale):
    # Step 1: Extract necessary data
    dataframe['MATERIAL']=dataframe['MATERIAL'].apply(lambda x: replace_material_labels(x))
    material_types = dataframe['MATERIAL'].unique()
    
    my_order = my_order = dataframe.groupby(by=["MATERIAL"])[field_name].median().sort_values().index[::-1]
    
    # Step 2: Create violin plot
    plt.style.use('ggplot')
    plt.figure(figsize=figsize)
    sns.boxplot(y='MATERIAL', x=field_name
                , data=dataframe,order=my_order,palette='RdYlGn')
    
    plt.title(title,fontsize=28, fontweight='bold')
    if log_scale:
        plt.xscale('log')
    plt.xlabel(label,fontsize=22, fontweight='bold')
    plt.ylabel(label,fontsize=22, fontweight='bold')
    plt.yticks(fontsize=18, fontweight='bold')
    plt.xticks(fontsize=18, fontweight='bold')
    plt.grid(True, color='black', linestyle='-', linewidth=1.0)
    plt.grid(True, color='k', linestyle='-', linewidth=0.2,which='minor')
    
    #plt.xscale('log')
    plt.savefig(f'{title}.png', bbox_inches='tight')
    plt.show()
    print('Hi')
def main ():
    T=[1,10,50,100,500]
    
    
    parser = argparse.ArgumentParser(description='Process simulation parameters')
    parser.add_argument('--Model', type=str, help='Flag for indicating if NoCorr or model1 or model4 or model5')
    args = parser.parse_args()
    ModelName=args.Model
    #NoCorr
    fail_rate=f'Results/Results_{ModelName}'
    NPGV=f'Results/Results_N_PGV_{ModelName}'
    NPGD=f'Results/Results_N_PGD_{ModelName}'
    RRPGV=f'Results/Results_RR_PGV_{ModelName}'
    RRPGD=f'Results/Results_RR_PGD_{ModelName}'
    
    keys=[fail_rate
          ,NPGV
          ,NPGD
          ,RRPGD
          ,RRPGV
          ]

    
    dataFrames=[gpd.read_file(key_i) for key_i in keys]
    i=0
    
    for T_i in T:
        dataFrames[0][f'Pf_PGD_T_{T_i}']=1-np.exp(-1*np.array(dataFrames[0]['failRate_PGD'])*T_i)
        dataFrames[0][f'Pf_PGV_T_{T_i}']=1-np.exp(-1*np.array(dataFrames[0]['failRate_PGV'])*T_i)
        dataFrames[0][f'Pf_T_{T_i}']=1-np.exp(-1*(np.array(dataFrames[0]['failRate_PGD'])+np.array(dataFrames[0]['failRate_PGV']))*T_i)
    
    export_df=dataFrames[0]
    
    for i in enumerate(dataFrames[1:]):
        export_df=export_df.merge(i[1],on=['Global_site_id','geometry'])
    #export_df['MATERIAL']=export_df['MATERIAL'].apply(lambda x: replace_material_labels(x))
    field_Name='failRate_PGD'    
    Corr_Model='No Correlation'
    flag='Fail rate PGD No corr'
    title=f'{flag} by Material'
    label='Fail Rate [fail/yr]'
    plot_attribute_distribution(export_df,field_Name,title, label, figsize=(10,8),log_scale=False)
    #------------ fail rate PGV
    field_Name='failRate_PGV'    
    Corr_Model='GH 2008'
    flag='Fail rate PGV No corr'
    title=f'{flag} by Material'
    label='Fail Rate [fail/yr]'
    plot_attribute_distribution(export_df,field_Name,title, label, figsize=(10,8),log_scale=False)
    #-----------RR PGD
    field_Name='RR_PGD_50.0_yr'    
    Corr_Model='B 2003'
    flag='RR PGD No corr 50 yr'
    title=f'{flag} by Material'
    label='Repair Rate [fail/km]'
    plot_attribute_distribution(export_df,field_Name,title, label, figsize=(10,8),log_scale=False)
    #-----------RR PGV
    field_Name='RR_PGV_50.0_yr'    
    Corr_Model='GH 2008'
    flag='RR PGV No corr 50 yr'
    title=f'{flag} by Material'
    label='Repair Rate [fail/km]'
    plot_attribute_distribution(export_df,field_Name,title, label, figsize=(10,8),log_scale=False)
    #load 
if __name__ == '__main__':

    main() 
    #'--NoCorr'