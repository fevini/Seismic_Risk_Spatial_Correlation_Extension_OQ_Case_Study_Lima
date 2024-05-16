# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 12:25:30 2024

@author: nivfe
"""

import pandas as pd
import numpy as np
import argparse

def computeDisplacementCorrection(M):
    K=0.0086*np.power(M,3)-0.0914*np.power(M,2)+0.4698*M-0.9835
    return K
def susceptibiltyCat(Vs30):
    #equation in ft per second- convert from m/
    z=30.0
    vs=Vs30
    PGAt=(0.0001*np.power(vs,2))/(1.2*z)
    return PGAt
def latSpreadingDisplacement(PGA, PGAt):
    x = np.divide(PGA,PGAt)
    # Initialize an array for displacement with the same shape as x
    disp = np.zeros_like(x)
    
    # Condition 1: x <= 1
    # disp is already initialized to 0.0, so no action is required for this condition
    
    # Condition 2: 1 < x <= 2
    mask1 = (x > 1) & (x <= 2)
    disp[mask1] = (12 * x[mask1] - 12) * 2.54
    
    # Condition 3: 2 < x <= 3
    mask2 = (x > 2) & (x <= 3)
    disp[mask2] = (18 * x[mask2] - 24) * 2.54
    
    # Condition 4: 3 < x <= 4
    mask3 = (x > 3) & (x <= 4)
    disp[mask3] = (70 * x[mask3] - 180) * 2.54
    
    # Condition 5: x > 4
    mask4 = x > 4
    disp[mask4] = 254.0
    
    return disp
def computePGD(disp,K):
    
    return np.multiply(disp,K)
def test():
    PGA_ex=np.vstack(np.array([[0.15154],[0.151539]]))
    
    Vs_ex=np.vstack(np.array([[353.55],[565.387]]))
    
    M_ex=np.vstack(np.array([[7.05],[6.55]]))
    
    K_ex=computeDisplacementCorrection(M_ex)
    
    computeDisplacementCorrection(7.05)
    
    
    PGAt_ex=susceptibiltyCat(Vs_ex)
    
    susceptibiltyCat(353.55)
    
    
    disp_ex=latSpreadingDisplacement(PGA_ex, PGAt_ex)
    
    latSpreadingDisplacement(0.15154, 0.03472155625000001)
    
    PGD_ex=computePGD(disp_ex,K_ex)
    print(PGD_ex)


def main():
    #load IM
    parser = argparse.ArgumentParser(description='Process simulation parameters')
    parser.add_argument('--eventFileName', type=str, help='Path to the events file')
    parser.add_argument('--rupFileName', type=str, help='Path to the rup file')
    parser.add_argument('--IMFileName', type=str, help='Path to the IM file')
    parser.add_argument('--siteFileName', type=str, help='Path to the Vs30 file')
    parser.add_argument('--siteMesh', type=str, help='Path to the site Mesh coordinates')
    parser.add_argument('--outputFileName', type=str, help='Path to the Output file')
    
    args=parser.parse_args()
    
    eventsFileName=args.eventFileName#"events_66.csv"
    rupturesFileName=args.rupFileName#"ruptures_66.csv"
    IMFileName=args.IMFileName#"gmf-data_66.csv"
    #IMFileName=args.IMFileName
    events_df=pd.read_csv(eventsFileName,header=1)
    IM_df=pd.read_csv(IMFileName,header=0)
    ruptures_df=pd.read_csv(rupturesFileName,header=1)
    #Load Sites
    VS30FileName=args.siteFileName#'Vs30_Lima.csv'
    VS30_df=pd.read_csv(VS30FileName)
    
    meshFileName=args.siteMesh#"sitemesh_66.csv"
    mesh_df=pd.read_csv(meshFileName,header=1)
    VS30_df=VS30_df.merge(mesh_df)
    merged_df = pd.merge(events_df, ruptures_df, on='rup_id', how='left')
    
    merged_df = merged_df[['event_id','mag','rup_id']]
    IM_merged_df = pd.merge(IM_df, merged_df[['event_id', 'mag']],
                     on='event_id', how='left')
    
    IM_merged_df = pd.merge(IM_merged_df,VS30_df[['site_id','vs30']],
                            on='site_id', how='left')
    
    PGA=np.vstack(np.array(IM_merged_df['gmv_PGA']))
    
    Mag=np.vstack(np.array(IM_merged_df['mag']))
    
    VS30=np.vstack(np.array(IM_merged_df['vs30']))
    
    K=computeDisplacementCorrection(Mag)
    PGAt=susceptibiltyCat(VS30)
    disp=latSpreadingDisplacement(PGA,PGAt)
    PGD=computePGD(disp,K)
    #IM_df=pd.merge(IM_df,mesh_df,on='site_id')
    
    IM_merged_df['gmv_PGD']=PGD
    
    # #Load_PGV_file
    # PGVFileName='gmf-data_corr_PGV_66_3.csv'
    # pgv=pd.read_csv(PGVFileName,header=0)
    # IM_merged_df=pd.merge(IM_merged_df,pgv,on=['site_id','event_id'])
    
    output_FileName=args.outputFileName#'gmf-data_66_complete.csv'
    IM_merged_df[['site_id','event_id','gmv_PGA','gmv_PGV','gmv_PGD','mag','vs30']].to_csv(output_FileName,index=False)
    print('Hello')
if __name__=='__main__':
    main()
    #
    #--eventFileName events_66.csv --rupFileName ruptures_66.csv --IMFileName gmf-data_corr_PGA_66_3.csv --siteFileName Vs30_lima.csv --siteMesh sitemesh_66.csv --outputFileName gmf-data_corr_complete_66_3.csv
    #--eventFileName events_16.csv --rupFileName ruptures_16.csv --IMFileName gmf-data_corr_PGA_16_3.csv --siteFileName Vs30_lima.csv --siteMesh sitemesh_16.csv --outputFileName gmf-data_corr_complete_16_3.csv
    #--eventFileName events_66.csv --rupFileName ruptures_66.csv --IMFileName gmf_corr/gmf-data_corr_66_PGA_Boore_2003_PGVBoore_2003.csv --siteFileName Vs30_lima.csv --siteMesh sitemesh_66.csv --outputFileName gmf_corr/complete_gmf-data_corr_66_PGA_Boore_2003_PGVBoore_2003.csv
    '''
    Correlation models combinations
    gmf-data_corr_16_PGA_Boore_2003_PGV_Goda_Hong_PGV.csv     
    gmf-data_corr_16_PGA_Goda_Atk_2009_PGV_Goda_Hong_PGV.csv 
    gmf-data_corr_16_PGA_J_Baker_2009_PGV_Goda_Hong_PGV.csv
    
    gmf-data_corr_16_PGA_Boore_2003_PGV_Wang_Takada_2005.csv 
    gmf-data_corr_16_PGA_Goda_Atk_2009_PGV_Wang_Takada_2005.csv 
    gmf-data_corr_16_PGA_J_Baker_2009_PGV_Wang_Takada_2005.csv
    '''
    '''
    Input example 
    #python3.8 compute_PGD_corr.py --eventFileName events_16.csv --rupFileName ruptures_16.csv --IMFileName gmf_corr/gmf-data_corr_16_PGA_Boore_2003_PGV_Goda_Hong_PGV.csv --siteFileName Vs30_lima.csv --siteMesh sitemesh_16.csv --outputFileName gmf_corr/complete_gmf-data_corr_16_PGA_Boore_2003_PGV_Goda_Hong_PGV.csv
    #python3.8 compute_PGD_corr.py --eventFileName events_16.csv --rupFileName ruptures_16.csv --IMFileName gmf_corr/gmf-data_corr_16_PGA_Goda_Atk_2009_PGV_Wang_Takada_2005.csv --siteFileName Vs30_lima.csv --siteMesh sitemesh_16.csv --outputFileName gmf_corr/complete_gmf-data_corr_16_PGA_Goda_Atk_2009_PGV_Wang_Takada_2005.csv
    #python3.8 compute_PGD_corr.py --eventFileName events_16.csv --rupFileName ruptures_16.csv --IMFileName gmf_corr/gmf-data_corr_16_PGA_Boore_2003_PGV_Wang_Takada_2005.csv --siteFileName Vs30_lima.csv --siteMesh sitemesh_16.csv --outputFileName gmf_corr/complete_gmf-data_corr_16_PGA_Boore_2003_PGV_Wang_Takada_2005.csv
    #python3.8 compute_PGD_corr.py --eventFileName events_16.csv --rupFileName ruptures_16.csv --IMFileName gmf_corr/gmf-data_corr_16_PGA_J_Baker_2009_PGV_Goda_Hong_PGV.csv --siteFileName Vs30_lima.csv --siteMesh sitemesh_16.csv --outputFileName gmf_corr/complete_gmf-data_corr_16_PGA_J_Baker_2009_PGV_Goda_Hong_PGV.csv
    #python3.8 compute_PGD_corr.py --eventFileName events_16.csv --rupFileName ruptures_16.csv --IMFileName gmf_corr/gmf-data_corr_16_PGA_Goda_Atk_2009_PGV_Goda_Hong_PGV.csv --siteFileName Vs30_lima.csv --siteMesh sitemesh_16.csv --outputFileName gmf_corr/complete_gmf-data_corr_16_PGA_Goda_Atk_2009_PGV_Goda_Hong_PGV.csv
    #python3.8 compute_PGD_corr.py --eventFileName events_16.csv --rupFileName ruptures_16.csv --IMFileName gmf_corr/gmf-data_corr_16_PGA_J_Baker_2009_PGV_Wang_Takada_2005.csv --siteFileName Vs30_lima.csv --siteMesh sitemesh_16.csv --outputFileName gmf_corr/complete_gmf-data_corr_16_PGA_J_Baker_2009_PGV_Wang_Takada_2005.csv
    '''
    
    