import numpy as np
import pandas as pd
import hazard as hz
import correlation_3 as cor
import argparse
import matplotlib.pyplot as plt
import gc
from scipy.stats import norm, lognorm, genextreme
from scipy.optimize import curve_fit

"""
Load the data of the intensity meassures for all the events that are in the
gmf (ground motion fields) using the pandas
"""

  

'extend to create a dict that can run all the simmulations'

# Define command-line arguments
parser = argparse.ArgumentParser(description='Process simulation parameters')
#parser.add_argument('--IMFileName', type=str, help='Path to the IM file')
#parser.add_argument('--eventFileName', type=str, help='Path to the events file')
#parser.add_argument('--meshFileName', type=str, help='Path to the mesh file')
parser.add_argument('--simID',type=str,help='simulation ID')
# parser.add_argument('--simmulationTimeSpan', type=float, help='Simulation time span')
parser.add_argument('--corrModelPGA', type=str, help='Spatial correlation model used')
parser.add_argument('--corrModelPGV', type=str, help='Spatial correlation model used')

#parser.add_argument('--csvName',type=str,help='Path or Filename of the desired Export')
#parser.add_argument('--ruptures',type=str,help='Path or Filename of the ruptures File')
# # Parse the command-line arguments
args = parser.parse_args()

# # Check if required arguments are provided
# if args.IMFileName is None or args.meshFileName is None or args.simmulationTimeSpan is None or args.timePeriod is None or args.csvName is None or args.IMtype is None:
#     print("Error: Please provide all required arguments.")
#     parser.print_help()
#     exit()


'_______________'



def main():
    print('Hello')
    simID=args.simID
    #load IM
    IMFileName=f'gmf-data_{simID}.csv'
    #IMFileName=args.IMFileName
    IM_df=pd.read_csv(IMFileName,header=1)
    #Load Sites
    meshFileName=f'sitemesh_{simID}.csv'
    #meshFileName=args.meshFileName
    mesh_df=pd.read_csv(meshFileName,header=1)
    corrModelName=[args.corrModelPGA,args.corrModelPGV]
    IM_LabelArray=['PGA','PGV']
    #IM_Label='PGA'
    
    eventsFileName=f'events_{simID}.csv'#args.eventFileName
    events_df=pd.read_csv(eventsFileName,header=1)
    
    ruptures=f'ruptures_{simID}.csv'#args.ruptures
    #sig_eps="sigma_epsilon_66.csv"
    ruptures_df=pd.read_csv(ruptures,header=1)
    
    
    #path=f'gmf-data_corr_{IM_Label}_16_2.csv'
    event_df=events_df[['event_id','rup_id']].merge(ruptures_df[['rup_id','trt']],on='rup_id')
    #This is highly dependent of the GMM, and needs to be addressed
    hazard_data = {
        "Subduction IntraSlab": {
            "Inter_PGA": 0.266,
            "Intra_PGA": 0.101232,
            "Inter_PGV": 0.239,
            "Intra_PGV": 0.0953
        },
        "Subduction Interface": {
            "Inter_PGA": 0.4513066,
            "Intra_PGA": 0.6539,
            "Inter_PGV": 0.317756,
            "Intra_PGV": 0.449
        },
        "Active Shallow Crust": {
            "Inter_PGA": 0.363266667,
            "Intra_PGA": 0.631783333,
            "Inter_PGV": 0.37368,
            "Intra_PGV": 0.601667
        }
    }
    
    df_array=[]
    #merge data frames and group them by id and lon and lat
    i=0
    for IM_Label in IM_LabelArray:
        print(f'Merging Data Frames...{IM_Label}')
        IM_col=f'gmv_{IM_Label}'
        hazard_col=f'{IM_Label}_Hazard'
        IM_merged_df=pd.merge(IM_df,mesh_df,on='site_id')
        IM_merged_df.set_index(['site_id','event_id'],inplace=True)
        
        
        IM_col=f'gmv_{IM_Label}'
        IM_pivot_event_df=IM_merged_df.pivot_table(index='site_id',columns='event_id',values=[IM_col])
        IM_pivot_event_df = IM_pivot_event_df.fillna(0)
        
        
        
        
        event_df[f'Intra_{IM_Label}']=event_df.apply(lambda x:cor.get_intra(x['trt'], IM_Label, hazard_data) ,axis=1)
        
        print(f'applying spatial correlation...{IM_Label}')
        correlation_matrix =np.load(f'cor_mat_{corrModelName[i]}.npy')
        corr_df_IM=cor.modify_IM_with_correlation(IM_pivot_event_df,ruptures_df,event_df,correlation_matrix,IM_Label)
        df_array.append(corr_df_IM)
        i=i+1
    corr_df=pd.merge(df_array[0],df_array[1],on=['site_id','event_id'])
    path=f'gmf_corr/gmf-data_corr_{simID}_PGA_{corrModelName[0]}_PGV_{corrModelName[1]}.csv'
    corr_df.to_csv(path,index=False)
if __name__=='__main__':
    main()
#python3.8 main_corr_data_frame_3.py --simID 16 --corrModelPGA Boore_2003 --corrModelPGV Goda_Hong_PGV
#python3.8 main_corr_data_frame_3.py --simID 16 --corrModelPGA Boore_2003 --corrModelPGV Wang_Takada_2005
#python3.8 main_corr_data_frame_3.py --simID 16 --corrModelPGA Goda_Atk_2009 --corrModelPGV Wang_Takada_2005
#python3.8 main_corr_data_frame_3.py --simID 16 --corrModelPGA Goda_Atk_2009 --corrModelPGV Goda_Hong_PGV
#python3.8 main_corr_data_frame_3.py --simID 16 --corrModelPGA J_Baker_2009 --corrModelPGV Wang_Takada_2005
#python3.8 main_corr_data_frame_3.py --simID 16 --corrModelPGA J_Baker_2009 --corrModelPGV Goda_Hong_PGV