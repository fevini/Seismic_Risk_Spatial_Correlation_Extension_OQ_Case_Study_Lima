import pandas as pd
import numpy as np
from scipy.stats import norm, lognorm, genextreme
from scipy.spatial.distance import pdist, squareform
from haversine import haversine

def numba_multivariate_normal(mean, cov, size=1):
    dim = len(mean)
    L = np.linalg.cholesky(cov)
    samples = np.empty((size, dim))
    for i in range(size):
        z = np.random.normal(size=dim)
        samples[i] = mean + L @ z
    return samples
def distances(data):
    distances = pdist(data[['lat', 'lon']], lambda u, v: haversine(u, v))
    distance_matrix = squareform(distances)
    return distance_matrix
def compute_correlation_matrix_Vanmarcke(distance_matrix):
    
    correlation_matrix =2* np.exp(-(0.3 * distance_matrix)**0.5)-np.exp(-2*(0.3 * distance_matrix)**0.5)
    return correlation_matrix

def compute_correlation_matrix_Boore_2003(distance_matrix):
    c=1.0
    b=0.5
    a=0.774
    correlation_matrix =(1-c)+c*np.exp(-a*np.power(distance_matrix,b))
    return correlation_matrix
def compute_correlation_matrix_Goda_Atkinson_2009_PGA(distance_matrix):
    c=2.6
    b=0.336
    a=0.095
    correlation_matrix =(1-c)+c*np.exp(-a*np.power(distance_matrix,b))
    correlation_matrix[correlation_matrix<0]=0.0
    return correlation_matrix
def compute_correlation_matrix_Jayaram_Baker_2009_PGA(distance_matrix):
    c=1.0
    b=1.0
    a=0.333
    correlation_matrix =(1-c)+c*np.exp(-a*np.power(distance_matrix,b))
    correlation_matrix[correlation_matrix<0]=0.0
    return correlation_matrix
def compute_correlation_matrix_Goda_Hong_2008_PGV(distance_matrix):
    c=1.0
    b=0.5
    a=0.58
    correlation_matrix =(1-c)+c*np.exp(-a*np.power(distance_matrix,b))
    correlation_matrix[correlation_matrix<0]=0.0
    return correlation_matrix
def compute_correlation_matrix_Wang_Takada_2005_PGV(distance_matrix):
    c=1.0
    b=1.0
    a=0.048
    correlation_matrix =(1-c)+c*np.exp(-a*np.power(distance_matrix,b))
    correlation_matrix[correlation_matrix<0]=0.0
    return correlation_matrix

def extract_inter_event_residual(data,IM_type):
    #Extract the intervent residuals loaded in the data frame of sigma_epsilon
        sigma_col=f'sig_inter_{IM_type}'
        eps_col=f'eps_inter_{IM_type}'
        sigma=data[sigma_col]
        epsilon=data[eps_col]
        Inter_residual=np.multiply(sigma,epsilon)
        data[f'Inter_res_{IM_type}']=Inter_residual
        return data
        
#function extractiong the list of all the events
def extract_all_events(df):
    # Use a set to store unique events
    unique_events = set()

    # Iterate through the dataframe to extract events
    for event_list in df['event_id']:
        unique_events.update(event_list)

    # Convert the set to a sorted list
    all_events = sorted(list(unique_events))
    return all_events
#function computing the median of the log of the IM
def compute_median_log_im_vector(data,IM_type):
    # Initialize a vector with zeros
    median_log_im_vector = np.zeros(len(data))
    #SET THE iM FLAG
    IM=f'gmv_{IM_type}'
    # Iterate through the dataframe to compute the median of the log of IM for each site
    for index, row in data.iterrows():
        im_list = row[IM]
        site_median_log_im = np.median(np.log(im_list))
        median_log_im_vector[row['site_id'] - 1] = site_median_log_im

    return median_log_im_vector
def get_intra(hazard_type, intensity_measure,hazard_data):
    """
    Returns the inter and intra values for a given hazard type and intensity measure.

    Parameters:
    - hazard_type (str): The type of the hazard (e.g., "Subduction IntraSlab").
    - intensity_measure (str): The intensity measure (e.g., "PGA" or "PGV").

    Returns:
    - tuple: (inter_value, intra_value)
    """
    # Construct the keys for inter and intra based on the intensity measure
    
    intra_key = f"Intra_{intensity_measure}"
    
    # Extract the values
    
    intra_value = hazard_data[hazard_type].get(intra_key, "Key not found")
    
    return intra_value

def modify_IM_with_correlation(data,data_inter,event_df,correlation_matrix,IM_type):
    # List of the event ids in the dataframe
    print('aggregating events')
    all_events_list = data.columns.get_level_values('event_id').unique().to_numpy()
    #Compute medianlogIM
    print(f'aggregating median of ln{IM_type}')
    logdata=data.transform(lambda x:np.log(x))
    logdata.replace(-np.inf, -10000, inplace=True)
    
        #Compute the inter event residual 
    print('Computing Inter_residual')
    mean=np.zeros(2350)
    
    
    Intra_mat=np.random.multivariate_normal(mean, correlation_matrix, size=len(all_events_list))
    
    
    
    event_Intra=np.array(event_df[f'Intra_{IM_type}'][all_events_list])
    Intra_res=Intra_mat.T*event_Intra
    new_IM=np.exp(logdata+Intra_res)
    return_data=new_IM.stack(level='event_id').reset_index()
    return return_data
            
#__________________________________________________________



