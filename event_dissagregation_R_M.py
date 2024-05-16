# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 08:46:18 2024

@author: nivfe
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from haversine import haversine,haversine_vector, Unit
import argparse



def main():
    parser = argparse.ArgumentParser(description='Process simulation parameters')
    parser.add_argument('--simmulationID', type=str, help='simmulation ID')
    parser.add_argument('--simmulationTime', type=int, help='Time period for hazard curves')
    args = parser.parse_args()
    simID=args.simmulationID#'66'
    simTime=args.simmulationTime#1000
    rupturesFileName=f"ruptures_{simID}.csv"
    eventsFileName=f"events_{simID}.csv"
    events_df=pd.read_csv(eventsFileName,header=1)
    ruptures_df=pd.read_csv(rupturesFileName,header=1)
    site_lon=-76.915602
    site_lat=-12.014
    df=events_df[['event_id','rup_id']].merge(ruptures_df,on='rup_id')
    df['distance']=df.apply(lambda row:
                            haversine((site_lat,site_lon), (row['centroid_lat'],row['centroid_lon'])),axis=1)
    bins_r=np.arange(min(df['distance']),max(df['distance']),3)
    bins_mag=np.arange(min(df['mag']),max(df['mag']),0.1)
    bins=[bins_r,bins_mag]
        # Create Pivot Table
    pivot_table = np.histogram2d(df['distance'], df['mag'], bins=bins)[0]
    
    # Create 3D Column Map
    fig = plt.figure(figsize=(10, 8))
    
    ax = fig.add_subplot(111, projection='3d')
    
    xpos, ypos = np.meshgrid(bins_r[:-1], bins_mag[:-1], indexing='ij')
    mask = pivot_table > 0
    xpos = xpos[mask]
    ypos = ypos[mask]
    zpos = np.zeros_like(xpos)
    
    dx = 3
    dy = 0.1
    dz = pivot_table[mask].ravel()
    
    # Color by dz
    norm = plt.Normalize(dz.min(), dz.max())
    colors = plt.cm.Reds(norm(dz))
    
    # Create 3D Column Map
    
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average', color=colors)
    
    ax.set_xlabel('Distance [km]')
    ax.set_ylabel('Mw')
    ax.set_zlabel('')  # No z-label
    
    # Invert y-axis label
    ax.invert_yaxis()
    
    # Set ticks
    ax.xaxis.set_major_locator(plt.MultipleLocator(100))
    ax.xaxis.set_minor_locator(plt.MultipleLocator(50))
    ax.yaxis.set_major_locator(plt.MultipleLocator(1))
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    
    # Set tick parameters
    # Minor ticks next to major ticks for x and y axes
    ax.xaxis.set_minor_formatter(plt.FuncFormatter(lambda x, _: f'{x:,.0f}\n'))
    ax.yaxis.set_minor_formatter(plt.FuncFormatter(lambda y, _: f'{y:.1f}\n'))
    # Set font properties for minor tick labels
    ax.xaxis.set_tick_params(which='minor', labelsize=8)
    ax.yaxis.set_tick_params(which='minor', labelsize=8)
    
    total_events = len(df)
    ax.text2D(0.95, 0.95, f'Total events in {simTime} yr = {total_events}',
              transform=ax.transAxes, ha='right', va='top',
              bbox=dict(facecolor='white', alpha=0.8, edgecolor='black'))
    
    # Title
    plt.title(f'Event Disaggregation by Magnitude and Distance {simTime} yr', fontsize=16)
    
    # Color bar
    sm = plt.cm.ScalarMappable(cmap='Reds', norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label('Count', rotation=270, labelpad=15)
    
    # Show the plot
    plt.style.use('ggplot')
    # Remove numbers from Z-axis
    ax.zaxis.set_major_formatter(plt.NullFormatter())
    plt.savefig(f'event_disaggregation_{simID}_{simTime}_yr.png', dpi=300)  # Save with 300 DPI
    plt.show()
    # Reset Z-axis ticks to default
    ax.zaxis.set_major_formatter(plt.ScalarFormatter(useOffset=False))
    # Histogram for Distance
    plt.figure(figsize=(6, 4))
    plt.hist(np.array(df['distance']), bins=bins_r, color='blue', alpha=0.7)
    plt.xlabel('Distance [km]')
    plt.ylabel('Count')
    plt.title('Distance Histogram')
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(100))
    plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(50))
    
    plt.gca().tick_params(axis='x', which='minor', labelsize=8)
    plt.gca().tick_params(axis='y', which='minor', labelsize=8)
    
    # Add annotation with total number of events
    ax = plt.gca()
    plt.text(0.95, 0.95, f'Total events in {simTime} yr = {total_events}',
             transform=plt.gca().transAxes, ha='right', va='top',
             bbox=dict(facecolor='white', alpha=0.8, edgecolor='black'))
    
    # Save Distance Histogram
    plt.savefig(f'distance_histogram_{simID}_{simTime}_yr.png', dpi=300)
    
    # Histogram for Magnitude
    plt.figure(figsize=(6, 4))
    plt.hist(np.array(df['mag']), bins=bins_mag, color='green', alpha=0.7)
    plt.xlabel('Mw')
    plt.ylabel('Count')
    plt.title('Magnitude Histogram')
    plt.gca().xaxis.set_major_locator(plt.MultipleLocator(1))
    plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(0.5))
    
    plt.gca().tick_params(axis='x', which='minor', labelsize=8)
    plt.gca().tick_params(axis='y', which='minor', labelsize=8)
    
    # Add annotation with total number of events
    ax = plt.gca()
    plt.text(0.95, 0.95, f'Total events in {simTime} yr = {total_events}',
             transform=plt.gca().transAxes, ha='right', va='top',
             bbox=dict(facecolor='white', alpha=0.8, edgecolor='black'))
    
    # Save Magnitude Histogram
    plt.savefig(f'magnitude_histogram_{simID}_{simTime}_yr.png', dpi=300)
    
    
    
    # Show plots
    plt.show()


if __name__:
    main() 
    '''
    USSAGE EXAMPLE
    #--simmulationID 66 --simmulationTime 1000
    #python3.8 event_dissagregationR_M.py --simmulationID 16 --simmulationTime 10000
    #python3.8 event_dissagregationR_M.py --simmulationID 17 --simmulationTime 100000
    '''