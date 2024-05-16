# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 01:27:54 2024

@author: nivfe
"""

#rutina para unir los resultados de el ana+àlisis de riesgo de las tuberias para mostrar 
import geopandas as gpd
import pandas as pd 

def mergeo_geo_df(geo,df,keys):
    return geo.merge(df[[keys]],on=[keys[0]])
    

fileNameGeoJson='extended_data.geojson'
fileName='Muestra_epf_corr.csv'
geo_df = gpd.read_file(fileNameGeoJson)
geo_df['Global_site_id']=geo_df.apply(
        lambda row: row['GLOBALID']+'_'+str(row['site_id']),axis=1)
epf_df=pd.read_csv(fileName,header=0)

keys=['Global_site_id','Epf_PGD','Epf_PGV']

merge_df=geo_df.merge(epf_df[keys],on=keys[0])

gpd.GeoDataFrame(merge_df[['GLOBALID','MATERIAL','DIAMETER','site_id','Epf_PGD','Epf_PGV','geometry']],
                           #geometry=merge_df[['geometry']],
                           crs="EPSG:4326").to_file('Example_ef_corr', driver='GeoJSON')

