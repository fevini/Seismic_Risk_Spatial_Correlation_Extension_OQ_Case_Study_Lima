# -*- coding: utf-8 -*-
"""
Created on Sun Feb 11 08:45:52 2024

@author: nivfe
"""
import geopandas as gpd
import pandas as pd
import json
import xml.etree.ElementTree as ET

# Function to parse the XML file and extract k1 and k2 coefficients
def parse_xml(xml_file):
    coefficients = {}
    tree = ET.parse(xml_file)
    root = tree.getroot()
    for material in root.findall('material'):
        name = material.find('name').text
        diameter = material.find('diameter').text
        k1 = float(material.find('k1').text)
        k2 = float(material.find('k2').text)
        coefficients[(name, diameter)] = (k1, k2)
    return coefficients

# Function to extend the loaded JSON file with k1 and k2 attributes
def extend_json(json_data, coefficients):
    for obj in json_data:
        material = obj['MATERIAL']
        diameter = obj['DIAMETER']
        if diameter <= 300:
            diameter_size = 'small'
        else:
            diameter_size = 'large'
        k1, k2 = coefficients.get((material, diameter_size), (0.0, 0.0))
        obj['k1'] = k1
        obj['k2'] = k2
    return gpd.GeoDataFrame(json_data,crs="EPSG:4326")


# Export DataFrame to JSON file with given name
def export_dataframe_to_json(df, file_name):
    df.to_json(file_name, orient='records', indent=4)

# Example usage
xml_file = 'Constants.xml'
json_file = 'RED_PRIM_python_Modified.geojson'

# Parse XML file to extract coefficients
coefficients = parse_xml(xml_file)

# Load JSON file into Pandas DataFrame
df = gpd.read_file(json_file)
#df_geo=json.
# Extend DataFrame with k1 and k2 attributes
Geo_df=extend_json(df.to_dict('records'), coefficients)
file_name='extended_data.geojson'
# Export DataFrame to JSON file with a given name
Geo_df.to_file(file_name, driver='GeoJSON')

