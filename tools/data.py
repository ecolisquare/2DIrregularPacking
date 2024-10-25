from tools.geofunc import GeoFunc  # Import the GeoFunc module for geometric functions
import pandas as pd  # Import pandas for data manipulation and analysis
import json  # Import JSON for parsing JSON formatted strings
from shapely.geometry import Polygon  # Import the Polygon class from Shapely for geometric operations

def normData(poly, num):
    '''Normalize the vertices of the polygon by a scaling factor.'''
    for ver in poly:
        ver[0] = ver[0] * num  # Scale the x-coordinate of the vertex
        ver[1] = ver[1] * num  # Scale the y-coordinate of the vertex

def get_data(file_path, scale=1):
    '''Read polygon data from a CSV file, normalize it, and sort by area.'''
    df = pd.read_csv(file_path)  # Read the CSV file into a DataFrame
    polygons = []  # Initialize a list to hold the polygons

    # Read and normalize polygon data
    for i in range(df.shape[0]):  # Iterate over each row in the DataFrame
        for j in range(df['num'][i]):  # Iterate for the number of polygons specified in the 'num' column
            poly = json.loads(df['polygon'][i])  # Load the polygon data from JSON format
            GeoFunc.normData(poly, scale)  # Normalize the polygon using the specified scale
            polygons.append(poly)  # Append the normalized polygon to the list

    # Calculate the area of each polygon and store it with the polygon
    polygons_with_area = [(Polygon(poly).area, poly) for poly in polygons]

    # Sort the polygons by area in descending order
    polygons_with_area.sort(key=lambda x: x[0], reverse=True)

    # Extract the sorted polygons, ignoring the area
    sorted_polygons = [poly for _, poly in polygons_with_area]

    return sorted_polygons  # Return the sorted list of polygons
