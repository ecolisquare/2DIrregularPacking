from scipy.spatial import ConvexHull
from shapely.geometry import Polygon
import numpy as np
from bottom_left_fill import BottomLeftFill
import pdb
from shapely.geometry import Polygon
import copy
from tools.geofunc import GeoFunc
from tools.double_pack import find_best_arrangement

def compute_merge_convex_hull(polygons):
    """
    Compute the convex hull that encloses all polygons.
    
    :param polygons: List of polygons, where each polygon is represented as a list of points.
    :return: A polygon object representing the convex hull.
    """
    # Merge all vertices from the polygons
    all_points = []
    for poly in polygons:
        all_points.extend(poly)
    
    # Convert the points to a numpy array
    all_points = np.array(all_points)
    
    # If there are fewer than 3 points, a convex hull cannot be computed
    if len(all_points) < 3:
        raise ValueError("Cannot compute a convex hull with fewer than 3 points.")
    
    # Compute the convex hull
    hull = ConvexHull(all_points)
    
    # Extract the vertices of the convex hull
    hull_points = all_points[hull.vertices]
    
    # Return a polygon object containing the convex hull vertices
    return Polygon(hull_points)

def compute_convex_hull(polygons):
    """
    Compute the convex hull for each individual polygon.
    
    :param polygons: List of polygons, where each polygon is represented as a list of points.
    :return: List of convex hulls, one for each polygon.
    """
    polygons_con = []
    
    # For each polygon, compute its convex hull
    for poly in polygons:
        hull = ConvexHull(poly)
        hull_points = [poly[vertex] for vertex in hull.vertices]
        polygons_con.append(hull_points)
    
    return polygons_con

def translate_polygon(poly, dx, dy):
    """
    Translate a polygon by a given distance in x and y directions.
    
    :param poly: Polygon represented as a list of points.
    :param dx: Translation distance along the x-axis.
    :param dy: Translation distance along the y-axis.
    :return: Translated polygon.
    """
    return [(x + dx, y + dy) for x, y in poly]

def make_groups(polygons, tolerance=0.01):
    """
    Group polygons with similar areas.
    
    :param polygons: List of polygons, where each polygon is represented as a list of points.
    :param tolerance: Tolerance for area differences, default is 0.01 (1%).
    :return: List of grouped polygons, where each group contains polygons with similar areas.
    """
    res = []
    n = len(polygons)
    
    # Track which polygons have already been grouped
    used = [False] * n
    
    # Loop over each polygon
    for i in range(n):
        if used[i]:
            continue
        group = []
        poly1 = polygons[i]
        
        # Compute the area of the first polygon in the group
        poly1_area = Polygon(poly1).area
        group.append(poly1)
        used[i] = True
        
        # Find other polygons with similar areas to group together
        for j in range(i + 1, n):
            if used[j]:
                continue
            poly2 = polygons[j]
            poly2_area = Polygon(poly2).area
            
            # Group polygons whose area difference is within the tolerance
            if abs(poly1_area - poly2_area) / min(poly1_area, poly2_area) <= tolerance:
                group.append(poly2)
                used[j] = True
        
        res.append(group)
    
    return res

def compute_convex_hull_with_rotation_and_translation(polygons):
    """
    Repeatedly pair similar polygons and compute the optimal packing until no more pairs can be made.
    
    :param polygons: List of polygons, each represented as a list of points.
    :return: Two lists - the resulting convex hulls after packing and the final polygon groups after packing.
    """
    # Group polygons with similar areas
    polygon_groups = make_groups(polygons)
    polygons_con_res = []
    polygon_groups_res = []
    
    # Process each group of polygons
    for polys in polygon_groups:
        
        # If the group contains only one polygon, add it directly to the results
        if len(polys) < 2:
            polygons_con_res.append(polys[0])
            polygon_groups_res.append([polys[0]])
            continue
        
        # Rotate the polygons to lower their center of gravity
        polys = [GeoFunc.rotateToLowestCenterOfGravity(poly) for poly in polys]
        polys_copy = copy.deepcopy(polys)
        
        # If the group has an odd number of polygons, add the last polygon directly to the result
        if (len(polys) % 2) == 1:
            polygons_con_res.append(polys[(len(polys) - 1)])
            polygon_groups_res.append([polys[(len(polys) - 1)]])
        
        # Pair and pack polygons two by two, finding the best arrangement
        for i in range(len(polys) // 2):
            best_hull = None
            best_polygons = None
            
            # Test both arrangements of the polygon pair and select the one with the best packing efficiency
            best_arrangement = min(
                [find_best_arrangement(polys[2 * i], polys[2 * i + 1]), 
                 find_best_arrangement(polys[2 * i + 1], polys[2 * i])],
                key=lambda x: x[1]
            )
            
            print('Best arrangement: ', best_arrangement[0])
            
            # Extract the packed polygons and compute the convex hull
            best_polygons = [list(best_arrangement[2][0].exterior.coords), list(best_arrangement[2][1].exterior.coords)]
            best_hull = compute_merge_convex_hull(best_polygons)
            
            # Add the resulting convex hull and packed polygons to the results
            polygons_con_res.append(best_hull)
            polygon_groups_res.append(best_polygons)
    
    return polygons_con_res, polygon_groups_res

# Example usage:
# import pandas as pd
# import json
# from tools.geofunc import GeoFunc

# def get_data(file_path, scale = 1):
#     df = pd.read_csv(file_path)
#     polygons = []
#     
#     # Read and normalize polygon data
#     for i in range(df.shape[0]):
#         for j in range(df['num'][i]):
#             poly = json.loads(df['polygon'][i])
#             GeoFunc.normData(poly, scale)
#             polygons.append(poly)
#     
#     # Compute the area of each polygon and store it along with the polygon
#     polygons_with_area = [(Polygon(poly).area, poly) for poly in polygons]
#     
#     # Sort polygons by area in descending order
#     polygons_with_area.sort(key=lambda x: x[0], reverse=True)
#     
#     # Extract the sorted polygons
#     sorted_polygons = [poly for _, poly in polygons_with_area]
#     return sorted_polygons

# polygons = get_data('/home/user8/QZ-iregularpacking/jndata/Test-Data/43.csv')

# polygons_con, polygons_pair_res = compute_convex_hull_with_rotation_and_translation(polygons)
# # pdb.set_trace()

# import matplotlib.pyplot as plt
# import matplotlib.patches as patches
# import numpy as np

# def plot_polygons(original_polygons, convex_hulls):
#     """
#     Plot the original polygons and their corresponding convex hulls.
#     
#     :param original_polygons: List of original polygons.
#     :param convex_hulls: List of convex hulls for the polygons.
#     """
#     for i, (polys, hull) in enumerate(zip(original_polygons, convex_hulls)):
#         # Ensure the hull is a NumPy array
#         hull = np.array(hull.exterior.coords)
#         
#         # Create a new plot
#         fig, ax = plt.subplots()
#         
#         # Set axis limits based on the convex hull size
#         mi = min(min(hull[:, 0]) - 10, min(hull[:, 1]) - 10)
#         ma = max(max(hull[:, 0]) + 10, max(hull[:, 1]) + 10)
#         ax.set_xlim(mi, ma)
#         ax.set_ylim(mi, ma)
#         
#         # Plot each polygon in the group
#         for poly in polys:
#             polygon = patches.Polygon(poly, closed=True, fill=True, edgecolor='black', facecolor='blue', alpha=0.5)
#             ax.add_patch(polygon)
#         
#         # Plot the convex hull
#         con_polygon = patches.Polygon(hull, closed=True, fill=False, edgecolor='green')
#         ax.add_patch(con_polygon)
#         
#         # Save the plot as an image file
#         plt.savefig(f'./e/e_{i}.png')
#         plt.close()

# # Call the plot function
# plot_polygons(polygons_pair_res, polygons_con)
