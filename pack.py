from shapely.geometry import Polygon
import pandas as pd
import json
import datetime
import bottom_left_fill
import pdb
import csv
from tools.geofunc import GeoFunc
from tools.data import get_data
from tools.plot import plot_res
import copy
from arg import get_args
from group_merge import compute_convex_hull_with_rotation_and_translation
from scipy.spatial import ConvexHull

def compute_convex_hull(polygons):
    """
    Computes the convex hull for a list of polygons.
    
    :param polygons: List of polygons (each polygon is a list of points)
    :return: List of convex hulls, each represented by a list of points
    """
    polygons_con = []
    for poly in polygons:
        hull = ConvexHull(poly)
        hull_points = [poly[vertex] for vertex in hull.vertices]
        polygons_con.append(hull_points)
    return polygons_con

def compute_mbr(polygons):
    """
    Computes the minimum bounding rectangle (MBR) for each polygon.
    
    :param polygons: List of polygons (each polygon is a list of points)
    :return: List of MBRs, each represented by 4 points of the rectangle
    """
    polygons_mbr = []
    for poly in polygons:
        polygon = Polygon(poly)
        minx, miny, maxx, maxy = polygon.bounds
        polygons_mbr.append([[minx, miny], [maxx, miny], [maxx, maxy], [minx, maxy]])
    return polygons_mbr

def translate_polygon(polygon, dx, dy):
    """ 
    Translates a polygon by a given distance along the x and y axes.
    
    :param polygon: List of points representing the polygon
    :param dx: Distance to translate along the x-axis
    :param dy: Distance to translate along the y-axis
    :return: Translated polygon (list of points)
    """
    return [(x + dx, y + dy) for x, y in polygon]

def convert_polygons_to_list(polygons):
    """
    Converts all Polygon objects to list format.
    
    :param polygons: List of polygons (each polygon is a Polygon object or already a list)
    :return: List of polygons, all in list format
    """
    converted_polygons = []
    
    for poly in polygons:
        if isinstance(poly, Polygon):
            # Convert Polygon object to list
            poly_list = list(poly.exterior.coords)
            converted_polygons.append(poly_list)
        else:
            # If it's already a list, just append
            converted_polygons.append(poly)
    
    return converted_polygons

def main():
    args = get_args()  # Parse command line arguments
    polygons = get_data(args.input_path)  # Read input data (polygons)
    polygons_area = sum(Polygon(poly).area for poly in polygons)  # Calculate total area
    
    starttime = datetime.datetime.now()  # Start timing

    if args.grouping:
        # Perform grouping and lower the center of gravity for polygons
        polygons, polygon_groups = compute_convex_hull_with_rotation_and_translation(polygons)
        polygons = convert_polygons_to_list(polygons)
        for i, (polygon, polygon_group) in enumerate(zip(polygons, polygon_groups)):
            if len(polygon_group) < 2:
                polygons[i] = GeoFunc.rotateToLowestCenterOfGravity(polygon)
                polygon_groups[i][0] = GeoFunc.rotateToLowestCenterOfGravity(polygon)
        polygons_copy = copy.deepcopy(polygons)  # Backup original polygons
    else:
        # Lower the center of gravity for all polygons
        polygons = [GeoFunc.rotateToLowestCenterOfGravity(polygon) for polygon in polygons]
        polygons_copy = copy.deepcopy(polygons)  # Backup original polygons

    if args.packing_method == 'ConvexHull':
        # Compute convex hull for polygons
        polygons_con = compute_convex_hull(polygons)
        polygons_con_copy = copy.deepcopy(polygons_con)  # Backup convex hulls
    elif args.packing_method == 'MBR':
        # Compute minimum bounding rectangle (MBR) for polygons
        polygons_con = compute_mbr(polygons)
        polygons_con_copy = copy.deepcopy(polygons_con)  # Backup MBRs

    # Expand convex hulls with part gap
    polygons_con = [GeoFunc.expand_convex_hull(con, args.part_gap / 2) for con in polygons_con]

    # Start BLF (Bottom Left Fill) packing algorithm
    blf = bottom_left_fill.BottomLeftFill(args.container_width - 2 * (args.board_gap - args.part_gap / 2),
                                          polygons_con, vertical=False, NFPAssistant=None)

    endtime = datetime.datetime.now()  # End timing

    # Restore polygons to their final positions and synchronize the packing results
    container_width = args.container_width
    container_height = blf.contain_length + 2 * (args.board_gap - args.part_gap / 2)
    polygons_con_res = [translate_polygon(poly, args.board_gap - args.part_gap / 2, args.board_gap - args.part_gap / 2) 
                        for poly in blf.polygons]
    polygons_con_res = [GeoFunc.restore_convex_hull(poly, args.part_gap / 2) for poly in polygons_con_res]
    polygons_res = []

    if args.grouping:
    # Check if the packing method is based on Convex Hull
        if args.packing_method == 'ConvexHull':
            # Iterate over each polygon, its corresponding convex hull, convex hull result after packing, rotation, and its group of similar polygons
            for poly, con, con_res, rotation, poly_group in zip(polygons, polygons_con, polygons_con_res, blf.rotation, polygon_groups):
                # Rotate all polygons in the group by the calculated rotation angle
                poly_group = GeoFunc.rotateGroupPoly(poly, poly_group, rotation)
                # Rotate the current polygon by the same rotation angle
                poly = GeoFunc.rotatePoly(poly, rotation)
                # Recompute the convex hull of the rotated polygon
                con = [poly[vertex] for vertex in ConvexHull(poly).vertices]
                
                # Compute the minimum x and y coordinates of the current convex hull before packing
                min_con_x = min(point[0] for point in con)
                min_con_y = min(point[1] for point in con)
                # Compute the minimum x and y coordinates of the convex hull after packing
                min_con_res_x = min(point[0] for point in con_res)
                min_con_res_y = min(point[1] for point in con_res)
                
                # Calculate the translation distances along the x and y axes based on the difference
                # between the current convex hull and the packed convex hull positions
                dx = min_con_res_x - min_con_x
                dy = min_con_res_y - min_con_y
                
                # Translate all polygons in the group by the calculated distances (dx, dy) 
                # to adjust their positions according to the packing result
                for polygon in poly_group:
                    translated_poly = translate_polygon(polygon, dx, dy)
                    # Append the translated polygon to the final results list
                    polygons_res.append(translated_poly)
        
        # If the packing method is based on the Minimum Bounding Rectangle (MBR)
        elif args.packing_method == 'MBR':
            # Iterate over each polygon, its MBR, the MBR result after packing, rotation, and group of similar polygons
            for poly, mbr, mbr_res, rotation, poly_group in zip(polygons, polygons_con, polygons_con_res, blf.rotation, polygon_groups):
                # Rotate all polygons in the group by the calculated rotation angle
                poly_group = GeoFunc.rotateGroupPoly(poly, poly_group, rotation)
                # Rotate the current polygon by the same rotation angle
                poly = GeoFunc.rotatePoly(poly, rotation)
                
                # Get the bounds (min x, min y, max x, max y) of the rotated polygon
                minx, miny, maxx, maxy = Polygon(poly).bounds
                # Construct the MBR using the bounds of the rotated polygon
                mbr = [[minx, miny], [maxx, miny], [maxx, maxy], [minx, maxy]]
                
                # Compute translation distances by comparing the packed MBR and original MBR
                dx = min(point[0] for point in mbr_res) - min(point[0] for point in mbr)
                dy = min(point[1] for point in mbr_res) - min(point[1] for point in mbr)
                
                # Translate all polygons in the group by the calculated distances (dx, dy)
                for polygon in poly_group:
                    translated_poly = translate_polygon(polygon, dx, dy)
                    # Append the translated polygon to the final results list
                    polygons_res.append(translated_poly)

    else:
        # If there is no grouping, handle each polygon individually
        if args.packing_method == 'ConvexHull':
            # Iterate over each polygon, its convex hull, convex hull result after packing, and rotation
            for poly, con, con_res, rotation in zip(polygons, polygons_con, polygons_con_res, blf.rotation):
                # Rotate the polygon by the calculated rotation angle
                poly = GeoFunc.rotatePoly(poly, rotation)
                # Recompute the convex hull of the rotated polygon
                con = [poly[vertex] for vertex in ConvexHull(poly).vertices]
                
                # Compute the minimum x and y coordinates of the current convex hull
                min_con_x = min(point[0] for point in con)
                min_con_y = min(point[1] for point in con)
                # Compute the minimum x and y coordinates of the convex hull after packing
                min_con_res_x = min(point[0] for point in con_res)
                min_con_res_y = min(point[1] for point in con_res)
                
                # Calculate the translation distances (dx, dy) to match the positions of the packed result
                dx = min_con_res_x - min_con_x
                dy = min_con_res_y - min_con_y
                
                # Translate the polygon by the calculated distances
                translated_poly = translate_polygon(poly, dx, dy)
                # Append the translated polygon to the final results list
                polygons_res.append(translated_poly)
        
        # If the packing method is based on Minimum Bounding Rectangle (MBR)
        elif args.packing_method == 'MBR':
            # Iterate over each polygon, its MBR, the MBR result after packing, and rotation
            for poly, mbr, mbr_res, rotation in zip(polygons, polygons_con, polygons_con_res, blf.rotation):
                # Rotate the polygon by the calculated rotation angle
                poly = GeoFunc.rotatePoly(poly, rotation)
                
                # Get the bounds (min x, min y, max x, max y) of the rotated polygon
                minx, miny, maxx, maxy = Polygon(poly).bounds
                # Construct the MBR using the bounds of the rotated polygon
                mbr = [[minx, miny], [maxx, miny], [maxx, maxy], [minx, maxy]]
                
                # Compute the translation distances (dx, dy) by comparing the packed and original MBRs
                dx = min(point[0] for point in mbr_res) - min(point[0] for point in mbr)
                dy = min(point[1] for point in mbr_res) - min(point[1] for point in mbr)
                
                # Translate the polygon by the calculated distances
                translated_poly = translate_polygon(poly, dx, dy)
                # Append the translated polygon to the final results list
                polygons_res.append(translated_poly)

    
    # Print placement result info
    placement_time = endtime - starttime
    print("Placement time:", placement_time)
    print("Container height:", container_height)
    print("Container width:", container_width)
    space_utilization = polygons_area / (container_width * container_height)
    print("Space utilization:", space_utilization)

    # Output placement result image
    plot_res(polygons_res, polygons_con_res, container_width, container_height, args.output_path)

    # Output placement result data to CSV
    with open(args.output_path + 'packing_result.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Polygon Index', 'Coordinates'])
        for idx, poly in enumerate(polygons_res):
            writer.writerow([idx, poly])

if __name__ == "__main__":
    main()
