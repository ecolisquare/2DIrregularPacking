from shapely.geometry import Polygon  # Import the Polygon class from Shapely for geometric operations
from shapely.affinity import translate, rotate, scale  # Import functions for geometric transformations
from scipy.spatial import ConvexHull  # Import ConvexHull for calculating the convex hull
import numpy as np  # Import NumPy for numerical operations

def find_minimal_translation(polygon1, polygon2, direction='horizontal', step=1.0, buffer_distance=15.0):
    """
    Find the minimal translation distance to ensure two polygons do not overlap.
    
    :param polygon1: The first polygon.
    :param polygon2: The second polygon.
    :param direction: The direction of translation ('horizontal' or 'vertical').
    :param step: The step size for translation.
    :param buffer_distance: The buffer distance to consider for non-overlapping.
    :return: The minimal translation distance and the translated polygon.
    """
    distance = 0  # Initialize the distance
    buffer_polygon1 = polygon1.buffer(buffer_distance)  # Create a buffered version of polygon1

    # Continue translating polygon2 until it no longer intersects with the buffered polygon1
    while not buffer_polygon1.intersects(polygon2):
        if direction == 'horizontal':
            polygon2 = translate(polygon2, xoff=-step)  # Translate horizontally
        elif direction == 'vertical':
            polygon2 = translate(polygon2, yoff=-step)  # Translate vertically
        distance += step  # Increment the distance

    return distance, polygon2  # Return the distance and the translated polygon

def compute_merge_convex_hull_area(polygons):
    """
    Calculate the area of the minimal convex hull that encloses all polygons.
    
    :param polygons: A list of polygons, where each polygon is a list of points.
    :return: The area of the minimal convex hull polygon.
    """
    all_points = []  # Initialize a list to hold all points from the polygons

    # Collect all vertices from the polygons
    for poly in polygons:
        all_points.extend(poly.exterior.coords)

    # Convert the list of points to a NumPy array
    all_points = np.array(all_points)

    # Check if there are enough points to compute a convex hull
    if len(all_points) < 3:
        raise ValueError("Cannot compute a convex hull with fewer than 3 points")

    # Calculate the convex hull
    hull = ConvexHull(all_points)
    hull_points = all_points[hull.vertices]  # Get the points of the convex hull

    # Calculate the bounding box for the convex hull
    minx, miny, maxx, maxy = Polygon(hull_points).bounds

    # Return the area of the convex hull plus a small buffer based on its dimensions
    return Polygon(hull_points).area + 0.25 * max((maxx - minx), (maxy - miny))

def align_polygons(polygon1, polygon2, direction='horizontal', buffer_distance=15.0):
    """
    Align two polygons either horizontally or vertically.
    
    :param polygon1: The first polygon.
    :param polygon2: The second polygon.
    :param direction: The alignment direction ('horizontal' or 'vertical').
    :param buffer_distance: The buffer distance to maintain between polygons.
    :return: The aligned polygons.
    """
    # Convert lists to Polygon objects if necessary
    if isinstance(polygon1, list):
        polygon1 = Polygon(polygon1)
    if isinstance(polygon2, list):
        polygon2 = Polygon(polygon2)

    # Translate polygon1 to the bottom left corner
    polygon1 = translate(polygon1, -Polygon(polygon1).bounds[0], -Polygon(polygon1).bounds[1])
    buffer_polygon1 = polygon1.buffer(buffer_distance)  # Create a buffer around polygon1

    # Align polygon2 next to polygon1 based on the specified direction
    if direction == 'horizontal':
        polygon2 = translate(polygon2, 
                              Polygon(buffer_polygon1).bounds[2] - Polygon(polygon2).bounds[0], 
                              Polygon(buffer_polygon1).bounds[1] - Polygon(polygon2).bounds[1])
    elif direction == 'vertical':
        polygon2 = translate(polygon2, 
                              Polygon(buffer_polygon1).bounds[0] - Polygon(polygon2).bounds[0], 
                              Polygon(buffer_polygon1).bounds[3] - Polygon(polygon2).bounds[1])

    return polygon1, polygon2  # Return the aligned polygons

def find_best_arrangement(poly1, poly2):
    """
    Arrange two polygons in various configurations to find the arrangement with the highest utilization.
    
    :param poly1: The first polygon object.
    :param poly2: The second polygon object.
    :return: The best arrangement and its corresponding area of the enclosing rectangle.
    """
    arrangements = []  # Initialize a list to hold different arrangements

    # Case 1: Vertical alignment, no rotation
    poly1, poly2 = align_polygons(poly1, poly2, 'vertical')
    distance, poly2_moved = find_minimal_translation(poly1, poly2, direction='vertical')
    arrangements.append(("Vertical, no rotation", compute_merge_convex_hull_area([poly1, poly2_moved]), [poly1, poly2_moved]))

    # Case 2: Vertical alignment, second polygon rotated 180 degrees
    poly2_rotated = rotate(poly2, 180)
    poly1, poly2 = align_polygons(poly1, poly2, 'vertical')
    distance, poly2_moved = find_minimal_translation(poly1, poly2_rotated, direction='vertical')
    arrangements.append(("Vertical, second rotated 180 degrees", compute_merge_convex_hull_area([poly1, poly2_moved]), [poly1, poly2_moved]))

    # Case 3: Horizontal alignment, no rotation
    poly1, poly2 = align_polygons(poly1, poly2, 'horizontal')
    distance, poly2_moved = find_minimal_translation(poly1, poly2, direction='horizontal')
    arrangements.append(("Horizontal, no rotation", compute_merge_convex_hull_area([poly1, poly2_moved]), [poly1, poly2_moved]))

    # Case 4: Horizontal alignment, second polygon rotated 180 degrees
    poly2_rotated = rotate(poly2, 180)
    poly1, poly2 = align_polygons(poly1, poly2, 'horizontal')
    distance, poly2_moved = find_minimal_translation(poly1, poly2_rotated, direction='horizontal')
    arrangements.append(("Horizontal, second rotated 180 degrees", compute_merge_convex_hull_area([poly1, poly2_moved]), [poly1, poly2_moved]))

    # Find the arrangement with the smallest area
    best_arrangement = min(arrangements, key=lambda x: x[1])
    
    return best_arrangement  # Return the best arrangement
