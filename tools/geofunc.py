from shapely.geometry import Polygon, Point, mapping, LineString
from shapely.ops import unary_union
from shapely import affinity
import pyclipper
import math
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt

bias = 0.00001  # Calculation precision bias

class GeoFunc(object):
    """
    Geometry-related functions.
    1. checkBottom, checkTop, checkLeft, checkRight do not consider multiple points for now.
    2. Both checkBottom and checkLeft consider the bottom-left corner.
    """

    def almostContain(line, point):
        # Due to int, there will be calculation bias!!!!!!!!!!
        pt1 = [line[0][0], line[0][1]]
        pt2 = [line[1][0], line[1][1]]
        point = [point[0], point[1]]

        # Horizontal line case: compare the two points and the middle point
        if abs(pt1[1] - point[1]) < bias and abs(pt2[1] - point[1]) < bias:
            # print("Horizontal case")
            if (pt1[0] - point[0]) * (pt2[0] - point[0]) < 0:
                return True
            else:
                return False

        # Exclude the vertical case
        if abs(pt1[0] - point[0]) < bias and abs(pt2[0] - point[0]) < bias:
            # print("Vertical case")
            if (pt1[1] - point[1]) * (pt2[1] - point[1]) < 0:
                return True
            else:
                return False

        if abs(pt1[0] - point[0]) < bias or abs(pt2[0] - point[0]) < bias or abs(pt1[0] - pt2[0]) < bias:
            return False

        # Normal case, calculate the difference in arc
        arc1 = np.arctan((line[0][1] - line[1][1]) / (line[0][0] - line[1][0]))
        arc2 = np.arctan((point[1] - line[1][1]) / (point[0] - line[1][0]))
        if abs(arc1 - arc2) < bias:  # Original value 0.03, dighe approximate parallel correction to 0.01
            if (point[1] - pt1[1]) * (pt2[1] - point[1]) > 0 and (point[0] - pt1[0]) * (pt2[0] - point[0]) > 0:
                # print("General case")
                return True
            else:
                return False
        else:
            return False

    def computeInterArea(orginal_inter):
        """
        Calculate the area of the intersection area.
        """
        inter = mapping(orginal_inter)
        # A polygon
        if inter["type"] == "Polygon":
            if len(inter["coordinates"]) > 0:
                poly = inter["coordinates"][0]
                return Polygon(poly).area
            else: return 0
        if inter["type"] == "MultiPolygon":
            area = 0
            for _arr in inter["coordinates"]:
                poly = _arr[0]
                area = area + Polygon(poly).area
            return area

        if inter["type"] == "GeometryCollection":
            area = 0
            for _arr in inter["geometries"]:
                if _arr["type"] == "Polygon":
                    poly = _arr["coordinates"][0]
                    area = area + Polygon(poly).area
            return area
        return 0

    def checkBottom(poly):
        polyP = Polygon(poly)
        min_y = polyP.bounds[1]
        for index, point in enumerate(poly):
            if point[1] == min_y:
                return index

    def checkTop(poly):
        polyP = Polygon(poly)
        max_y = polyP.bounds[3]
        for index, point in enumerate(poly):
            if point[1] == max_y:
                return index

    def checkLeft(poly):
        polyP = Polygon(poly)
        min_x = polyP.bounds[0]
        for index, point in enumerate(poly):
            if point[0] == min_x:
                return index

    def checkRight(poly):
        polyP = Polygon(poly)
        max_x = polyP.bounds[2]
        for index, point in enumerate(poly):
            if point[0] == max_x:
                return index

    def checkBound(poly):
        return GeoFunc.checkLeft(poly), GeoFunc.checkBottom(poly), GeoFunc.checkRight(poly), GeoFunc.checkTop(poly)

    def checkBoundPt(poly):
        """Get the boundary points"""
        left, bottom, right, top = poly[0], poly[0], poly[0], poly[0]
        for i, pt in enumerate(poly):
            if pt[0] < left[0]:
                left = pt
            if pt[0] > right[0]:
                right = pt
            if pt[1] > top[1]:
                top = pt
            if pt[1] < bottom[1]:
                bottom = pt
        return left, bottom, right, top

    def checkBoundValue(poly):
        """Get the boundary values"""
        left, bottom, right, top = poly[0][0], poly[0][1], poly[0][0], poly[0][1]
        for i, pt in enumerate(poly):
            if pt[0] < left:
                left = pt[0]
            if pt[0] > right:
                right = pt[0]
            if pt[1] > top:
                top = pt[1]
            if pt[1] < bottom:
                bottom = pt[1]
        return left, bottom, right, top

    def slideToPoint(poly, pt1, pt2):
        GeoFunc.slidePoly(poly, pt2[0] - pt1[0], pt2[1] - pt1[1])

    def getSlide(poly, x, y):
        """
        Get the situation after translation.
        """
        new_vertex = []
        for point in poly:
            new_point = [point[0] + x, point[1] + y]
            new_vertex.append(new_point)
        return new_vertex

    @staticmethod
    def boundsContain(bounds, pt):
        if pt[0] >= bounds[0] and pt[0] <= bounds[2] and pt[1] >= bounds[1] and pt[1] <= bounds[3]:
            return True
        return False

    def slidePoly(poly, x, y):
        for point in poly:
            point[0] = point[0] + x
            point[1] = point[1] + y

    def polyToArr(inter):
        res = mapping(inter)
        _arr = []
        if res["type"] == "MultiPolygon":
            for poly in res["coordinates"]:
                for point in poly[0]:
                    _arr.append([point[0], point[1]])
        elif res["type"] == "GeometryCollection":
            for item in res["geometries"]:
                if item["type"] == "Polygon":
                    for point in item["coordinates"][0]:
                        _arr.append([point[0], point[1]])
        else:
            if res["coordinates"][0][0] == res["coordinates"][0][-1]:
                for point in res["coordinates"][0][0:-1]:
                    _arr.append([point[0], point[1]])
            else:
                for point in res["coordinates"][0]:
                    _arr.append([point[0], point[1]])
        return _arr

    def normData(poly, num):
        for ver in poly:
            ver[0] = ver[0] * num
            ver[1] = ver[1] * num

    """
    Approximate calculation
    """
    def crossProduct(vec1, vec2):
        res = vec1[0] * vec2[1] - vec1[1] * vec2[0]
        # Simplest calculation
        if abs(res) < bias:
            return 0
        # Some cases have a large cross product but are still basically parallel
        if abs(vec1[0]) > bias and abs(vec2[0]) > bias:
            if abs(vec1[1] / vec1[0] - vec2[1] / vec2[0]) < bias:
                return 0
        return res

    '''Used for touching calculation of intersection points, can be merged with another intersection point calculation function'''
    def intersection(line1, line2):
        # If an intersection point can be directly calculated
        Line1 = LineString(line1)
        Line2 = LineString(line2)
        inter = Line1.intersection(Line2)
        if inter.is_empty == False:
            mapping_inter = mapping(inter)
            if mapping_inter["type"] == "LineString":
                inter_coor = mapping_inter["coordinates"][0]
            else:
                inter_coor = mapping_inter["coordinates"]
            return inter_coor

        # Check if all vertices are the same
        res = []
        for pt1 in line1:
            for pt2 in line2:
                if GeoFunc.almostEqual(pt1, pt2) == True:
                    # print("pt1,pt2:", pt1, pt2)
                    res = pt1
        if res != []:
            return res

        # Calculate if there is an almostContain
        for pt in line1:
            if GeoFunc.almostContain(line2, pt) == True:
                return pt
        for pt in line2:
            if GeoFunc.almostContain(line1, pt) == True:
                return pt
        return []

    '''Mainly used to determine if there is a coincidence of straight lines, too complex and needs to be refactored'''
    def newLineInter(line1, line2):
        vec1 = GeoFunc.lineToVec(line1)
        vec2 = GeoFunc.lineToVec(line2)
        vec12_product = GeoFunc.crossProduct(vec1, vec2)
        Line1 = LineString(line1)
        Line2 = LineString(line2)
        inter = {
            "length": 0,
            "geom_type": None
        }
        # Only parallel lines can overlap
        if vec12_product == 0:
            # Copy to avoid affecting the original value
            new_line1 = GeoFunc.copyPoly(line1)
            new_line2 = GeoFunc.copyPoly(line2)
            if vec1[0] * vec2[0] < 0 or vec1[1] * vec2[1] < 0:
                new_line2 = GeoFunc.reverseLine(new_line2)
            # If there are equal vertices, choose one of them
            if GeoFunc.almostEqual(new_line1[0], new_line2[0]) or GeoFunc.almostEqual(new_line1[1], new_line2[1]):
                inter["length"] = min(Line1.length, Line2.length)
                inter["geom_type"] = 'LineString'
                return inter
            # Exclude the case where only vertices intersect
            if GeoFunc.almostEqual(new_line1[0], new_line2[1]):
                inter["length"] = new_line2[1]
                inter["geom_type"] = 'Point'
                return inter
            if GeoFunc.almostEqual(new_line1[1], new_line2[0]):
                inter["length"] = new_line1[1]
                inter["geom_type"] = 'Point'
                return inter
            # Otherwise, determine if one line is contained within the other
            line1_contain_line2_pt0 = GeoFunc.almostContain(new_line1, new_line2[0])
            line1_contain_line2_pt1 = GeoFunc.almostContain(new_line1, new_line2[1])
            line2_contain_line1_pt0 = GeoFunc.almostContain(new_line2, new_line1[0])
            line2_contain_line1_pt1 = GeoFunc.almostContain(new_line2, new_line1[1])
            # Line1 directly contains Line2
            if line1_contain_line2_pt0 and line1_contain_line2_pt1:
                inter["length"] = Line1.length
                inter["geom_type"] = 'LineString'
                return inter
            # Line2 directly contains Line1
            if line1_contain_line2_pt0 and line1_contain_line2_pt1:
                inter["length"] = Line2.length
                inter["geom_type"] = 'LineString'
                return inter
            # Mutually containing intersection points
            if line1_contain_line2_pt0 and line2_contain_line1_pt1:
                inter["length"] = LineString([line2[0], line1[1]]).length
                inter["geom_type"] = 'LineString'
                return inter
            if line1_contain_line2_pt1 and line2_contain_line1_pt0:
                inter["length"] = LineString([line2[1], line1[0]]).length
                inter["geom_type"] = 'LineString'
                return inter
        return inter

    def reverseLine(line):
        pt0 = line[0]
        pt1 = line[1]
        return [[pt1[0], pt1[1]], [pt0[0], pt0[1]]]

    '''Approximate calculation'''
    def almostEqual(point1, point2):
        if abs(point1[0] - point2[0]) < bias and abs(point1[1] - point2[1]) < bias:
            return True
        else:
            return False

    def extendLine(line):
        '''
        Extend a line
        '''
        pt0 = line[0]
        pt1 = line[1]
        vect01 = [pt1[0] - pt0[0], pt1[1] - pt0[1]]
        vect10 = [-vect01[0], -vect01[1]]
        multi = 40
        new_pt1 = [pt0[0] + vect01[0] * multi, pt0[1] + vect01[1] * multi]
        new_pt0 = [pt1[0] + vect10[0] * multi, pt1[1] + vect10[1] * multi]
        return [new_pt0, new_pt1]

    def getArc(line):
        if abs(line[0][0] - line[1][0]) < 0.01:  # Vertical case
            if line[0][1] - line[1][1] > 0:
                return 0.5 * math.pi
            else:
                return -0.5 * math.pi
        k = (line[0][1] - line[1][1]) / (line[0][0] - line[1][0])
        arc = np.arctan(k)
        return arc

    def extendInter(line1, line2):
        '''
        Get the intersection point of the extended lines
        '''
        line1_extend = GeoFunc.extendLine(line1)
        line2_extend = GeoFunc.extendLine(line2)
        # Check for parallel case
        k1 = GeoFunc.getArc(line1_extend)
        k2 = GeoFunc.getArc(line2_extend)
        if abs(k1 - k2) < 0.01:
            return [line1[1][0], line1[1][1]]
        inter = mapping(LineString(line1_extend).intersection(LineString(line2_extend)))
        if inter["type"] == "GeometryCollection" or inter["type"] == "LineString":
            return [line1[1][0], line1[1][1]]
        return [inter["coordinates"][0], inter["coordinates"][1]]

    def twoDec(poly):
        for pt in poly:
            pt[0] = round(pt[0], 2)
            pt[1] = round(pt[1], 2)

    def similarPoly(poly):
        '''
        Solve the approximate polygon of a convex polygon, the concave part of the concave polygon needs additional processing
        '''
        change_len = 10
        extend_poly = poly + poly
        Poly = Polygon(poly)
        new_edges = []
        # Calculate the translation of the straight line
        for i in range(len(poly)):
            line = [extend_poly[i], extend_poly[i + 1]]
            new_line = GeoFunc.slideOutLine(line, Poly, change_len)
            new_edges.append(new_line)

        # Calculate the extension of the straight line
        new_poly = []
        new_edges.append(new_edges[0])
        for i in range(len(new_edges) - 1):
            inter = GeoFunc.extendInter(new_edges[i], new_edges[i + 1])
            new_poly.append(inter)

        GeoFunc.twoDec(new_poly)

        return new_poly

    def slideOutLine(line, Poly, change_len):
        '''
        Translates a line outward by a specified length while checking if the new line intersects a given polygon.
        
        Parameters:
            line: A list containing two points defining the line (e.g., [[x1, y1], [x2, y2]]).
            Poly: A Polygon object to check against.
            change_len: The length by which to translate the line outward.
        
        Returns:
            new_line: A new line that has been translated outward.
        '''
        pt0 = line[0]  # Starting point of the line
        pt1 = line[1]  # Ending point of the line
        mid = [(pt0[0] + pt1[0]) / 2, (pt0[1] + pt1[1]) / 2]  # Midpoint of the line segment
        
        # Check if the line is vertical
        if pt0[1] != pt1[1]:
            # Calculate the slope's negative reciprocal for the perpendicular line
            k = -(pt0[0] - pt1[0]) / (pt0[1] - pt1[1])
            theta = math.atan(k)  # Angle of the line
            delta_x = 1 * math.cos(theta)  # Change in x
            delta_y = 1 * math.sin(theta)  # Change in y
            
            # Check if the new point is inside the polygon
            if Poly.contains(Point([mid[0] + delta_x, mid[1] + delta_y])):
                delta_x = -delta_x  # Invert direction if inside
                delta_y = -delta_y
            
            # Create the new line after translation
            new_line = [[pt0[0] + change_len * delta_x, pt0[1] + change_len * delta_y],
                        [pt1[0] + change_len * delta_x, pt1[1] + change_len * delta_y]]
            return new_line  # Return the translated line
        else:
            delta_y = 1  # Default change in y for horizontal lines
            if Poly.contains(Point([mid[0], mid[1] + delta_y])):
                delta_y = -delta_y  # Invert if inside
            
            # Create the new line for a horizontal case
            return [[pt0[0], pt0[1] + change_len * delta_y],
                    [pt1[0], pt1[1] + change_len * delta_y]]

    def copyPoly(poly):
        '''
        Creates a deep copy of a polygon represented as a list of points.
        
        Parameters:
            poly: A list of points defining the polygon.
        
        Returns:
            new_poly: A new list of points that is a copy of the original polygon.
        '''
        new_poly = []  # Initialize an empty list for the new polygon
        for pt in poly:
            new_poly.append([pt[0], pt[1]])  # Append each point to the new polygon
        return new_poly  # Return the copied polygon

    def pointLineDistance(point, line):
        '''
        Calculates the shortest distance from a point to a line defined by two points.
        
        Parameters:
            point: A list containing the coordinates of the point (e.g., [x, y]).
            line: A list containing two points that define the line (e.g., [[x1, y1], [x2, y2]]).
        
        Returns:
            distance: The shortest distance from the point to the line.
            direction: A vector indicating the direction from the point to the nearest point on the line.
        '''
        point_x = point[0]  # x-coordinate of the point
        point_y = point[1]  # y-coordinate of the point
        line_s_x = line[0][0]  # Starting point x-coordinate of the line
        line_s_y = line[0][1]  # Starting point y-coordinate of the line
        line_e_x = line[1][0]  # Ending point x-coordinate of the line
        line_e_y = line[1][1]  # Ending point y-coordinate of the line
        
        # Handle vertical line case
        if line_e_x - line_s_x == 0:
            return abs(point_x - line_s_x), [line_s_x - point_x, 0]  # Return distance and direction
        
        # Handle horizontal line case
        if line_e_y - line_s_y == 0:
            return abs(point_y - line_s_y), [0, line_s_y - point_y]  # Return distance and direction

        # Calculate the slope of the line
        k = (line_e_y - line_s_y) / (line_e_x - line_s_x)
        # Extend the line to find the intersection point
        extend_line = [[point_x - 1000, point_y - 1000 * (-1 / k)],
                    [point_x + 1000, point_y + 1000 * (-1 / k)]]
        
        # Calculate the intersection between the extended line and the given line
        inter = LineString(line).intersection(LineString(extend_line))
        
        if inter.is_empty:
            # Calculate the distance to both endpoints of the line
            dis1 = math.sqrt((point_x - line_s_x) ** 2 + (point_y - line_s_y) ** 2)
            dis2 = math.sqrt((point_x - line_e_x) ** 2 + (point_y - line_e_y) ** 2)
            # Return the smaller distance and corresponding direction vector
            if dis1 > dis2:
                return dis2, [line_e_x - point_x, line_e_y - point_y]
            else:
                return dis1, [line_s_x - point_x, line_s_y - point_y]
        else:
            pt = GeoFunc.getPt(inter)  # Get coordinates of the intersection point
            # Calculate distance to the intersection point
            dis = math.sqrt((point_x - pt[0]) ** 2 + (point_y - pt[1]) ** 2)
            # Return the distance and the direction vector
            return dis, [pt[0] - point[0], pt[1] - point[1]]

    def getPt(point):
        '''
        Converts a point object to a list of coordinates.
        
        Parameters:
            point: A geometric point object.
        
        Returns:
            A list containing the coordinates of the point.
        '''
        mapping_result = mapping(point)  # Map the point to its coordinates
        return [mapping_result["coordinates"][0], mapping_result["coordinates"][1]]  # Return the coordinates

    def getPolyEdges(poly):
        '''
        Retrieves the edges of a polygon defined by a list of points.
        
        Parameters:
            poly: A list of points defining the polygon.
        
        Returns:
            edges: A list of edges where each edge is defined by a pair of points.
        '''
        edges = []  # Initialize an empty list for edges
        for index, point in enumerate(poly):
            # Connect each point to the next, wrapping around to the first point
            if index < len(poly) - 1:
                edges.append([poly[index], poly[index + 1]])
            else:
                edges.append([poly[index], poly[0]])  # Connect last point to the first
        return edges  # Return the list of edges


    def pointPrecisionChange(pt, num):
        # Change the precision of the point coordinates to the specified number of decimal places.
        return [round(pt[0], num), round(pt[1], num)]

    def linePrecisionChange(line, num):
        # Change the precision of the line endpoints' coordinates to the specified number of decimal places.
        return [GeoFunc.pointPrecisionChange(line[0], num), GeoFunc.pointPrecisionChange(line[1], num)]

    def lineToVec(edge):
        # Convert a line segment (edge) to a vector.
        return [edge[1][0] - edge[0][0], edge[1][1] - edge[0][1]]

    '''May need to encapsulate with approximate calculations'''
    def judgePosition(edge1, edge2):
        # Determine the relative position of two edges.
        x1 = edge1[1][0] - edge1[0][0]
        y1 = edge1[1][1] - edge1[0][1]
        x2 = edge2[1][0] - edge2[0][0]
        y2 = edge2[1][1] - edge2[0][1]
        res = x1 * y2 - x2 * y1
        right = False
        left = False
        parallel = False
        # print("res:", res)
        if res == 0:
            parallel = True
        elif res > 0:
            left = True
        else:
            right = True 
        return right, left, parallel

    def getSlideLine(line, x, y):
        # Translate a line by (x, y) offsets.
        new_line = []
        for pt in line:
            new_line.append([pt[0] + x, pt[1] + y])
        return new_line

    def getCentroid(poly):
        # Get the centroid of the polygon.
        return GeoFunc.getPt(Polygon(poly).centroid)

    # Add rotation calculation functionality
    def rotatePoly(polygon, angle):
        """
        Rotate a polygon around its geometric center.

        :param polygon: The polygon to be rotated, in the format [[x1, y1], [x2, y2], ...]
        :param angle: The angle of rotation, in degrees.
        :return: The rotated polygon, in the format [[x1', y1'], [x2', y2'], ...]
        """
        # Convert angle to radians
        radians = math.radians(angle)
        
        # Calculate the geometric center (centroid)
        centroid_x = sum([point[0] for point in polygon]) / len(polygon)
        centroid_y = sum([point[1] for point in polygon]) / len(polygon)
        
        # Rotate each vertex of the polygon
        rotated_polygon = []
        for x, y in polygon:
            qx = centroid_x + math.cos(radians) * (x - centroid_x) - math.sin(radians) * (y - centroid_y)
            qy = centroid_y + math.sin(radians) * (x - centroid_x) + math.cos(radians) * (y - centroid_y)
            rotated_polygon.append([qx, qy])
        
        return rotated_polygon

    def rotateGroupPoly(polygon, polygon_group, angle):
        """
        Rotate all geometries in polygon_group around the geometric center of polygon.

        :param polygon: The reference polygon (in list format, containing coordinate points)
        :param polygon_group: The group of polygons (in list format, where each element is a polygon)
        :param angle: The angle of rotation (in degrees)
        :return: The rotated group of polygons (in list format)
        """
        # Convert angle to radians
        radians = math.radians(angle)
        
        # Calculate the geometric center (centroid)
        centroid_x = sum([point[0] for point in polygon]) / len(polygon)
        centroid_y = sum([point[1] for point in polygon]) / len(polygon)
        
        # Rotate each polygon in the group
        rotated_polygons = []
        for poly in polygon_group:
            rotated_polygon = []
            for x, y in poly:
                qx = centroid_x + math.cos(radians) * (x - centroid_x) - math.sin(radians) * (y - centroid_y)
                qy = centroid_y + math.sin(radians) * (x - centroid_x) + math.cos(radians) * (y - centroid_y)
                rotated_polygon.append([qx, qy])
            rotated_polygons.append(rotated_polygon)
        
        return rotated_polygons

    def rotateToLowestCenterOfGravity(polygon, step=1):
        # Rotate the polygon to find the orientation with the lowest center of gravity height.
        best_polygon = polygon
        min_cg_height = float('inf')

        for rotation in range(0, 360, step):
            rotated_polygon = GeoFunc.rotatePoly(polygon, rotation)
            
            # Calculate the coordinates of the centroid
            cg_x = sum(point[0] for point in rotated_polygon) / len(rotated_polygon)
            cg_y = sum(point[1] for point in rotated_polygon) / len(rotated_polygon)
            
            # Calculate the y-coordinate of the lowest point of the geometry
            min_y = min(point[1] for point in rotated_polygon)
            
            # Calculate the height of the centroid
            cg_height = cg_y - min_y
            
            # Find the rotation with the lowest centroid height
            if cg_height < min_cg_height:
                min_cg_height = cg_height
                best_polygon = rotated_polygon

        return best_polygon

    def get_edge_vector(p1, p2):
        """Calculate the direction vector from point p1 to p2."""
        return np.array([p2[0] - p1[0], p2[1] - p1[1]])

    def get_unit_vector(vector):
        """Return the unit vector."""
        norm = np.linalg.norm(vector)
        return vector / norm if norm != 0 else vector

    def angle_bisector_vector(edge1_vector, edge2_vector):
        """Calculate the direction vector of the angle bisector between two edges."""
        unit_edge1 = GeoFunc.get_unit_vector(edge1_vector)
        unit_edge2 = GeoFunc.get_unit_vector(edge2_vector)
        bisector = unit_edge1 + unit_edge2
        return GeoFunc.get_unit_vector(bisector)

    def expand_convex_hull(polygon, distance=7.5):
        """
        Translate each vertex of the polygon using the direction of the angle bisector of adjacent edges.
        
        :param polygon: The polygon, in the form of a Polygon object or a list of its points.
        :param distance: The distance to translate the vertices, default is 7.5mm.
        :return: The expanded polygon.
        """
        if isinstance(polygon, list):
            polygon = Polygon(polygon)
        
        points = list(polygon.exterior.coords)
        n = len(points) - 1  # Exclude the last point that is repeated with the first
        expanded_points = []

        for i in range(n):
            # TODO: Needs modification
            
            p1 = points[(i + len(points) - 1) % len(points)]  # Previous vertex
            p2 = points[i]      # Current vertex
            p3 = points[(i + len(points) + 1) % len(points)]  # Next vertex

            # Calculate the direction vectors of the two adjacent edges
            edge1_vector = GeoFunc.get_edge_vector(p1, p2)
            edge2_vector = GeoFunc.get_edge_vector(p3, p2)

            # Calculate the angle at the vertex
            dot_product = np.dot(GeoFunc.get_unit_vector(edge1_vector), GeoFunc.get_unit_vector(edge2_vector))
            angle = np.arccos(dot_product)  # Angle in radians

            # Calculate the translation distance for each vertex based on the formula
            move_distance = distance / np.cos(np.pi / 2 - angle / 2)

            # Calculate the direction of the angle bisector
            bisector_vector = GeoFunc.angle_bisector_vector(edge1_vector, edge2_vector)

            # Move the current vertex along the direction of the angle bisector by the specified distance
            new_point = [float(p2[0] + bisector_vector[0] * move_distance), float(p2[1] + bisector_vector[1] * move_distance)]
            expanded_points.append(new_point)

        return expanded_points

    def restore_convex_hull(expanded_polygon, distance=7.5):
        """
        Restore the expanded polygon by contracting it using a negative distance.
        
        :param expanded_polygon: The expanded polygon, in the form of a Polygon object or a list of its points.
        :param distance: The distance to translate the vertices, default is 7.5mm.
        :return: The restored polygon.
        """
        return GeoFunc.expand_convex_hull(expanded_polygon, -distance)
