from tools.geofunc import GeoFunc
import tools.packing as packing
from tools.nfp import NFP
from shapely.geometry import Polygon, mapping
from shapely import affinity
import numpy as np, random, operator, pandas as pd, matplotlib.pyplot as plt
import json
import csv
import time
import multiprocessing
import datetime
import random
import copy

class BottomLeftFill(object):
    def __init__(self, width, original_polygons, **kw):
        """
        Initialize the BottomLeftFill algorithm.
        
        Args:
            width (int): The width of the container for nesting polygons.
            original_polygons (list): List of polygons to be placed.
            **kw: Additional arguments, including NFPAssistant for NFP calculation and vertical for vertical layout.
        """
        self.choose_nfp = False
        self.width = width
        self.length = 150000  # The container length
        self.contain_length = 2000
        self.polygons = original_polygons
        self.NFPAssistant = None
        self.rotation = []
        if 'NFPAssistant' in kw:
            self.NFPAssistant = kw["NFPAssistant"]
        self.vertical = False
        if 'vertical' in kw:
            self.vertical = kw['vertical']
        
        print("Total Num:", len(original_polygons))
        self.placeFirstPoly()
        for i in range(1, len(self.polygons)):
            print("############################## Place the ", i + 1, "th shape #################################")
            self.placePoly(i)
        
        self.getLength()

    def placeFirstPoly(self):
        """
        Place the first polygon at the best position with the optimal rotation.
        """
        best_height = float('inf')
        best_rotation = None
        best_position = None
        errors = 0

        for rotation in range(0, 360, 90):
            try:
                rotated_poly = GeoFunc.rotatePoly(self.polygons[0], rotation)
                # Check if rotated polygon exceeds the width limit
                min_x = min([point[1] for point in rotated_poly])
                max_x = max([point[1] for point in rotated_poly])
                if max_x - min_x > self.width:
                    raise Exception(f"Rotation {rotation} degrees exceeds width limit, skipping.")

                left_index, bottom_index, right_index, top_index = GeoFunc.checkBound(rotated_poly)  # Get the polygon's bounds        
                GeoFunc.slidePoly(rotated_poly, -rotated_poly[left_index][0], -rotated_poly[bottom_index][1])  # Move to bottom-left corner
                
                height = max([point[0] for point in rotated_poly])
                if height < best_height:
                    best_height = height
                    best_rotation = rotation
                    best_position = rotated_poly
            except Exception as e:
                print(f"Error at rotation {rotation} degrees: {e}")
                errors += 1

        if errors == 4:
            raise Exception("Failed to place the first polygon in any rotation.")
        
        # Place the polygon with the best rotation
        print('best rotation:', best_rotation)
        self.rotation.append(best_rotation)
        self.polygons[0] = best_position

    def placePoly(self, index):
        """
        Place the polygon at the given index in the optimal position with the best rotation.
        
        Args:
            index (int): The index of the polygon to place.
        """
        best_height = float('inf')
        best_position = None
        best_rotation = None
        errors = 0
        
        for rotation in range(0, 360, 90):
            try:
                rotated_poly = GeoFunc.rotatePoly(self.polygons[index], rotation)  # Rotate the polygon
                # Check if rotated polygon exceeds the width limit
                min_x = min([point[1] for point in rotated_poly])
                max_x = max([point[1] for point in rotated_poly])
                if max_x - min_x > self.width:
                    raise Exception(f"Rotation {rotation} degrees exceeds width limit, skipping.")
                adjoin = rotated_poly
                differ_region = self.getDifferRegion(index, adjoin)
                differ_index = self.getBottomLeft(differ_region)
                refer_pt_index = GeoFunc.checkTop(adjoin)
                GeoFunc.slideToPoint(adjoin, adjoin[refer_pt_index], differ_region[differ_index])
                
                height = max([point[0] for point in rotated_poly])
                if height < best_height:
                    best_height = height
                    best_position = differ_region[differ_index]
                    best_rotation = rotation
                    best_poly = adjoin
            except Exception as e:
                print(f"Error at rotation {rotation} degrees for polygon {index + 1}: {e}")
                errors += 1

        if errors == 4:
            raise Exception(f"Failed to place polygon {index + 1} in any rotation.")
        
        # Place the polygon with the best rotation
        print('best rotation: ', best_rotation)
        self.rotation.append(best_rotation)
        self.polygons[index] = best_poly

    def getDifferRegion(self, index, adjoin):
        """
        Get the difference region where the polygon can be placed.
        
        Args:
            index (int): The index of the polygon being placed.
            adjoin (Polygon): The polygon being placed.
        
        Returns:
            list: The coordinates of the difference region.
        """
        # Check for vertical layout
        if self.vertical:
            ifr = packing.PackingUtil.getInnerFitRectangle(adjoin, self.width, self.length)
        else:
            ifr = packing.PackingUtil.getInnerFitRectangle(adjoin, self.length, self.width)            
        differ_region = Polygon(ifr)
        
        for main_index in range(0, index):
            main = self.polygons[main_index]
            if self.NFPAssistant is None:
                nfp = NFP(main, adjoin).nfp
            else:
                nfp = self.NFPAssistant.getDirectNFP(main, adjoin)
            nfp_poly = Polygon(nfp)
            try:
                differ_region = differ_region.difference(nfp_poly)
            except:
                print('NFP failure, areas of polygons are:')
                self.showAll()
                for poly in main, adjoin:
                    print(Polygon(poly).area)
                self.showPolys([main] + [adjoin] + [nfp])
                print('NFP loaded from: ', self.NFPAssistant.history_path)
                
        return GeoFunc.polyToArr(differ_region)

    def getBottomLeft(self, poly):
        """
        Get the bottom-left corner point of a polygon. If multiple points have the same x-coordinate, 
        select the one with the lowest y-coordinate.
        
        Args:
            poly (list): The polygon coordinates.
        
        Returns:
            int: The index of the bottom-left point.
        """
        bl=[]  # List of bottom-left points
        _min = 999999
        # Find the leftmost point
        for i, pt in enumerate(poly):
            pt_object = {
                    "index": i,
                    "x": pt[0],
                    "y": pt[1]
            }
            if self.vertical:
                target = pt[1]
            else:
                target = pt[0]
            if target < _min:
                _min = target
                bl = [pt_object]
            elif target == _min:
                bl.append(pt_object)
        if len(bl) == 1:
            return bl[0]["index"]
        else:
            if self.vertical:
                target = "x"                
            else:
                target = "y"
            _min = bl[0][target]
            one_pt = bl[0]
            for pt_index in range(1, len(bl)):
                if bl[pt_index][target] < _min:
                    one_pt = bl[pt_index]
                    _min = one_pt["y"]
            return one_pt["index"]

    def getLength(self):
        """
        Calculate and return the current length of the container based on the placed polygons.
        
        Returns:
            int: The maximum length of the container.
        """
        _max = 0
        for i in range(0, len(self.polygons)):
            if self.vertical:
                extreme_index = GeoFunc.checkTop(self.polygons[i])
                extreme = self.polygons[i][extreme_index][1]
            else:
                extreme_index = GeoFunc.checkRight(self.polygons[i])
                extreme = self.polygons[i][extreme_index][0]
            if extreme > _max:
                _max = extreme
        self.contain_length = _max
        return _max

class TOPOS(object):
    """
    TOPOS heuristic algorithm: Iteratively places shapes and dynamically adjusts their positions.
    This algorithm is based on Bennell's TOPOS Revised.
    """
    def __init__(self, original_polys, width):
        """
        Initialize the TOPOS algorithm.
        
        Args:
            original_polys (list): List of original polygons to be placed.
            width (int): The width of the container for nesting.
        """
        self.polys = original_polys
        self.cur_polys = []
        self.width = width
        self.NFPAssistant = packing.NFPAssistant(self.polys, store_nfp=False, get_all_nfp=True, load_history=True)
        
        self.run()

    def run(self):
        # Slide the first polygon to position (1000, 1000) and add it to the current polygon list
        self.cur_polys.append(GeoFunc.getSlide(self.polys[0], 1000, 1000)) 
        
        # Initialize the four boundaries of the enclosing rectangle
        self.border_left, self.border_right, self.border_bottom, self.border_top = 0, 0, 0, 0 
        self.border_height, self.border_width = 0, 0
        
        # Iterate through the remaining polygons and place them
        for i in range(1, len(self.polys)):
            # Update the current boundary information
            self.updateBound()

            # Initialize feasible border with the first polygon
            feasible_border = Polygon(self.cur_polys[0])
            
            # Calculate the NFP (No-Fit Polygon) between the current polygon and the placed polygons
            for fixed_poly in self.cur_polys:
                # Get the NFP between the fixed polygon and the current polygon, and update the feasible boundary
                nfp = self.NFPAssistant.getDirectNFP(fixed_poly, self.polys[i])
                feasible_border = feasible_border.union(Polygon(nfp))  # Union the NFP into the feasible boundary
            
            # Get all feasible placement points
            feasible_point = self.chooseFeasiblePoint(feasible_border)
            
            # Calculate the left and right boundary widths of the current polygon
            poly_left_pt, poly_bottom_pt, poly_right_pt, poly_top_pt = GeoFunc.checkBoundPt(self.polys[i])
            poly_left_width, poly_right_width = poly_top_pt[0] - poly_left_pt[0], poly_right_pt[0] - poly_top_pt[0]

            # Iterate through all feasible points and select the position with the smallest change
            min_change = float('inf')  # Initialize min_change to a large number
            target_position = []  # Variable to hold the best target position
            
            for pt in feasible_point:
                change = min_change  # Store the current min_change for comparison
                
                # Check if the shape does not exceed the boundaries
                if pt[0] - poly_left_width >= self.border_left and pt[0] + poly_right_width <= self.border_right:
                    # The shape does not exceed the boundary; min_change is negative
                    change = min(self.border_left - pt[0], self.border_left - pt[0])
                elif min_change > 0:
                    # The shape exceeds the left or right boundary; choose the larger change value
                    change = max(self.border_left - pt[0] + poly_left_width, pt[0] + poly_right_width - self.border_right)
                else:
                    # If exceeding and min_change <= 0, no need to change
                    pass

                # Update the target position if a smaller change is found
                if change < min_change:
                    min_change = change
                    target_position = pt  # Update target_position to the current feasible point
                
            # Slide the polygon to the final position based on the selected target position
            reference_point = self.polys[i][GeoFunc.checkTop(self.polys[i])]
            self.cur_polys.append(GeoFunc.getSlide(self.polys[i], target_position[0] - reference_point[0], target_position[1] - reference_point[1]))

        # Move all polygons to the bottom left corner
        self.slideToBottomLeft()
        # Display the result of the arrangement
        self.showResult()

        
    def updateBound(self):
        '''
        Update the enclosing rectangle boundaries based on the current polygons
        '''
        # Get the boundary values of the last polygon in the current polygons
        border_left, border_bottom, border_right, border_top = GeoFunc.checkBoundValue(self.cur_polys[-1])
        
        # Update the overall boundaries if the new polygon exceeds or is less than current boundaries
        if border_left < self.border_left:
            self.border_left = border_left
        if border_bottom < self.border_bottom:
            self.border_bottom = border_bottom
        if border_right > self.border_right:
            self.border_right = border_right
        if border_top > self.border_top:
            self.border_top = border_top
        
        # Update the height and width of the enclosing rectangle
        self.border_height = self.border_top - self.border_bottom
        self.border_width = self.border_right - self.border_left
        
    def chooseFeasiblePoint(self, border):
        '''Select feasible points based on the given border'''
        res = mapping(border)  # Map the border to a geometric representation
        _arr = []  # Array to hold feasible points
        
        # If the result is a MultiPolygon, extract feasible points from each polygon
        if res["type"] == "MultiPolygon":
            for poly in res["coordinates"]:
                _arr = _arr + self.feasiblePoints(poly)
        else:
            _arr = _arr + self.feasiblePoints(res["coordinates"][0])  # Single polygon case
        
        return _arr  # Return the array of feasible points
        
    def feasiblePoints(self, poly):
        '''
        1. Convert the Polygon object to points
        2. Exclude points exceeding the width range
        3. Include intersection points with the boundary
        '''
        result = []  # Initialize result list to hold feasible points
        for pt in poly:
            # Check if the point is feasible based on the upper boundary
            feasible1 = pt[1] - self.border_top > 0 and pt[1] - self.border_top + self.border_height <= self.width
            # Check if the point is feasible based on the lower boundary
            feasible2 = self.border_bottom - pt[1] > 0 and self.border_bottom - pt[1] + self.border_height <= self.width
            # Check if the point is within the top and bottom boundaries
            feasible3 = pt[1] <= self.border_top and pt[1] >= self.border_bottom
            
            # If any of the feasibility conditions are met, add the point to the result
            if feasible1 or feasible2 or feasible3:
                result.append([pt[0], pt[1]])
        return result  # Return the list of feasible points

    def slideToBottomLeft(self):
        '''Move all polygons to the bottom left corner'''
        for poly in self.cur_polys:
            GeoFunc.slidePoly(poly, -self.border_left, -self.border_bottom)  # Slide each polygon
