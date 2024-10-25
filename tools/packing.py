from tools.nfp import NFP
from shapely.geometry import Polygon, Point, mapping, LineString
from shapely.ops import unary_union
from shapely import affinity
from tools.geofunc import GeoFunc
import pyclipper 
import math
import numpy as np
import pandas as pd
import json
import matplotlib.pyplot as plt
import csv
import logging
import random
import copy
import os

class PackingUtil(object):
    
    @staticmethod
    def getInnerFitRectangle(poly, x, y):
        # Get the boundaries of the polygon
        left_index, bottom_index, right_index, top_index = GeoFunc.checkBound(poly)  
        # Get the translated result of the polygon by sliding it to the origin
        new_poly = GeoFunc.getSlide(poly, -poly[left_index][0], -poly[bottom_index][1])  

        # Define a reference point from the top edge of the translated polygon
        refer_pt = [new_poly[top_index][0], new_poly[top_index][1]]
        # Calculate the width and height for the inner fitting rectangle
        ifr_width = x - new_poly[right_index][0]
        ifr_height = y - new_poly[top_index][1]

        # Define the inner fitting rectangle using the reference point and calculated dimensions
        IFR = [refer_pt, 
               [refer_pt[0] + ifr_width, refer_pt[1]], 
               [refer_pt[0] + ifr_width, refer_pt[1] + ifr_height], 
               [refer_pt[0], refer_pt[1] + ifr_height]]
        return IFR  # Return the coordinates of the inner fitting rectangle
    
class NFPAssistant(object):
    def __init__(self, polys, **kw):
        # Remove redundancy from the input polygons and store them
        self.polys = PolyListProcessor.deleteRedundancy(copy.deepcopy(polys))  
        # Initialize lists to store area, first vector, and centroid of each polygon
        self.area_list, self.first_vec_list, self.centroid_list = [], [], []  
        for poly in self.polys:
            P = Polygon(poly)  # Create a Polygon object from the polygon coordinates
            # Get the centroid of the polygon and store it
            self.centroid_list.append(GeoFunc.getPt(P.centroid))  
            # Store the area of the polygon
            self.area_list.append(int(P.area))  
            # Store the first vector (direction) of the polygon
            self.first_vec_list.append([poly[1][0] - poly[0][0], poly[1][1] - poly[0][1]])  
        # Initialize a list to store the NFPs (no-fit polygons) for each polygon pair
        self.nfp_list = [[0] * len(self.polys) for i in range(len(self.polys))]  
        self.load_history = False  # Flag to check if history should be loaded
        self.history_path = None  # Path to the history file
        self.history = None  # History data

        # Check if a history path is provided in the keyword arguments
        if 'history_path' in kw:
            self.history_path = kw['history_path']  

        # Check if history should be loaded
        if 'load_history' in kw:
            if kw['load_history'] == True:
                # Load history from memory directly to reduce I/O time
                if 'history' in kw:
                    self.history = kw['history']
                self.load_history = True
                self.loadHistory()  # Load the historical data
        
        self.store_nfp = False  # Flag to check if NFPs should be stored
        # Check if NFPs should be stored
        if 'store_nfp' in kw:
            if kw['store_nfp'] == True:
                self.store_nfp = True
        
        self.store_path = None  # Path to store NFPs
        # Check if a store path is provided
        if 'store_path' in kw:
            self.store_path = kw['store_path']

        # Check if all NFPs should be calculated if not loading history
        if 'get_all_nfp' in kw:
            if kw['get_all_nfp'] == True and self.load_history == False:
                self.getAllNFP()  # Calculate all NFPs
        
    def loadHistory(self):
        # Load NFP history from a CSV file or from the given history
        if not self.history:
            if not self.history_path:
                path = "record/nfp.csv"  # Default path for the history file
            else:
                path = self.history_path
            df = pd.read_csv(path, header=None)  # Read the CSV file into a DataFrame
        else:
            df = self.history  # Use the provided history

        # Populate the NFP list with historical data
        for index in range(df.shape[0]):
            i = self.getPolyIndex(json.loads(df[0][index]))  # Get the index of the first polygon
            j = self.getPolyIndex(json.loads(df[1][index]))  # Get the index of the second polygon
            if i >= 0 and j >= 0:
                self.nfp_list[i][j] = json.loads(df[2][index])  # Assign the NFP data

    # Get the index of a shape based on its geometry
    def getPolyIndex(self, target):
        area = int(Polygon(target).area)  # Calculate the area of the target shape
        # Calculate the first vector (direction) of the target shape
        first_vec = [target[1][0] - target[0][0], target[1][1] - target[0][1]]  
        # Get the area index from the stored area list
        area_index = PolyListProcessor.getIndexMulti(area, self.area_list)  
        if len(area_index) == 1:  # If there is only one matching area
            return area_index[0]  # Return its index
        else:
            # Get the vector index from the stored first vector list
            vec_index = PolyListProcessor.getIndexMulti(first_vec, self.first_vec_list)  
            index = [x for x in area_index if x in vec_index]  # Get the common indices
            if len(index) == 0:
                return -1  # Return -1 if no match is found
            return index[0]  # Usually, there will be only one match
    
    # Calculate all the NFPs for the polygon pairs
    def getAllNFP(self):
        nfp_multi = False  # Flag to enable multi-threading for NFP calculation
        if nfp_multi == True:
            # Generate tasks for all pairs of polygons for multi-threaded processing
            tasks = [(main, adjoin) for main in self.polys for adjoin in self.polys]
            res = pool.starmap(NFP, tasks)  # Execute tasks in parallel
            for k, item in enumerate(res):
                i = k // len(self.polys)  # Determine the index of the main polygon
                j = k % len(self.polys)  # Determine the index of the adjacent polygon
                # Slide the NFP to the centroid of the main polygon
                self.nfp_list[i][j] = GeoFunc.getSlide(item.nfp, -self.centroid_list[i][0], -self.centroid_list[i][1])  
        else:
            # Calculate NFP for each pair of polygons in a single-threaded manner
            for i, poly1 in enumerate(self.polys):
                for j, poly2 in enumerate(self.polys):
                    nfp = NFP(poly1, poly2).nfp  # Calculate the NFP
                    # Store the NFP slid to the centroid of the main polygon
                    self.nfp_list[i][j] = GeoFunc.getSlide(nfp, -self.centroid_list[i][0], -self.centroid_list[i][1])  
        if self.store_nfp == True:
            self.storeNFP()  # Store the calculated NFPs if required
    
    def storeNFP(self):
        # Determine the path to store NFPs
        if self.store_path == None:
            path = "record/nfp.csv"  # Default path for storing NFPs
        else:
            path = self.store_path
        # Open the CSV file in append mode to save the NFPs
        with open(path, "a+") as csvfile:
            writer = csv.writer(csvfile)
            for i in range(len(self.polys)):
                for j in range(len(self.polys)):
                    # Write the polygon pair and their corresponding NFP to the CSV file
                    writer.writerows([[self.polys[i], self.polys[j], self.nfp_list[i][j]]])

    # Get the NFP directly between two shapes
    def getDirectNFP(self, poly1, poly2, **kw):
        if 'index' in kw:
            i = kw['index'][0]  # Get index of the first polygon if provided
            j = kw['index'][1]  # Get index of the second polygon if provided
            centroid = GeoFunc.getPt(Polygon(self.polys[i]).centroid)  # Get the centroid of the first polygon
        else:
            # First get the IDs of poly1 and poly2
            i = self.getPolyIndex(poly1)  # Get index of poly1
            j = self.getPolyIndex(poly2)  # Get index of poly2
            centroid=GeoFunc.getPt(Polygon(poly1).centroid)
        # have nfp been calculated
        if self.nfp_list[i][j]==0:
            nfp=NFP(poly1,poly2).nfp
            #self.nfp_list[i][j]=GeoFunc.getSlide(nfp,-centroid[0],-centroid[1])
            if self.store_nfp==True:
                with open("record/nfp.csv","a+") as csvfile:
                    writer = csv.writer(csvfile)
                    writer.writerows([[poly1,poly2,nfp]])
            return nfp
        else:
            return GeoFunc.getSlide(self.nfp_list[i][j],centroid[0],centroid[1])

class PolyListProcessor(object):
    @staticmethod
    def getPolyObjectList(polys, allowed_rotation):
        '''
        Converts the given polys and allowed rotation angles into a list of Poly objects.
        '''
        poly_list = []
        for i, poly in enumerate(polys):
            poly_list.append(Poly(i, poly, allowed_rotation))  # Create a Poly object and add it to the list
        return poly_list

    @staticmethod
    def getPolysVertices(_list):
        '''Returns the vertices of the polygons; the sorting will affect the result.'''
        polys = []
        for i in range(len(_list)):
            polys.append(_list[i].poly)  # Append the polygon's vertices to the list
        return polys
    
    @staticmethod
    def getPolysVerticesCopy(_list):
        '''Returns a deep copy of the vertices of the polygons; does not affect the shapes in the list.'''
        polys = []
        for i in range(len(_list)):
            polys.append(copy.deepcopy(_list[i].poly))  # Deep copy each polygon's vertices
        return polys

    @staticmethod
    def getPolyListIndex(poly_list):
        index_list = []
        for i in range(len(poly_list)):
            index_list.append(poly_list[i].num)  # Append the index of each polygon to the index list
        return index_list
    
    @staticmethod
    def getIndex(item, _list):
        for i in range(len(_list)):
            if item == _list[i]:
                return i  # Return the index if the item matches
        return -1  # Return -1 if the item is not found
    
    @staticmethod
    def getIndexMulti(item, _list):
        index_list = []
        for i in range(len(_list)):
            if item == _list[i]:
                index_list.append(i)  # Append the index if the item matches
        return index_list  # Return the list of indices that matched

    @staticmethod
    def randomSwap(poly_list, target_id):
        new_poly_list = copy.deepcopy(poly_list)  # Create a deep copy of the polygon list

        swap_with = int(random.random() * len(new_poly_list))  # Select a random index to swap with
        
        item1 = new_poly_list[target_id]  # Get the target polygon
        item2 = new_poly_list[swap_with]  # Get the polygon to swap with
            
        new_poly_list[target_id] = item2  # Swap the polygons
        new_poly_list[swap_with] = item1
        return new_poly_list  # Return the new list with swapped polygons

    @staticmethod
    def randomRotate(poly_list, min_angle, target_id):
        new_poly_list = copy.deepcopy(poly_list)  # Create a deep copy of the polygon list

        index = random.randint(0, len(new_poly_list) - 1)  # Select a random polygon to rotate
        RatotionPoly(min_angle).rotation(new_poly_list[index].poly)  # Rotate the selected polygon
        return new_poly_list  # Return the list with the rotated polygon

    @staticmethod
    def showPolyList(width, poly_list):
        blf = BottomLeftFill(width, PolyListProcessor.getPolysVertices(poly_list))  # Create a BottomLeftFill object
        blf.showAll()  # Display all polygons

    @staticmethod
    def deleteRedundancy(_arr):
        new_arr = []
        for item in _arr:
            if not item in new_arr:  # Check if the item is already in the new array
                new_arr.append(item)  # Append unique items to the new array
        return new_arr  # Return the new array without duplicates

    @staticmethod
    def getPolysByIndex(index_list, poly_list):
        choosed_poly_list = []
        for i in index_list:
            choosed_poly_list.append(poly_list[i])  # Append selected polygons based on the index list
        return choosed_poly_list

class RatotionPoly():
    def __init__(self, angle):
        self.angle = angle  # Set the rotation angle
        self._max = 360 / angle  # Calculate the maximum number of possible rotations

    def rotation(self, poly):
        if self._max > 1:
            # print("Rotating the shape")
            rotation_res = random.randint(1, self._max - 1)  # Randomly select a rotation value
            Poly = Polygon(poly)  # Create a Polygon object from the provided vertices
            new_Poly = affinity.rotate(Poly, rotation_res * self.angle)  # Rotate the polygon
            mapping_res = mapping(new_Poly)  # Map the new polygon back to coordinates
            new_poly = mapping_res["coordinates"][0]  # Extract the new coordinates
            for index in range(0, len(poly)):
                poly[index] = [new_poly[index][0], new_poly[index][1]]  # Update the original polygon with new coordinates
        else:
            pass
            # print("Rotation is not allowed")

    def rotation_specific(self, poly, angle=-1):
        '''
        Rotate the polygon by a specific angle.
        '''
        Poly = Polygon(poly)  # Create a Polygon object from the provided vertices
        if angle == -1: 
            angle = self.angle  # Use the default angle if none is provided
        elif len(angle) > 0:
            angle = np.random.choice(angle)  # Randomly choose an angle from the given options
            # print('Rotating {}°'.format(angle))
        new_Poly = affinity.rotate(Poly, angle)  # Rotate the polygon by the specified angle
        mapping_res = mapping(new_Poly)  # Map the new polygon back to coordinates
        new_poly = mapping_res["coordinates"][0]  # Extract the new coordinates
        for index in range(0, len(poly)):
            poly[index] = [new_poly[index][0], new_poly[index][1]]  # Update the original polygon with new coordinates

