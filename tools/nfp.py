from tools.geofunc import GeoFunc
from shapely.geometry import Polygon, Point, mapping, LineString
from shapely.ops import unary_union
import pandas as pd
import json
import copy

class NFP(object):
    def __init__(self, poly1, poly2, **kw):
        self.stationary = copy.deepcopy(poly1)  # Make a deep copy of the first polygon
        self.sliding = copy.deepcopy(poly2)      # Make a deep copy of the second polygon
        start_point_index = GeoFunc.checkBottom(self.stationary)  # Check the bottom point of the stationary polygon
        self.start_point = [poly1[start_point_index][0], poly1[start_point_index][1]]  # Starting point
        self.locus_index = GeoFunc.checkTop(self.sliding)  # Check the top point of the sliding polygon
        # Convert to list to avoid referencing original top point
        self.original_top = list(self.sliding[self.locus_index])  
        GeoFunc.slideToPoint(self.sliding, self.sliding[self.locus_index], self.start_point)  # Slide to starting point
        self.start = True  # Flag to indicate if it's the initial state
        self.nfp = []  # List to store non-overlapping polygon
        self.rectangle = False  # Flag for rectangle detection
        if 'rectangle' in kw:
            if kw["rectangle"] == True:
                self.rectangle = True  # Set rectangle flag if specified
        self.error = 1  # Error tracking variable
        self.main()  # Execute main logic
        # Slide back to the original position after calculations
        GeoFunc.slideToPoint(self.sliding, self.sliding[self.locus_index], self.original_top)

    def main(self):
        i = 0
        if self.rectangle:  # If it's a rectangle, perform quick computation
            width = self.sliding[1][0] - self.sliding[0][0]  # Calculate width
            height = self.sliding[3][1] - self.sliding[0][1]  # Calculate height
            self.nfp.append([self.stationary[0][0], self.stationary[0][1]])  # Add starting point
            self.nfp.append([self.stationary[1][0] + width, self.stationary[1][1]])  # Add next point
            self.nfp.append([self.stationary[2][0] + width, self.stationary[2][1] + height])  # Add next point
            self.nfp.append([self.stationary[3][0], self.stationary[3][1] + height])  # Add final point
        else:
            while self.judgeEnd() == False and i < 75:  # Exit if over 75 iterations, generally due to errors
                # while i < 7:
                # print("########第",i,"轮##########")  # Debugging output for iteration number
                touching_edges = self.detectTouching()  # Detect touching edges
                all_vectors = self.potentialVector(touching_edges)  # Get potential vectors
                if len(all_vectors) == 0:
                    print("没有可行向量")  # No feasible vector
                    self.error = -2  # No feasible vector error
                    break

                vector = self.feasibleVector(all_vectors, touching_edges)  # Select a feasible vector
                if vector == []:
                    print("没有计算出可行向量")  # No feasible vector computed
                    self.error = -5  # No feasible vector error
                    break
                
                self.trimVector(vector)  # Trim the vector
                if vector == [0, 0]:
                    print("未进行移动")  # No movement occurred
                    self.error = -3  # No movement error
                    break

                GeoFunc.slidePoly(self.sliding, vector[0], vector[1])  # Slide the sliding polygon
                self.nfp.append([self.sliding[self.locus_index][0], self.sliding[self.locus_index][1]])  # Append new position
                i = i + 1  # Increment iteration count
                inter = Polygon(self.sliding).intersection(Polygon(self.stationary))  # Check intersection
                if GeoFunc.computeInterArea(inter) > 1:  # If the area of intersection is significant
                    # print("出现相交区域")  # Intersection area detected
                    raise Exception(f'出现相交区域')  # Raise an exception for intersection
                    self.error = -4  # Intersection error
                    break                

        if i == 75:
            print("超出计算次数")  # Exceeded calculation iterations
            self.error = -1  # Exceeded calculation error
    
    # Detect mutual touching edges
    def detectTouching(self):
        touch_edges = []  # List to store touching edges
        stationary_edges, sliding_edges = self.getAllEdges()  # Get all edges of both polygons
        for edge1 in stationary_edges:  # Iterate over stationary edges
            for edge2 in sliding_edges:  # Iterate over sliding edges
                inter = GeoFunc.intersection(edge1, edge2)  # Check for intersection between edges
                if inter != []:  # If there is an intersection
                    pt = [inter[0], inter[1]]  # Intersection point
                    edge1_bound = (GeoFunc.almostEqual(edge1[0], pt) or GeoFunc.almostEqual(edge1[1], pt))  # Check if it's a boundary point
                    edge2_bound = (GeoFunc.almostEqual(edge2[0], pt) or GeoFunc.almostEqual(edge2[1], pt))  # Check if it's a boundary point
                    stationary_start = GeoFunc.almostEqual(edge1[0], pt)  # Check if it starts from the stationary edge
                    orbiting_start = GeoFunc.almostEqual(edge2[0], pt)  # Check if it starts from the sliding edge
                    touch_edges.append({
                        "edge1": edge1,
                        "edge2": edge2,
                        "vector1": self.edgeToVector(edge1),  # Convert edge1 to vector
                        "vector2": self.edgeToVector(edge2),  # Convert edge2 to vector
                        "edge1_bound": edge1_bound,  # Boundary status of edge1
                        "edge2_bound": edge2_bound,  # Boundary status of edge2
                        "stationary_start": stationary_start,  # Starting status of stationary edge
                        "orbiting_start": orbiting_start,  # Starting status of sliding edge
                        "pt": [inter[0], inter[1]],  # Intersection point
                        "type": 0  # Type identifier for touching edges
                    })
        return touch_edges  # Return list of touching edges

    # Obtain potential transferable vectors
    def potentialVector(self, touching_edges):
        all_vectors = []  # List to store all potential vectors
        for touching in touching_edges:  # Iterate over touching edges
            # print("touching:", touching)
            aim_edge = []  # Initialize the aim edge variable
            # Situation 1
            if touching["edge1_bound"] == True and touching["edge2_bound"] == True:  # Both edges are boundaries
                right, left, parallel = GeoFunc.judgePosition(touching["edge1"], touching["edge2"])  # Judge their position
                # print("right,left,parallel:", right, left, parallel)
                if touching["stationary_start"] == True and touching["orbiting_start"] == True:
                    touching["type"] = 0  # Set type
                    if left == True:
                        aim_edge = [touching["edge2"][1], touching["edge2"][0]]  # Reverse direction
                    if right == True:
                        aim_edge = touching["edge1"]  # Keep direction
                if touching["stationary_start"] == True and touching["orbiting_start"] == False:
                    touching["type"] = 1  # Set type
                    if left == True:
                        aim_edge = touching["edge1"]  # Keep direction
                if touching["stationary_start"] == False and touching["orbiting_start"] == True:
                    touching["type"] = 2  # Set type
                    if right == True:
                        aim_edge = [touching["edge2"][1], touching["edge2"][0]]  # Reverse direction
                if touching["stationary_start"] == False and touching["orbiting_start"] == False:
                    touching["type"] = 3  # Set type
    
            # Situation 2
            if touching["edge1_bound"] == False and touching["edge2_bound"] == True:  # Edge1 is not a boundary, edge2 is
                aim_edge = [touching["pt"], touching["edge1"][1]]  # Aim towards edge1
                touching["type"] = 4  # Set type
            
            # Situation 3
            if touching["edge1_bound"] == True and touching["edge2_bound"] == False:  # Edge1 is a boundary, edge2 is not
                aim_edge = [touching["edge2"][1], touching["pt"]]  # Aim towards edge2
                touching["type"] = 5  # Set type

            if aim_edge != []:  # If aim_edge has been set
                vector = self.edgeToVector(aim_edge)  # Convert aim_edge to vector
                if self.detectExisting(all_vectors, vector) == False:  # Check if the vector already exists
                    all_vectors.append(vector)  # Add new vector

        return all_vectors  # Return list of potential vectors

    
    def detectExisting(self,vectors,judge_vector):
        # Check if the judge_vector already exists in the vectors
        for vector in vectors:
            # Compare each vector with the judge_vector using a nearly equal function
            if GeoFunc.almostEqual(vector,judge_vector):
                return True  # Return True if a match is found
        return False  # Return False if no match is found
    
    def edgeToVector(self,edge):
        # Convert an edge (defined by two points) into a vector representation
        return [edge[1][0]-edge[0][0],edge[1][1]-edge[0][1]]
    
    # Select feasible vectors based on various conditions
    def feasibleVector(self,all_vectors,touching_edges):
        '''
        This section of code needs refactoring; it is too complex.
        '''
        res_vector=[]  # Initialize the result vector as an empty list
        # print("\nall_vectors:",all_vectors)
        for vector in all_vectors:
            feasible=True  # Assume the vector is feasible initially
            # print("\nvector:",vector,"\n")
            for touching in touching_edges:
                vector1=[]  # Initialize vector1
                vector2=[]  # Initialize vector2
                # Determine the direction of vector1 based on the touching edge's stationary state
                if touching["stationary_start"]==True:
                    vector1=touching["vector1"]
                else:
                    vector1=[-touching["vector1"][0],-touching["vector1"][1]]  # Reverse direction if not stationary
                
                # Determine the direction of vector2 based on the touching edge's orbiting state
                if touching["orbiting_start"]==True:
                    vector2=touching["vector2"]
                else:
                    vector2=[-touching["vector2"][0],-touching["vector2"][1]]  # Reverse direction if not orbiting
                
                # Calculate the cross products to determine relative orientations
                vector12_product=GeoFunc.crossProduct(vector1,vector2) # Cross product of vector1 and vector2
                vector_vector1_product=GeoFunc.crossProduct(vector1,vector) # Cross product of vector1 and current vector
                vector_vector2_product=GeoFunc.crossProduct(vector2,vector) # Cross product of vector2 and current vector
                
                # Handle special cases based on the type of touching edge
                if touching["type"]==4 and (vector_vector1_product*vector12_product)<0:
                    feasible=False  # Infeasible if the conditions for type 4 are not met
                if touching["type"]==5 and (vector_vector2_product*(-vector12_product))>0:
                    feasible=False  # Infeasible if the conditions for type 5 are not met
                
                # Check the normal case for feasibility based on orientations
                if vector12_product>0:
                    if vector_vector1_product<0 and vector_vector2_product<0:
                        feasible=False  # Both vectors pointing in opposite directions
                if vector12_product<0:
                    if vector_vector1_product>0 and vector_vector2_product>0:
                        feasible=False  # Both vectors pointing in the same direction
                
                # In case of parallel vectors, perform additional checks
                if vector12_product==0:
                    inter=GeoFunc.newLineInter(touching["edge1"],touching["edge2"])  # Check intersection of edges
                    if inter["geom_type"]=="LineString":
                        if inter["length"]>0.01:  # If the intersection line has length greater than threshold
                            # If there is an intersection, it needs to be on the left side
                            if (touching["orbiting_start"]==True and vector_vector2_product<0) or (touching["orbiting_start"]==False and vector_vector2_product>0):
                                feasible=False  # Mark as infeasible if conditions are not satisfied
                    else:
                        # If the directions are the same and the transformed line is also parallel, it cannot take direction a
                        if touching["orbiting_start"]==True != touching["stationary_start"]==False and vector_vector1_product==0:
                            if touching["vector1"][0]*vector[0]>0: # Check if the directions are the same
                                feasible=False  # Mark as infeasible if conditions are met

            if feasible==True:  # If a feasible vector is found
                res_vector=vector  # Set the result vector
                break  # Exit the loop after finding the first feasible vector
        return res_vector  # Return the feasible vector found
        
    # Trim overly long vectors to fit within constraints
    def trimVector(self,vector):
        stationary_edges,sliding_edges=self.getAllEdges()  # Get all edges from stationary and sliding polygons
        new_vectors=[]  # Initialize a list to hold new vectors
        for pt in self.sliding:  # Iterate over all sliding points
            for edge in stationary_edges:  # Check against all stationary edges
                line_vector=LineString([pt,[pt[0]+vector[0],pt[1]+vector[1]]])  # Create a line vector
                end_pt=[pt[0]+vector[0],pt[1]+vector[1]]  # Calculate the endpoint of the line vector
                line_polygon=LineString(edge)  # Convert the edge to a LineString
                inter=line_vector.intersection(line_polygon)  # Find intersection between the line vector and the edge
                if inter.geom_type=="Point":  # Check if intersection is a point
                    inter_mapping=mapping(inter)  # Get mapping of the intersection point
                    inter_coor=inter_mapping["coordinates"]  # Get coordinates of the intersection
                    # Ensure the intersection point is significantly different from both the endpoint and the original point
                    if (abs(end_pt[0]-inter_coor[0])>0.01 or abs(end_pt[1]-inter_coor[1])>0.01) and (abs(pt[0]-inter_coor[0])>0.01 or abs(pt[1]-inter_coor[1])>0.01):
                        new_vectors.append([inter_coor[0]-pt[0],inter_coor[1]-pt[1]])  # Append the vector from pt to the intersection

        for pt in self.stationary:  # Iterate over all stationary points
            for edge in sliding_edges:  # Check against all sliding edges
                line_vector=LineString([pt,[pt[0]-vector[0],pt[1]-vector[1]]])  # Create a line vector in the opposite direction
                end_pt=[pt[0]-vector[0],pt[1]-vector[1]]  # Calculate the endpoint of the line vector
                line_polygon=LineString(edge)  # Convert the edge to a LineString
                inter=line_vector.intersection(line_polygon)  # Find intersection between the line vector and the edge
                if inter.geom_type=="Point":  # Check if intersection is a point
                    inter_mapping=mapping(inter)  # Get mapping of the intersection point
                    inter_coor=inter_mapping["coordinates"]  # Get coordinates of the intersection
                    # Ensure the intersection point is significantly different from both the endpoint and the original point
                    if (abs(end_pt[0]-inter_coor[0])>0.01 or abs(end_pt[1]-inter_coor[1])>0.01) and (abs(pt[0]-inter_coor[0])>0.01 or abs(pt[1]-inter_coor[1])>0.01):
                        new_vectors.append([pt[0]-inter_coor[0],pt[1]-inter_coor[1]])  # Append the vector from pt to the intersection
        
        # print(new_vectors)
        for vec in new_vectors:  # Iterate over newly generated vectors
            if abs(vec[0])<abs(vector[0]) or abs(vec[1])<abs(vector[1]):  # Check if the new vector is shorter
                # print(vec)
                vector[0]=vec[0]  # Update the vector with the new shorter vector
                vector[1]=vec[1]  # Update the vector with the new shorter vector
    
    # Get all edges of the two polygons
    def getAllEdges(self):
        return GeoFunc.getPolyEdges(self.stationary),GeoFunc.getPolyEdges(self.sliding)  # Return edges of both polygons
    
    # Check if the process has ended
    def judgeEnd(self):
        sliding_locus=self.sliding[self.locus_index]  # Get the current sliding locus
        main_bt=self.start_point  # Reference to the starting point
        # Check if the sliding locus is close to the starting point
        if abs(sliding_locus[0]-main_bt[0])<0.1 and abs(sliding_locus[1]-main_bt[1])<0.1:
            if self.start==True:  # If starting point is active
                self.start=False  # Mark starting point as inactive
                # print("Check if ended: No")
                return False  # Not ended yet
            else:
                # print("Check if ended: Yes")
                return True  # Ended
        else:
            # print("Check if ended: No")
            return False  # Not ended

    # Calculate penetration depth
    def getDepth(self):
        '''
        Calculate the distance from poly2's checkTop to the NFP.
        '''
        d1=Polygon(self.nfp).distance(Point(self.original_top))
        # if point is inside polygon, d1=0
        # d2: distance from the point to nearest boundary
        if d1==0:
            d2=Polygon(self.nfp).boundary.distance(Point(self.original_top))
            # print('d2:',d2)
            return d2
        else: 
            return 0