import os
import csv
import json

def read_polygon_from_txt(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()  # Read all lines from the file
        num_polygons = 1  # Set the number of polygons (currently hardcoded to 1)
        points_per_polygon = int(lines[0].strip())  # Get the number of points for the polygon from the first line
        polygons = []  # Initialize an empty list to store polygons

        for i in range(num_polygons):
            polygon = []  # Initialize an empty list for the current polygon
            for j in range(1, 1 + points_per_polygon):  # Loop through the number of points
                x, y = map(float, lines[j].strip().split(','))  # Extract x and y coordinates from the line
                polygon.append([x, y])  # Add the point to the current polygon
            polygons.append(polygon)  # Add the current polygon to the list of polygons

    return num_polygons, polygons  # Return the number of polygons and the list of polygons

def write_polygons_to_csv(polygons, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)  # Create a CSV writer object
        writer.writerow(["num", "polygon"])  # Write the header row
        for num, polygon_list in polygons:
            for polygon in polygon_list:
                writer.writerow([num, json.dumps(polygon)])  # Write the polygon data to the CSV file

def process_folder(folder_path):
    all_polygons = []  # Initialize a list to store all polygons from all files
    for filename in os.listdir(folder_path):  # Loop through each file in the folder
        if filename.endswith('.txt'):  # Check if the file is a .txt file
            file_path = os.path.join(folder_path, filename)  # Get the full path of the file
            num, polygons = read_polygon_from_txt(file_path)  # Read polygons from the file
            all_polygons.append((num, polygons))  # Append the polygons to the all_polygons list

    output_file = f"{folder_path}.csv"  # Define the output CSV file name
    write_polygons_to_csv(all_polygons, output_file)  # Write all polygons to the CSV file

def main():
    input_folder = './Test-Data2'  # Set the path to the input folder containing text files
    
    for foldername in os.listdir(input_folder):  # Loop through each folder in the input directory
        folder_path = os.path.join(input_folder, foldername)  # Get the full path of the folder
        if os.path.isdir(folder_path):  # Check if the path is a directory
            process_folder(folder_path)  # Process the folder to read and write polygons

if __name__ == "__main__":
    main()  # Execute the main function
