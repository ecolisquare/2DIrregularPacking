# arg.py
import argparse  # Import the argparse module for command-line argument parsing

def get_args():
    # Initialize the argument parser with a description for the script
    parser = argparse.ArgumentParser(description="Parameters for layout arrangement")

    # Add an argument for the container width, which is a required parameter
    parser.add_argument('--container_width', type=float, required=True, 
                        help="Width of the container in millimeters (mm)")

    # Add an optional argument for the gap between parts, with a default value of 15 mm
    parser.add_argument('--part_gap', type=float, default=15, 
                        help="Gap distance between parts, default is 15mm")

    # Add an optional argument for the gap between the parts and the board, with a default value of 10 mm
    parser.add_argument('--board_gap', type=float, default=10, 
                        help="Distance between parts and the board, default is 10mm")

    # Add an argument to select the packing method with a default value of 'MBR'
    parser.add_argument('--packing_method', type=str, default='MBR', choices=['ConvexHull', 'MBR'], 
                        help="Choice of packing method: either 'ConvexHull' or 'MBR', default is MBR")

    # Add a flag to indicate whether to group polygons; if specified, this option will be set to True
    parser.add_argument('--grouping', action='store_true', default=False,
                        help="Indicate whether to group polygons; this option enables grouping")

    # Add an argument for the input path of the data to be processed, which is required
    parser.add_argument('--input_path', type=str, required=True, 
                        help="Path to the input data for arrangement; must be a valid file path")
    
    # Add an argument for the output path where the results will be saved, which is also required
    parser.add_argument('--output_path', type=str, required=True, 
                        help="Path to save the arrangement results; must be a valid file path")

    # Parse the command-line arguments and return the resulting object
    return parser.parse_args()
