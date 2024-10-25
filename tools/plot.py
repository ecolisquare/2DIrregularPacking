import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plot_res(polygons_res, polygons_con_res, container_width, container_height, output_dir):
    fig, ax = plt.subplots()
    ax.set_xlim(0, max(container_height, container_width) + 50)  # Set x-axis limits
    ax.set_ylim(0, max(container_height, container_width) + 50)  # Set y-axis limits

    # Draw a horizontal line at y = container_width
    ax.axhline(y=container_width, color='red', linestyle='--', label=f'y = {container_width}')
    # Draw a vertical line at x = container_height
    ax.axvline(x=container_height, color='red', linestyle='--', label=f'y = {container_height}')

    for polygon in polygons_res:
        polygon = patches.Polygon(polygon, closed=True, fill=True, edgecolor='black', facecolor='blue', alpha=0.5)
        ax.add_patch(polygon)  # Add each filled polygon to the plot
    
    for polygons_con in polygons_con_res:
        polygons_con = patches.Polygon(polygons_con, closed=True, fill=False, edgecolor='green')
        ax.add_patch(polygons_con)  # Add each outline polygon to the plot
    
    plt.title('Polygon Placement Layout')  # Set the title of the plot
    plt.xlabel('Height')  # Label for x-axis
    plt.ylabel('Width')   # Label for y-axis
    plt.savefig(output_dir + 'packing_image.png')  # Save the plot as an image
