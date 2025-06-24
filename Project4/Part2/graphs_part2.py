import pandas as pd
import matplotlib.pyplot as plt

def plot(filename, x_label, y_label, title):

    #read file into dataframe
    data_frame = pd.read_csv(filename, delim_whitespace=True, usecols=[2,3,4,5], header=None)
    data_frame.columns = ['x1', 'y1', 'x2', 'y2']
    
    #Make a scatter plot for each object 
    plt.scatter(data_frame['x1'], data_frame['y1'], color='blue', marker='o', label='Object1')
    plt.scatter(data_frame['x2'], data_frame['y2'], color='red', marker='o', label='Object2')

    #Label and name 
    plt.xlabel('X position')
    plt.ylabel('Y position')
    plt.title(title) 

    #Make gridlines 
    plt.grid(True, which='both', ls="--", linewidth=0.5)

    # Save the plot
    safe_title = title.replace(" ", "_").replace(":", "").lower()  # Make title safe for filenames
    plt.savefig(f'{safe_title}.png')
    
    # Display the plot
    plt.show()
    plt.legend(loc='upper right')
    plt.close()  # Close the plot to free up memory

plot('orbitadvanced2.dat','x','y','Orbit of 2 Objects Around a Primary: 2nd Run') #Change name/number for each run of the program
