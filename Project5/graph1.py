import pandas as pd
import matplotlib.pyplot as plt

def plot2(filename, x_label, y_label, title):

    #read file into dataframe
    data_frame = pd.read_csv(filename, delim_whitespace=True, usecols=[0,1,2,3,4], header=None)
    data_frame.columns = ['r', 'E0', 'E1', 'E2', 'E3']
    
    #Make a scatter plot for each object 
    plt.scatter(data_frame['r'], data_frame['E0'], color='blue', marker='o', label='Energy0')
    plt.scatter(data_frame['r'], data_frame['E1'], color='red', marker='o', label='Energy1')
    plt.scatter(data_frame['r'], data_frame['E2'], color='green', marker='o', label='Energy2')
    plt.scatter(data_frame['r'], data_frame['E3'], color='yellow', marker='o', label='Energy3')

    #Label and name 
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title) 

    #Make gridlines 
    plt.grid(True, which='both', ls="--", linewidth=0.5)

    #Add a legend
    plt.legend(loc='upper right')

    # Save the plot
    safe_title = title.replace(" ", "_").replace(":", "").lower()  # Make title safe for filenames
    plt.savefig(f'{safe_title}.png')
    
    # Display the plot
    plt.show()

    plt.close()  # Close the plot to free up memory

    
plot2('r_vs_energy.dat','r (fm)','E (MeV)','Radius vs. Energy')



