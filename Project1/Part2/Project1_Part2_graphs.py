#!/usr/bin/env python
# coding: utf-8

# In[3]:


import pandas as pd
import matplotlib.pyplot as plt

def plot(filename, x_label, y_label, title):
    
    #read file into dataframe
    data_frame = pd.read_csv(filename, delim_whitespace=True, header=None)
    data_frame.columns = [x_label, y_label]
    
    #Make a scatter plot
    plt.scatter(data_frame[x_label], data_frame[y_label], color='blue', marker='o')

    #Label and name 
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title) 

    #Make gridlines 
    plt.grid(True, which='both', ls="--", linewidth=0.5)

    # Save the plot
    safe_title = title.replace(" ", "_").replace(":", "").lower()  # Make title safe for filenames
    plt.savefig(f'{safe_title}.png')
    
    # Display the plot
    plt.show()
    plt.close()  # Close the plot to free up memory

    

#File1 
plot('LJ_original_function.dat', 'r(angstroms)', 'V(kJ/mole)', 'Leonard Jones Potential')
plot('LJ_derivative1.dat', 'r(angstroms)', 'V (KJ/mole/angstrom)', 'Derivative 1: 2nd Order Central Difference')
plot('LJ_derivative2.dat', 'r(angstroms)', 'V (KJ/mole/angstrom^2)', 'Derivative 2: 4th Order Central Difference')


# In[ ]:





# In[ ]:




