#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd
import matplotlib.pyplot as plt

def loglog_plot(filename, x_label, y_label, title):
    
    #read file into dataframe
    data_frame = pd.read_csv(filename, delim_whitespace=True, header=None)
    data_frame.columns = [x_label, y_label]

    #Take the absolute value of the x,y axis data to avoid taking the log of negative values 
    data_frame[x_label] = data_frame[x_label].abs()
    data_frame[y_label] = data_frame[y_label].abs()
    
    #Make a scatter plot
    plt.scatter(data_frame[x_label], data_frame[y_label], color='blue', marker='o')

    #Set to log-log scale 
    plt.xscale('log')
    plt.yscale('log')

    #Label and name 
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title) 

    #Make gridlines 
    plt.grid(True, which='both', ls="--", linewidth=0.5)
    plt.show()


# In[3]:


#File1 
loglog_plot('derivative1.dat', 'h (delta_x)', '|f_exact_d1(x0)-f_central_d1(x0)|', 'First Derivative: 2nd Order Central Difference x0 =0.7')

#File2
loglog_plot('derivative2.dat', 'h (delta_x)', '|f_exact_d2(x0)-f_central_d2(x0)|', 'Second Derivative: 4th Order Central Difference x0 =0.7')


# In[ ]:
