import numpy as np
import matplotlib.pyplot as plt

rgroup = np.linspace(0,10000, 100,dtype = float)

r_location = [2,3]
for i in r_location:
    exp = rgroup**i
    plt.plot(rgroup, exp, label = i)
    

plt.legend(title = "Number of R-group Location")       
plt.xlabel("Number of R-groups")
plt.ylabel("Total Number of Molecules Possible")
plt.title("Exponential Behavior of Adding R-group Locations")
 
plt.show()