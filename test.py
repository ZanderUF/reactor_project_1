## Test
import numpy as np
import matplotlib.pyplot as plt

x= np.linspace(0,25,1000)

y = ((0.125)*(2.3 + 0.23*x))*(np.exp(-0.125*(2.38*x+0.238*np.power(x,2))))

bounding = ((-.3*x)/15) + 0.3

ee = 2*(np.exp(-(np.power(x,2)))) 

#e = np.power(x,2)
test = 1-(3*np.power(x,2))/2 + (1*np.power(x,3))/3   

test2 = 0.2875 - 0.0539*x - 0.004649*np.power(x,2) + 0.0016* np.power(x,3) - 0.000017*np.power(x,4) - 0.00002*np.power(x,5)+1.031*np.power(10,-6)*np.power(x,6)

actual = (0.05)*(2.3 + .23*x)*np.exp((-2.3*x -.23*np.power(x,2))*0.05)

x4 = np.power(x,4)

#plt.plot(x,y)
#plt.plot(x,bounding)
#plt.plot(x,e)
plt.plot(x,actual, label='actual')
plt.plot(x,test2, label='approx')
plt.axis([0,25,0,1])
#plt.plot(x,x4, label='x^4')
plt.legend()
plt.grid(True)
plt.show()