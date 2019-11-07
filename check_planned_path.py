import numpy as np
import matplotlib.pyplot as plt


pathnew = np.load('./planned_path.npy')
dim = 4
plt.axis([pathnew[0]-10, pathnew[-4]+20, -100, 100],)
plt.plot(pathnew[2*0 :  : dim], pathnew[2*0+1 :  : dim], 'bo')
plt.plot(pathnew[2*1 :  : dim], pathnew[2*1+1 :  : dim], 'rs')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Cozmo Trajectory Generation!')
plt.show()
































