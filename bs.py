import numpy as np
import matplotlib.pyplot as plt

bs = np.loadtxt('groundstate.txt')
plt.title('Bright Soliton')
plt.xlim([-64,64])
plt.plot(bs[:,0], bs[:,1],linewidth=2,label=r"$|\psi(x)|^{2}$")
plt.xlabel("x")
plt.legend(loc=1, prop={'size': 10})
plt.savefig("brightsoliton",format='jpg')
