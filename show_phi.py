import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import sys

f1 = open("phi.txt", "r")
M1 = int(f1.readline())
TN1 = int(f1.readline())
u1 = np.zeros((M1, TN1))
for j in range(TN1):
    for i in range(M1):
        u1[i, j] = (float(f1.readline()))
f1.close()

f2 = open("phi_solution.txt", "r")
M2 = int(f2.readline())
TN2 = int(f2.readline())
u2 = np.zeros((M2, TN2))
for j in range(TN2):
    for i in range(M2):
        u2[i, j] = (float(f2.readline()))
f2.close()

x = np.linspace(0, 101, TN1)
y = np.linspace(0, 101, M1)
X, Y = np.meshgrid(x, y)

ax = plt.axes(projection='3d')

plt.xlabel("t")
plt.ylabel("x")
#ax.contour3D(X, Y, u1, 50, cmap='binary')
ax.plot_surface(X, Y, u1)
ax.plot_surface(X, Y, u2)

plt.show()

'''
plt.title("Compare u_c")
plt.plot(u1, label="Numerical solution")
#plt.plot(u2, label="Analytical solution")
plt.xlabel("j")
plt.ylabel("u_c(j)")
plt.grid()
plt.legend()
plt.show()
'''