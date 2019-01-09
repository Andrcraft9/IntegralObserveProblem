import numpy as np
import matplotlib.pyplot as plt
import sys

f1 = open("u_c.txt", "r")
TN1 = int(f1.readline())
u1 = np.zeros(TN1)
for i in range(TN1):
    u1[i] = (float(f1.readline()))
f1.close()
    
f2 = open("u_c_solution.txt", "r")
TN2 = int(f2.readline())
u2 = np.zeros(TN2)
for i in range(TN2):
    u2[i] = (float(f2.readline()))
f2.close()

plt.title("Compare u_c")
plt.plot(u1, label="Numerical solution")
plt.plot(u2, label="Analytical solution")
plt.xlabel("j")
plt.ylabel("u_c(j)")
plt.grid()
plt.legend()
plt.show()
