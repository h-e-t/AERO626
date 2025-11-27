import numpy as np
import matplotlib.pyplot as plt

p0 = 101325
L = 0.0065
T0 = 288.15
g = 9.80655
R = 8.314
M = 0.0289644

altitudes = np.array([_ for _ in np.linspace(10000, 11000, 100)])
pressures = p0 * (1 - L * altitudes / T0) ** (g * M / (R * L))
plt.plot(altitudes, pressures/1000,color='black',linewidth=3)
plt.xlabel("Altitude (m)")
plt.ylabel("Pressure (kPa)")
plt.show()
