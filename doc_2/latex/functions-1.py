import numpy as np
import matplotlib.pyplot as plt
import astrolyze.functions.astro_functions as astFunc

frequency_range = np.arange(1e4, 1e7, 1e4)
temperature_1 = 6000
temperature_2 = 12000  # Kelvin
blackbody_1 = astFunc.black_body(frequency_range, temperature_1)
blackbody_2 = astFunc.black_body(frequency_range, temperature_2)
figure = plt.figure()
axis = figure.add_subplot(111)
pl = axis.loglog(frequency_range, blackbody_1, label='T = 6000 K')
pl = axis.loglog(frequency_range, blackbody_2, label='T = 12000 K')
pl = axis.legend()
plt.savefig('black_body.eps')