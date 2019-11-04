#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np

speed = [2.0, 5.0, 10.0, 20.0]
volts = [0.32, 0.80, 1.59, 3.18]
pervolt = volts[-1]/speed[-1]
print(pervolt)
print(1.0/pervolt)
speeds = np.arange(0,np.max(speed)+1)
plt.plot(speeds, speeds*pervolt, 'b-')
plt.plot(speed, volts, 'r+')
plt.savefig('rotationvals.png')
