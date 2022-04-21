#!/usr/bin/env python3.3

import numpy as np
import matplotlib.pyplot as plt
import os


number = 5000
step = 10
data_path = "pisoTest"
pics_path = data_path + "/" + "pics"

if os.path.exists(pics_path) == False:
    os.mkdir(pics_path)

def mplot(n, u0, v0, p0):
    plt.subplot(3, 1, 1)
    plt.pcolor(u0)
    plt.colorbar()
    plt.subplot(3, 1, 2)
    plt.pcolor(v0)
    plt.colorbar()
    plt.subplot(3, 1, 3)
    plt.pcolor(p0)
    plt.colorbar()
    plt.savefig("./" + pics_path + "/" + str(n) + ".png")
    # wait = input("PRESS ENTER TO CONTINUE.")
    plt.close()

for item in range(0, number):
    print("Current: ", item)
    mplot(item,
          np.load("./" + data_path + "/u0_" + str(item * step) + ".npy"),
          np.load("./" + data_path + "/v0_" + str(item * step) + ".npy"),
          np.load("./" + data_path + "/p0_" + str(item * step) + ".npy")
         )
