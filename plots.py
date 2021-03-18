# -*- coding: utf-8 -*-
"""
Created on Thu Mar  4 02:36:27 2021

@author: Elif
"""

import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np


import pandas as p
import os
path = os.path.join(os.getcwd(), 'C:\\Users\\Elif\\Desktop\\reaction_full.csv')
df = p.read_csv(path)

x = df['pairPackage_id']
y = df['reactionEnergy']

points= plt.scatter(x, y,cmap="jet", lw=0);
plt.colorbar(points)
plt.xscale('log', basex=2)
plt.xticks(x)