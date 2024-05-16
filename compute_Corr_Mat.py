# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 12:45:49 2024

@author: nivfe
"""

import numpy as np
import pandas as pd
import hazard as hz
import correlation_3 as cor
import argparse
import matplotlib.pyplot as plt
import gc
from scipy.stats import norm, lognorm, genextreme
from scipy.optimize import curve_fit




print('Hello')


#Load Sites

#distance_matrix=np.linspace(0,150)

meshFileName="sitemesh_66.csv"

mesh_df=pd.read_csv(meshFileName,header=1)
distance_matrix=cor.distances(mesh_df)



#Val_mat=correlation_matrix =np.load('Correlation_mat.npy')
#Vanmarcke=cor.compute_correlation_matrix_Vanmarcke(distance_matrix)

#Boore_2003=cor.compute_correlation_matrix_Boore_2003(distance_matrix)


Goda_Atk_2009=cor.compute_correlation_matrix_Goda_Atkinson_2009_PGA(distance_matrix)
#J_Baker_2009=cor.compute_correlation_matrix_Jayaram_Baker_2009_PGA(distance_matrix)


#Goda_Hong_PGV=cor.compute_correlation_matrix_Goda_Hong_2008_PGV(distance_matrix)
#Wang=cor.compute_correlation_matrix_Wang_Takada_2005_PGV(distance_matrix)
plt.style.use('ggplot')
legends=['Boore 2003 PGA a=0.774 b=0.5 c=1'
         ,'Goda Atkinson 2009 PGA a=0.095 b=0.336 c=2.6'
         ,'Jayaram Baker 2010 PGA a=0.33 b=1 c=1'
         ,'GodaHong 2008 PGV a=0.58 b=0.5 c=1'
         ,'Wang Takada 2005 PGV a=0.048 b=1 c=1'
         ]
#fig, axs = plt.subplots(2, 3, figsize=(18, 10))

# Plot each correlation matrix
# axs[0, 0].imshow(Boore_2003, cmap='inferno', interpolation='nearest')
# axs[0, 0].set_title(legends[0])
# axs[0, 0].set_xticks([])
# axs[0, 0].set_yticks([])
# axs[0, 0].set_xticklabels([])
# axs[0, 0].set_yticklabels([])

# axs[0, 1].imshow(Goda_Atk_2009, cmap='inferno', interpolation='nearest')
# axs[0, 1].set_title(legends[1])
# axs[0, 1].set_xticks([])
# axs[0, 1].set_yticks([])
# axs[0, 1].set_xticklabels([])
# axs[0, 1].set_yticklabels([])
# axs[0, 2].imshow(J_Baker_2009, cmap='inferno', interpolation='nearest')
# axs[0, 2].set_title(legends[2])
# axs[0, 2].set_xticks([])
# axs[0, 2].set_yticks([])
# axs[0, 2].set_xticklabels([])
# axs[0, 2].set_yticklabels([])

# axs[1, 0].imshow(Goda_Hong_PGV, cmap='inferno', interpolation='nearest')
# axs[1, 0].set_title(legends[3])
# axs[1, 0].set_xticks([])
# axs[1, 0].set_yticks([])
# axs[1, 0].set_xticklabels([])
# axs[1, 0].set_yticklabels([])
# axs[1, 1].imshow(Wang, cmap='inferno', interpolation='nearest')
# axs[1, 1].set_title(legends[4])
# axs[1, 1].set_xticks([])
# axs[1, 1].set_yticks([])
# axs[1, 1].set_xticklabels([])
# axs[1, 1].set_yticklabels([])
# axs[1, 2].imshow(np.zeros_like(Boore_2003), cmap='inferno', interpolation='nearest')
# axs[1, 2].set_title('No Correlation')
# axs[1, 2].set_xticks([])
# axs[1, 2].set_yticks([])
# axs[1, 2].set_xticklabels([])
# axs[1, 2].set_yticklabels([])

plt.figure(figsize=(12,6.75))
#plt.plot(distance_matrix,Boore_2003, marker='o', linestyle='-')
plt.plot(distance_matrix,Goda_Atk_2009, marker='o', linestyle='-')
#plt.plot(distance_matrix,J_Baker_2009, marker='o', linestyle='-')
#plt.plot(distance_matrix,Goda_Hong_PGV, marker='.', linestyle=':')
#plt.plot(distance_matrix,Wang, marker='.', linestyle=':')
plt.title('Spatial Correlation models', fontsize=28, fontweight='bold')
plt.xlabel('Distance [km]', fontsize=22, fontweight='bold')
plt.ylabel('Correlation [-]', fontsize=22, fontweight='bold')
plt.legend(legends,loc='upper right')
plt.legend(legends,fontsize=14)
plt.savefig('corr_models.png', dpi=100)