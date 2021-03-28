import matplotlib.pyplot as plt
from .slider_plot import *

def mpl_animation(pos, vel):
    
    plt.close()
    fig = Interactive_3D_Scatter(steps, pos, xlim=[0, 1/sigma], ylim=[0, 1/sigma], zlim=[0, 1/sigma], markersize=150, stepsize=500, subplot=131)
    ax1 = fig.add_subplot(132)
    ax1.plot(range(steps),E_tot,lw=2)
    ax1.set_xlim([0,steps-1])
    ax1.set_ylim([-0.15,0.15])
    ax1.set_title("Total energy.")
    ax1.grid(True)

    ax2 = fig.add_subplot(133)
    ax2.plot(range(steps), rel_dist_0, lw=2)
    ax2.plot(range(steps), rel_dist_1, lw=2, ls='dashed')
    ax2.set_title("Relative Distance of 2 particles.")
    ax2.grid(True)
    plt.show()