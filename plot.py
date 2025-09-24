import subprocess
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os
from scipy.interpolate import CubicSpline

EXECUTABLE_NAME = "./kg_sim" 
if os.name == 'nt':
    EXECUTABLE_NAME = "kg_sim.exe"

fig, ax = plt.subplots()
plt.style.use('seaborn-v0_8-darkgrid')

line = None
points = None
time_text = None 
sites = []
sites_smooth = []

ax.set_ylim(-0.1, 1.1)
ax.set_xlabel("Lattice Site")
ax.set_ylabel("Probability")
fig.suptitle("Klein-Gordon Particle Propagation", fontsize=16)

try:
    process = subprocess.Popen(EXECUTABLE_NAME, stdout=subprocess.PIPE, text=True)
except FileNotFoundError:
    print(f"Error: The executable '{EXECUTABLE_NAME}' was not found.")
    print("Please make sure you have compiled the C++ code first.")
    exit()

def init():
    """Initializes the plot before the animation starts."""
    global line, points, sites, sites_smooth, time_text
    
    first_line = process.stdout.readline()
    if not first_line:
        return []
    
    data = [float(val) for val in first_line.strip().split(',')]
    probabilities = data[1:]
    num_sites = len(probabilities)
    sites = np.arange(num_sites)
    sites_smooth = np.linspace(sites.min(), sites.max(), 150)
    
    spline = CubicSpline(sites, probabilities)
    probs_smooth = spline(sites_smooth)
    
    ax.set_xlim(sites.min(), sites.max())

    line, = ax.plot(sites_smooth, probs_smooth, '-', lw=2.5, color='deepskyblue', label='Smoothed Wavepacket')
    points, = ax.plot(sites, probabilities, 'o', color='navy', markersize=5, label='Discrete Probabilities')
    
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes, fontsize=12,
                       verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.legend(loc='upper right')
    
    return line, points, time_text

def update(frame):
    line_data = process.stdout.readline()
    
    if not line_data:
        ani.event_source.stop()
        return ()
        
    data = [float(val) for val in line_data.strip().split(',')]
    current_time = data[0]
    probabilities = data[1:]
    
    spline = CubicSpline(sites, probabilities)
    probs_smooth = spline(sites_smooth)
    
    line.set_ydata(probs_smooth)
    points.set_ydata(probabilities)
    
    time_text.set_text(f'Time = {current_time:.2f} s')
        
    return line, points, time_text

ani = animation.FuncAnimation(fig, update, init_func=init, blit=True, interval=75, repeat=False, cache_frame_data=False)

plt.show()

process.terminate()

