# helper functions to streamline notebooks when the matplotlib stuff isn't relevant
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import rc, animation, cm
import numpy as np
from scipy.integrate import quad
from IPython.display import HTML, Image

#change the color/style parameters to your liking
rc("figure",figsize=(4,4))
sqrt2pi = np.sqrt(2*np.pi)
def unot(x):
    return np.exp(-0.5*np.power(2*x,2)/sqrt2pi)
def simple_advection_figure(speed=2):
    x = np.linspace(-8,8,1600)
    y1 = unot(x)
    y2 = unot(x-speed)
    fig, ax = plt.subplots()
    ax.set_title("Simple Linear Advection")
    ax.set_xlim(-1, 4)
    ax.set_ylim(0,2)
    ax.plot(x, y1,label="t=0")
    ax.plot(x, y2,label="t=1")
    ax.legend()
    return fig, ax

def simple_characteristic_figure(speed=0.5):
    x = np.linspace(-2, 2, 400)
    xnot = 0
    fig, ax = plt.subplots()
    ax.set_title("Characteristic Lines for Linear Advection")
    ax.set_xlim(0 ,2)
    ax.set_ylim(0, 1)
    while xnot < 2:
        ax.plot(x, 1/speed * (x-xnot), color=cm.tab10.colors[0])
        xnot += 0.5
    return fig, ax

if __name__=="__main__":
    fig, ax = simple_advection_figure()
    fig.savefig("simple_advection_figure.png")
    fig, ax = simple_characteristic_figure()
    fig.savefig("simple_characteristic_figure.png")