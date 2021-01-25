import time
import numpy as np
import matplotlib.pyplot as plt

# https://stackoverflow.com/questions/51922480/javascript-error-ipython-is-not-defined-in-jupyterlab
# put '%matplotlib widget' in notebook to update image interactively
# image widget for jupyter-lab
# pip install ipympl
# conda install nodejs
# jupyter labextension install @jupyter-widgets/jupyterlab-manager
# jupyter labextension install jupyter-matplotlib
class ImageViewer():
    def __init__(self, size, vmin=0, vmax=1):
        self.size = size
        self.vmin = vmin
        self.vmax = vmax
    def init(self):
        self.fig, self.ax = plt.subplots(1,1)
        self.im = self.ax.imshow(np.zeros(self.size), cmap='gist_gray_r', vmin=self.vmin, vmax=self.vmax)
    def show(self):
        return self.im
    def update(self, data):
        self.im.set_data(data.reshape(self.size))
        self.fig.canvas.draw()
        time.sleep(1e-2)