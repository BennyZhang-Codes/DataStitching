import os

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation



img1 = np.random.rand(100, 100)
img2 = np.random.rand(100, 100)
imgs = [img1, img2]

fig, ax = plt.subplots(1, 1)

images = []
for i in range(2):
    ax.set_title('Frame {}'.format(i))
    im = ax.plot(img1[i, :])
    images.append(im)


# Writer = animation.writers['ffmpeg']
# # Writer = animation.writers['pillow']
# writer = Writer(fps=5, metadata=dict(artist='张金源'), bitrate=500000)
# ani = animation.ArtistAnimation(fig, images, interval=5, repeat_delay=1000)
# print('Saving...')
# ani.save('Loss_30fps-.mp4', writer=writer)

Writer = animation.writers['pillow']
writer = Writer(fps=5, metadata=dict(artist='张金源'), bitrate=500000)
ani = animation.ArtistAnimation(fig, images, interval=500, repeat_delay=0)
print('Saving...')
ani.save('Loss_30fps-.gif', writer=writer)


plt.show()