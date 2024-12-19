import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

img1 = np.random.randn(100, 100)
img2 = np.random.randn(100, 100)
imgs = [img1, img2]


img_w  = 'E:/pulseq/20241104_ABDL/syn/1p0_Stitched_w.png'
img_wo = 'E:/pulseq/20241104_ABDL/syn/1p0_Stitched_wo.png'
img1 = plt.imread(img_w)
img2 = plt.imread(img_wo)
imgs = [img1, img2]
titles = ['w/ Synchronization ', 'w/o Synchronization']

fig, ax = plt.subplots()
l = ax.set_title(titles[0])
im1 = ax.imshow(img1, vmin=-1, vmax=1, cmap='gray')
ax.axis('off')
def update(frame):
    print(frame)

    l.set_text(titles[frame])
    im1.set_data(imgs[frame])
    return [im1, l]

# 创建动画
ani = FuncAnimation(fig, update, frames=2, interval=500, blit=False, repeat=True, repeat_delay=0)

# # 显示动画
plt.show()