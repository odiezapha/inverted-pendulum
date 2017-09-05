import os
from math import sin, cos, pi

import matplotlib.pyplot as plt

import pendulum_control_1_axis

dt = 0.01
x = [0.1, 0., 0, 0.]
dur = 15
pendulum = pendulum_control_1_axis.Pendulum(
    dt,
    x,
    dur,
)
data = pendulum.integrate()

fig = plt.figure(0)

"""
x1_timeline = plt.subplot2grid(
    (40, 12),
    (0, 0),
    colspan=12,
    rowspan=8
)
x1_timeline.axis([
    0,
    dur,
    min(data[:,1])*1.1,
    max(data[:,1])*1.1+.1,
])
x1_timeline.set_ylabel('Theta_p')
x1_timeline.plot(data[:,0],data[:,1],'r-')

x2_timeline = plt.subplot2grid(
    (40, 12),
    (10, 0),
    colspan=12,
    rowspan=8
)
x2_timeline.axis([
    0,
    dur,
    min(data[:,2])*1.1,
    max(data[:,2])*1.1+.1,
])
x2_timeline.set_ylabel('W_p')
x2_timeline.plot(data[:,0],data[:,2],'c-')

x3_timeline = plt.subplot2grid(
    (40, 12),
    (20, 0),
    colspan=12,
    rowspan=8
)
x3_timeline.axis([
    0,
    dur,
    min(data[:,5])*1.1,
    max(data[:,5])*1.1+.1,
])
x3_timeline.set_ylabel('V')
x3_timeline.plot(data[:,0],data[:,5],'k-')

x4_timeline = plt.subplot2grid(
    (40, 12),
    (30, 0),
    colspan=12,
    rowspan=8
)
x4_timeline.axis([
    0,
    dur,
    min(data[:,4])*1.1,
    max(data[:,4])*1.1+.1,
])
x4_timeline.set_xlabel('time (s)')
x4_timeline.set_ylabel('W_w')
x4_timeline.plot(data[:,0],data[:,4],'b-')

plt.show()

"""
x1_timeline = plt.subplot2grid(
    (18, 12),
    (11, 0),
    colspan=12,
    rowspan=5
)
x1_timeline.axis([
    0,
    dur,
    min(data[:,1])*1.1,
    max(data[:,1])*1.1+.1,
])
x1_timeline.set_ylabel('Theta (Rad)')
x1_timeline.set_xlabel('Seconds')
x1_timeline.yaxis.set_label_position("right")
x1_timeline.plot(data[:,0],data[:,1],'r-')

cart_plot = plt.subplot2grid(
    (18, 12),
    (0,3),
    rowspan=10,
    colspan=6
)
cart_plot.axes.get_yaxis().set_visible(False)
cart_plot.axes.get_xaxis().set_visible(False)

time_bar, = x1_timeline.plot([0,0], [10, -10], lw=2)

def draw_point(point):
    time_bar.set_xdata([t, t])
    cart_plot.cla()
    cart_plot.axis([-0.7,0.7,0,0.7])
    cart_plot.plot([.4*sin(point[1]),(.4*sin(point[1])+.1*sin(point[3]))],[.4*cos(point[1]),(.4*cos(point[1])+.1*cos(point[3]))],'r-',lw=2)
    cart_plot.plot([0,.4*sin(point[1])],[0,.4*cos(point[1])],'g-', lw=3)

t = 0
fps = 25.
frame_number = 1
for point in data:
    if point[0] >= t + 1./fps or not t:
        draw_point(point)
        t = point[0]
        fig.savefig('img/_tmp%03d.png' % frame_number)
        frame_number += 1

print os.system("ffmpeg -framerate 25 -i img/_tmp%03d.png  -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4")
