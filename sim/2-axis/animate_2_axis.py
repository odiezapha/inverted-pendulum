import matplotlib.pyplot as plt
import numpy
import os

import pendulum_control_2_axis

dt = 0.01
x = [0.003, 0.001, 0., 0., 0., 0., 0., 0., 0.]
dur = 5
pendulum = pendulum_control_2_axis.Pendulum(
    dt, x, dur,
)
data = pendulum.integrate()

fig = plt.figure(0)

"""
x1_timeline = plt.subplot2grid(
    (60, 12),
    (0, 0),
    colspan=12,
    rowspan=7
)
x1_timeline.axis([
    0,
    dur,
    min(data[:,12])*1.1,
    max(data[:,12])*1.1+.1,
])
x1_timeline.set_ylabel('Theta_x')
x1_timeline.yaxis.set_label_position("right")
x1_timeline.plot(data[:,0],data[:,12],'c-')

x2_timeline = plt.subplot2grid(
    (60, 12),
    (10, 0),
    colspan=12,
    rowspan=7
)
x2_timeline.axis([
    0,
    dur,
    min(data[:,13])*1.1,
    max(data[:,13])*1.1+.1,
])
x2_timeline.set_ylabel('Theta_y')
x2_timeline.yaxis.set_label_position("right")
x2_timeline.plot(data[:,0],data[:,13],'c-')

x3_timeline = plt.subplot2grid(
    (60, 12),
    (20, 0),
    colspan=12,
    rowspan=7
)
x3_timeline.axis([
    0,
    dur,
    min(data[:,10])*1.1,
    max(data[:,10])*1.1+.1,
])
x3_timeline.set_ylabel('V_x')
x3_timeline.yaxis.set_label_position("right")
x3_timeline.plot(data[:,0],data[:,10],'r')

x4_timeline = plt.subplot2grid(
    (60, 12),
    (30, 0),
    colspan=12,
    rowspan=7
)
x4_timeline.axis([
    0,
    dur,
    min(data[:,11])*1.1,
    max(data[:,11])*1.1+.1,
])
x4_timeline.set_ylabel('V_y')
x4_timeline.yaxis.set_label_position("right")
x4_timeline.plot(data[:,0],data[:,11],'r')

x5_timeline = plt.subplot2grid(
    (60, 12),
    (40, 0),
    colspan=12,
    rowspan=7
)
x5_timeline.axis([
    0,
    dur,
    min(data[:,7])*1.1,
    max(data[:,7])*1.1+.1,
])
x5_timeline.set_ylabel('W_x')
x5_timeline.yaxis.set_label_position("right")
x5_timeline.plot(data[:,0],data[:,7],'b-')

x6_timeline = plt.subplot2grid(
    (60, 12),
    (50, 0),
    colspan=12,
    rowspan=7
)
x6_timeline.axis([
    0,
    dur,
    min(data[:,8])*1.1,
    max(data[:,8])*1.1+.1,
])
x6_timeline.set_xlabel('time (s)')
x6_timeline.yaxis.set_label_position("right")
x6_timeline.set_ylabel('W_y')
x6_timeline.plot(data[:,0],data[:,8],'b-')

plt.show()

"""
x1_timeline = plt.subplot2grid(
    (21, 12),
    (10, 0),
    colspan=12,
    rowspan=3
)
x1_timeline.axis([
    0,
    dur,
    min(data[:,12])*1.1,
    max(data[:,12])*1.1+.1,
    ])
x1_timeline.set_ylabel('Theta_x')
x1_timeline.yaxis.set_label_position("right")
x1_timeline.plot(data[:,0],data[:,12],'g-')

x2_timeline = plt.subplot2grid(
    (21, 12),
    (15, 0),
    colspan=12,
    rowspan=3
)
x2_timeline.axis([
    0,
    dur,
    min(data[:,13])*1.1,
    max(data[:,13])*1.1+.1,
    ])
x2_timeline.set_ylabel('Theta_y')
x2_timeline.set_xlabel('Seconds')
x2_timeline.yaxis.set_label_position("right")
x2_timeline.plot(data[:,0],data[:,13],'g-')

cart_plot = plt.subplot2grid(
    (21, 12),
    (0,3),
    rowspan=8,
    colspan=6
)
cart_plot.axes.get_yaxis().set_visible(False)
cart_plot.axes.get_xaxis().set_visible(False)

def draw_point(point):
    cart_plot.cla()
    cart_plot.axis([-1,1,-0.5,0.5])
    cart_plot.plot([0.01+2*point[12],2*point[12]],[0.01+2*point[13],2*point[13]],'g-', lw=8)

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
