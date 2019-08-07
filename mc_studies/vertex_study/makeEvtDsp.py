import matplotlib.pyplot as plt
from matplotlib  import cm

_x, _y, _z, _e = [],[],[],[]
_sub_x, _sub_y, _sub_z, _sub_e = [],[],[],[]

with open('ides_for_viewing.txt') as f:
    while True:
        line = f.readline()
        if not line:
            break
        line_vec = line.split()
        _x.append(float(line_vec[0]))
        _y.append(float(line_vec[1]))
        _z.append(float(line_vec[2]))
        _e.append(float(line_vec[3]))
 
with open('ides_for_viewing_sub.txt') as f:
    while True:
        line = f.readline()
        if not line:
            break
        line_vec = line.split()
        print float(line_vec[3])
        if float(line_vec[3]) <0.001:
            continue
        _sub_x.append(float(line_vec[0]))
        _sub_y.append(float(line_vec[1]))
        _sub_z.append(float(line_vec[2]))
        _sub_e.append(float(line_vec[3]))

fig0 = plt.figure(0)
ax0 = fig0.add_subplot(111)
s0 = plt.scatter(_z, _x, c=_e, cmap=cm.gnuplot, vmin=0, vmax=0.5)
fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
s1 = plt.scatter(_sub_z, _sub_x, c=_sub_e, cmap=cm.gnuplot, vmin=0, vmax=0.5)
fig0.colorbar(s0)
fig1.colorbar(s1)
ax0.set_xlim([0,100])
ax0.set_ylim([0,47])
ax1.set_xlim([0,100])
ax1.set_ylim([0,47])
plt.show()