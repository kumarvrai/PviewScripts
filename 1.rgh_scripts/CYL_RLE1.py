"""
Created on Sun Feb 17 17:53:23 2019

@author: dpastrana
"""
import sys
import math as m
import matplotlib.pylab as plt
import numpy as np

plt.close('all')


def car_pol(x, y):
    r = m.sqrt(x**2.0 + y**2.0)
    ang = m.atan2(y,x)
    return r, ang




def pol_car(r, ang):
    x = r * m.cos(ang * m.pi / 180.0)
    y = r * m.sin(ang * m.pi / 180.0)
    return x, y


def unwrap_cy(x, y, z=0.0):
    s = m.sqrt(x**2.0+y**2.0)
    t = m.atan2(y,x)
    z = z
    return s, t, z

def wrap_cy(s, t, z=0.0):
    x = s * m.cos(t)
    y = s * m.sin(t)
    z = z

    return x, y, z


def arc_length(r, t):
    arc = m.pi*2.0*r * (t/360)
    return arc




geofile_in = './cylinder.geo.dat'
geofile_out= './cylinder_roughness.geo.dat'
coords = 'cylinder_roughness.coords'

f = open(geofile_in, 'r')
g = open(geofile_out, 'w')
o = open(coords,'w')



lineB = True
inCoords = False

#PARAMETERS
D = 0.08
h = 0.008
#D = 0.1
#h = 0.05
H = 1.0
fact = 10.0 # 50


b_radio = 0.555

Roughness = -1
delta_z = 0.2 #distance in z between spheres
delta_x = 0.15  #distance in x between spheres
initial_pos = 17.45760312372209

x1 = unwrap_cy(pol_car(0.5, 180)[0], pol_car(0.5,   180)[1])[1]
x2 = unwrap_cy(pol_car(0.5, -180)[0], pol_car(0.5,  -180)[1])[1]
n_lines = 16 

Xpos = []
for i in range(1, n_lines+1):
    Xpos.append(x1-i*delta_x)
    Xpos.append(x2+i*delta_x)
    

Zpos = []
for i in range(0, n_lines):
    if i%2 == 0:
        Zpos.append([x for x in np.linspace(0, 1, 6)])
    else:
        Zpos.append([x for x in np.linspace(0.1, 0.9, 5)])



center_y =  (h**2.0 + (0.5*D)**2.0)/(2.0*h) - h
rad = center_y+h


center = []
for i in range(0, len(Zpos)):
    for x in Xpos[2*i:2*(i+1)]:
        for z in Zpos[i]:
            center.append((x,z))
            #print(center[i])

Roughness = len(center)

ii = 0
done = False

xx = []
yy = []
tt = []
ss = []

print('INI...')

for line in f:
    if inCoords:
        ii = ii + 1
        if 'END' in line:
            print('END WRITING...')
            print(line)


            inCoords = False
            lineB = True
            done = True

        else:
            line_elem = list(map(float, line.split()))
            x = line_elem[1]
            y = line_elem[2]
            z = line_elem[3]

            if m.sqrt( line_elem[1]**2.0 + line_elem[2]**2.0 ) <= b_radio:
                t = unwrap_cy(x, y)[1] #+ m.pi
                s = unwrap_cy(x, y)[0] - 0.5


                for n in range(0,Roughness):
                    rad_loc = m.sqrt((t - center[n][0])**2.0 + (z - center[n][1])**2.0)
#                y_sphere = s

                    if rad_loc <=D*0.5:
                        y_sphere = m.exp(-s*fact)*(-((center_y - m.sqrt(rad**2.0 - (t-center[n][0])**2.0 - (z-center[n][1])**2.0)) - s)) + (1.0-m.exp(-s*fact))*s
                        s = ((H-s)/H)*y_sphere + (1.0 - (H-s)/H)*s


                g.write("%i %.17f %.17f %.17f\n"%(ii,wrap_cy(s+0.5, t)[0],wrap_cy(s+0.5, t)[1],z))
                o.write("%i %.17f %.17f %.17f\n"%(ii,wrap_cy(s+0.5, t)[0],wrap_cy(s+0.5, t)[1],z))

                if z == 0.0:
                    xx.append(t)
                    yy.append(s)
                    tt.append(wrap_cy(s+0.5, t)[0])
                    ss.append(wrap_cy(s+0.5, t)[1])

            else:
                g.write("%i %.17f %.17f %.17f\n"%(ii,line_elem[1],line_elem[2],z))
                o.write("%i %.17f %.17f %.17f\n"%(ii,line_elem[1],line_elem[2],z))


    if lineB:
        g.write(line)
        if 'COORD' in line:
            if done == False:
                inCoords = True
                lineB = False
                print('START WRITING...')
                print(line)



f.close()
g.close()
o.close()


plt.figure()
plt.plot(xx[::1], yy[::1], marker=".", linewidth=0,  markersize=1)
plt.xlabel('$theta$')
plt.ylabel('$R-R_0$')
plt.title('Un-wrapped cylinder view at $z=0.5$')
plt.savefig('1.eps', format='eps')
plt.show()
#
fig, axs = plt.subplots()
plt.plot(tt[::1], ss[::1], marker=".", linewidth=0,  markersize=1)
plt.plot(pol_car(0.5,175)[0], pol_car(0.5, 175)[1], '*')
plt.plot(pol_car(0.5,-175)[0], pol_car(0.5, -175)[1], '*')
axs.axis('equal')
plt.xlabel('$x$')
plt.ylabel('$y$')
plt.title('Wrapped cylinder view at $z=0.5$')
plt.savefig('2.eps', format='eps')
plt.show()
