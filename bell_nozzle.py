import math
import numpy as np
import matplotlib.pyplot as plt

from matplotlib.patches import Arc
from bisect import bisect_left


"""
 Implemented from the following technical notes 
 The thrust optimised parabolic nozzle
 http://www.aspirespace.org.uk/downloads/Thrust%20optimised%20parabolic%20nozzle.pdf
 .................................................................

The radius of the nozzle exit: 
Re = √ε * Rt							[Eqn. 2]
and nozzle length 
LN = 0.8 ((√∈−1) * Rt )/ tan(15)		[Eqn. 3]
.................................................................
For the throat entrant section:
x = 1.5 Rt cosθ
y = 1.5 Rt sinθ + 1.5 Rt + Rt			[Eqn. 4]
where: −135 ≤ θ ≤ −90
(The initial angle isn’t defined and is up to the
combustion chamber designer, -135 degrees is typical.)
.................................................................
For the throat exit section:
x = 0.382 Rt cosθ
y = 0.382 Rt sinθ + 0.382 Rt + Rt		[Eqn. 5]
where: −90 ≤ θ ≤ (θn − 90)
.................................................................
The bell is a quadratic Bézier curve, which has equations:
x(t) = (1 − t)^2 * Nx + 2(1 − t)t * Qx + t^2 * Ex, 0≤t≤1
y(t) = (1 − t)^2 * Ny + 2(1 − t)t * Qy + t^2 * Ey, 0≤t≤1 [Eqn. 6]
.................................................................
Selecting equally spaced divisions between 0 and 1 produces 
the points described earlier in the graphical method, 
for example 0.25, 0.5, and 0.75.
.................................................................
Equations 6 are defined by points N, Q, and E (see the graphical method 
earlier for the locations of these points).

Point N is defined by equations 5 setting the angle to (θn – 90).
Nx = 0.382 Rt cos(θn – 90)
Ny = 0.382 Rt sin(θn – 90) + 0.382 Rt + Rt
.................................................................
Coordinate Ex is defined by equation 3, and coordinate Ey is defined by equation 2.
Ex = 0.8*(((√ε−1)-1)*Rt)/(tan(15)) # degrees in rad
Ey = √ε * Rt
.................................................................
Point Q is the intersection of the lines: ⃗⃗⃗⃗⃗⃗
NQ = m1 x + C1 and: ⃗⃗⃗⃗⃗
QE = m2 x + C2 			[Eqn. 7]

where: gradient 
m1 = tan(θn ) , m2 = tan(θe )	[Eqn. 8]

and: intercept 
C1 = Ny − m1 Nx
C2 = Ey − m2 Ex		[Eqn. 9]
.................................................................
The intersection of these two lines (at point Q) is given by:
Qx = (C2 − C1 ) /(m1 − m2 )
Qy = (m1 C2 − m2 C1 ) / (m1 − m2 ) [Eqn. 10]
.................................................................
"""
	
# sp.heat, area_ratio, throat_radius, length percentage, 
def bell_nozzle(k, aratio, Rt, l_percent):
	# upto the nozzle designer, usually -135
	entrant_angle  	= -135
	ea_radian 		= math.radians(entrant_angle)

	# nozzle length percntage
	if l_percent == 60:		Lnp = 0.6
	elif l_percent == 80:	Lnp = 0.8
	elif l_percent == 90:	Lnp = 0.9	
	else:					Lnp = 0.8
	# find wall angles (theta_n, theta_e) for given aratio (ar)		
	angles = find_wall_angles(aratio, throat_radius, l_percent)
	# wall angles
	nozzle_length = angles[0]; theta_n = angles[1]; theta_e = angles[2];

	data_intervel  	= 100
	# entrant functions
	ea_start 		= ea_radian
	ea_end 			= -math.pi/2	
	angle_list 		= np.linspace(ea_start, ea_end, data_intervel)
	xe = []; ye = [];
	for i in angle_list:
		xe.append( 1.5 * Rt * math.cos(i) )
		ye.append( 1.5 * Rt * math.sin(i) + 2.5 * Rt )

	#exit section
	ea_start 		= -math.pi/2
	ea_end 			= theta_n - math.pi/2
	angle_list 		= np.linspace(ea_start, ea_end, data_intervel)	
	xe2 = []; ye2 = [];
	for i in angle_list:
		xe2.append( 0.382 * Rt * math.cos(i) )
		ye2.append( 0.382 * Rt * math.sin(i) + 1.382 * Rt )

	# bell section
	# Nx, Ny-N is defined by [Eqn. 5] setting the angle to (θn – 90)
	Nx = 0.382 * Rt * math.cos(theta_n - math.pi/2)
	Ny = 0.382 * Rt * math.sin(theta_n - math.pi/2) + 1.382 * Rt 
	# Ex - [Eqn. 3], and coordinate Ey - [Eqn. 2]
	Ex = Lnp * ( (math.sqrt(aratio) - 1) * Rt )/ math.tan(math.radians(15) )
	Ey = math.sqrt(aratio) * Rt 
	# gradient m1,m2 - [Eqn. 8]
	m1 = math.tan(theta_n);  m2 = math.tan(theta_e);
	# intercept - [Eqn. 9]
	C1 = Ny - m1*Nx;  C2 = Ey - m2*Ex;
	# intersection of these two lines (at point Q)-[Eqn.10]
	Qx = (C2 - C1)/(m1 - m2)
	Qy = (m1*C2 - m2*C1)/(m1 - m2)	
	
	# Selecting equally spaced divisions between 0 and 1 produces 
	# the points described earlier in the graphical method
	# The bell is a quadratic Bézier curve, which has equations:
	# x(t) = (1 − t)^2 * Nx + 2(1 − t)t * Qx + t^2 * Ex, 0≤t≤1
	# y(t) = (1 − t)^2 * Ny + 2(1 − t)t * Qy + t^2 * Ey, 0≤t≤1 [Eqn. 6]		
	int_list = np.linspace(0, 1, data_intervel)
	xbell = []; ybell = [];
	for t in int_list:		
		xbell.append( ((1-t)**2)*Nx + 2*(1-t)*t*Qx + (t**2)*Ex )
		ybell.append( ((1-t)**2)*Ny + 2*(1-t)*t*Qy + (t**2)*Ey )
	
	# create negative values for the other half of nozzle
	nye 	= [ -y for y in ye]
	nye2  	= [ -y for y in ye2]
	nybell  = [ -y for y in ybell]
	# return
	return angles, (xe, ye, nye, xe2, ye2, nye2, xbell, ybell, nybell)

# find wall angles (theta_n, theta_e) in radians for given aratio (ar)
def find_wall_angles(ar, Rt, l_percent = 80 ):
	# wall-angle empirical data
	aratio 		= [ 4,    5,    10,   20,   30,   40,   50,   100]
	theta_n_60 	= [20.5, 20.5, 16.0, 14.5, 14.0, 13.5, 13.0, 11.2]
	theta_n_80 	= [21.5, 23.0, 26.3, 28.8, 30.0, 31.0, 31.5, 33.5]
	theta_n_90 	= [20.0, 21.0, 24.0, 27.0, 28.5, 29.5, 30.2, 32.0]
	theta_e_60 	= [26.5, 28.0, 32.0, 35.0, 36.2, 37.1, 35.0, 40.0]
	theta_e_80 	= [14.0, 13.0, 11.0,  9.0,  8.5,  8.0,  7.5,  7.0]
	theta_e_90 	= [11.5, 10.5,  8.0,  7.0,  6.5,  6.0,  6.0,  6.0]	

	# nozzle length
	f1 = ( (math.sqrt(ar) - 1) * Rt )/ math.tan(math.radians(15) )
	
	if l_percent == 60:
		theta_n = theta_n_60; theta_e = theta_e_60;
		Ln = 0.8 * f1
	elif l_percent == 80:
		theta_n = theta_n_80; theta_e = theta_e_80;
		Ln = 0.8 * f1		
	elif l_percent == 90:
		theta_n = theta_n_90; theta_e = theta_e_90;	
		Ln = 0.9 * f1	
	else:
		theta_n = theta_n_80; theta_e = theta_e_80;		
		Ln = 0.8 * f1

	# find the nearest ar index in the aratio list
	x_index, x_val = find_nearest(aratio, ar)
	# if the value at the index is close to input, return it
	if round(aratio[x_index], 1) == round(ar, 1):
		return Ln, theta_n[x_index], theta_e[x_index]

	# check where the index lies, and slice accordingly
	if (x_index>2):
		# slice couple of middle values for interpolation
		ar_slice = aratio[x_index-2:x_index+2]		
		tn_slice = theta_n[x_index-2:x_index+2]
		te_slice = theta_e[x_index-2:x_index+2]
		# find the tn_val for given ar
		tn_val = interpolate(ar_slice, tn_slice, ar)	
		te_val = interpolate(ar_slice, te_slice, ar)	
	elif( (len(aratio)-x_index) <= 1):
		# slice couple of values initial for interpolation
		ar_slice = aratio[x_index-2:len(x_index)]		
		tn_slice = theta_n[x_index-2:len(x_index)]
		te_slice = theta_e[x_index-2:len(x_index)]
		# find the tn_val for given ar
		tn_val = interpolate(ar_slice, tn_slice, ar)	
		te_val = interpolate(ar_slice, te_slice, ar)	
	else:
		# slice couple of end values for interpolation
		ar_slice = aratio[0:x_index+2]		
		tn_slice = theta_n[0:x_index+2]
		te_slice = theta_e[0:x_index+2]
		# find the tn_val for given ar
		tn_val = interpolate(ar_slice, tn_slice, ar)	
		te_val = interpolate(ar_slice, te_slice, ar)						

	return Ln, math.radians(tn_val), math.radians(te_val)

# simple linear interpolation
def interpolate(x_list, y_list, x):
	if any(y - x <= 0 for x, y in zip(x_list, x_list[1:])):
		raise ValueError("x_list must be in strictly ascending order!")
	intervals = zip(x_list, x_list[1:], y_list, y_list[1:])
	slopes = [(y2 - y1) / (x2 - x1) for x1, x2, y1, y2 in intervals]

	if x <= x_list[0]:
		return y_list[0]
	elif x >= x_list[-1]:
		return y_list[-1]
	else:
		i = bisect_left(x_list, x) - 1
		return y_list[i] + slopes[i] * (x - x_list[i])

# find the nearest index in the list for the given value
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]  

# nozzle contour plot
def plot_nozzle(ax, title, Rt, angles, contour):
	# wall angles
	nozzle_length = angles[0]; theta_n = angles[1]; theta_e = angles[2];
	
	# contour values
	xe = contour[0];   	ye = contour[1];   	nye = contour[2];
	xe2 = contour[3]; 	ye2 = contour[4];  	nye2 = contour[5];
	xbell = contour[6]; ybell = contour[7]; nybell = contour[8];
	
	# plot

	# set correct aspect ratio
	ax.set_aspect('equal')

	# throat enterant
	ax.plot(xe, ye, linewidth=2.5, color='g')
	ax.plot(xe, nye, linewidth=2.5, color='g')
	
	# throat inlet line
	x1 = xe[0]; y1 = 0;
	x2 = xe[0]; y2 = nye[0];
	dist = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
	# draw arrow, inlet radial line [x1, y1] to [x2, y2] 
	text = ' Ri = '+ str(round(dist,1))
	ax.plot(xe[0], 0, '+' )
	# draw dimension from [x1, y1] to [x2, y2] 
	ax.annotate( "", [x1, y1], [x2, y2] , arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text((x1+x2)/2, (y1+y2)/2, text, fontsize=9 )	

	# nozzle inlet length line [0,0] to [xe[0], 0]
	text = ' Li = ' + str( round( abs(xe[0]), 1) ) 
	ax.plot(0,0, '+' )
	# draw dimension from [0,0] to [xe[0], 0]
	ax.annotate( "", [0,0], [xe[0], 0], arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text( xe[0], 0, text, fontsize=9 )	
		
	# find mid point and draw arc radius
	i = int(len(xe)/2)
	xcenter = 0; 	ycenter = 2.5 * Rt;  
	xarch = xe[i];  yarch = ye[i]
	# draw arrow, enterant radial line [xcenter, ycenter] to [xarch, yarch] 
	text = ' 1.5 * Rt = '+ str( round( 1.5 * Rt, 1 ) ) 
	ax.plot(xcenter, ycenter, '+' )
	# draw dimension from [xcenter, ycenter] to [xarch, yarch]
	ax.annotate( "", [xcenter, ycenter], [xarch, yarch], arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text((xarch+xcenter)/2, (yarch+ycenter)/2, text, fontsize=9 )	
		
	# throat radius line [0,0] to [xe[-1], ye[-1]]
	text = ' Rt = '+ str(Rt)
	# draw dimension from [0,0] to [xe[-1], ye[-1]]
	ax.annotate( "", [0,0], [xe[-1], ye[-1]], arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text( xe[-1]/2, ye[-1]/2, text, fontsize=9 )	

	# throat exit
	ax.plot(xe2, ye2, linewidth=2.5, color='r')
	ax.plot(xe2, nye2, linewidth=2.5, color='r')
	# find mid point and draw arc radius
	i = int(len(xe2)/2)
	xcenter2 = 0; 	ycenter2 = 1.382 * Rt;  
	xarch2 = xe2[i];  yarch2 = ye2[i]
	# draw arrow, exit radial line from [xcenter2,ycenter2] to [xarch2, yarch2]
	text = ' 0.382 * Rt = '+ str( round(0.382 * Rt,1) ) 
	ax.plot(xcenter2, ycenter2, '+' )
	# draw dimension from [xcenter2,ycenter2] to [xarch2, yarch2]
	ax.annotate( "", [xcenter2,ycenter2], [xarch2, yarch2], arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text((xarch2+xcenter2)/2, (yarch2+ycenter2)/2, text, fontsize=9 )

	# draw theta_n, throat inflexion angle
	adj_text = 2
	origin	= [ xe2[-1], nye2[-1]-adj_text ]
	degree_symbol = r'$\theta$n'	
	draw_angle_arc(ax, theta_n, origin, degree_symbol )

	# bell section
	ax.plot(xbell, ybell, linewidth=2.5, color='b')
	ax.plot(xbell, nybell, linewidth=2.5, color='b')

	# throat radius line [0,0] to [xe[-1], ye[-1]]
	text = ' Re = ' + str( round( (math.sqrt(aratio) * Rt), 1) ) 
	ax.plot(xbell[-1],0, '+' )
	# draw dimension from [0,0] to [xe[-1], ye[-1]]
	ax.annotate( "", [xbell[-1],0], [xbell[-1], ybell[-1]], arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text( xbell[-1], ybell[-1]/2, text, fontsize=9 )	

	# draw theta_n, throat exit angle
	origin	= [ xbell[-1], nybell[-1] ]
	degree_symbol = r'$\theta$e'	
	draw_angle_arc(ax, theta_e, origin, degree_symbol )

	# nozzle length line [0,0] to [xe[-1], ye[-1]]
	text = ' Ln = ' + str( round( nozzle_length, 1) ) 
	ax.plot(0,0, '+' )
	# draw dimension from [0,0] to [xbell[-1], 0]
	ax.annotate( "", [0,0], [xbell[-1], 0], arrowprops=dict(lw=0.5, arrowstyle='<-') )
	ax.text( xbell[-1]/2, 0, text, fontsize=9 )	
				
	# axis
	ax.axhline(color='black', lw=0.5, linestyle="dashed")
	ax.axvline(color='black', lw=0.5, linestyle="dashed")		
	
	# grids
	ax.grid()
	ax.minorticks_on()
	ax.grid(which='major', linestyle='-', linewidth='0.5') # , color='red'
	ax.grid(which='minor', linestyle=':', linewidth='0.5') # , color='black'	
	
	# show
	plt.title(title, fontsize=9)
	return

# theta_n in rad,  origin =[startx, starty], degree symbol
def draw_angle_arc(ax, theta_n, origin, degree_symbol=r'$\theta$'):
	length = 50
	# start point
	startx = origin[0]; starty = origin[1];
	# find the end point
	endx = startx + np.cos(-theta_n) * length * 0.5
	endy = starty + np.sin(-theta_n) * length * 0.5
	# draw the angled line
	ax.plot([startx,endx], [starty,endy], linewidth=0.5, color='k')
	# horizontal line
	# ax.hlines(y=starty, xmin=startx, xmax=length, linewidth=0.5, color='k')
	# angle
	arc_obj = Arc([startx, starty], 1, 1, 0, 0, math.degrees(theta_n), color='k' )
	ax.add_patch(arc_obj)
	ax.text(startx+0.5, starty+0.5, degree_symbol + ' = ' + str(round(theta_n,1)) + u"\u00b0")	
	return

# ring of radius r, height h, base point a
def ring(r, h, a=0, n_theta=30, n_height=10):
    theta = np.linspace(0, 2*np.pi, n_theta)
    v = np.linspace(a, a+h, n_height )
    theta, v = np.meshgrid(theta, v)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    z = v
    return x, y, z

# Set 3D plot axes to equal scale. 
# Required since `ax.axis('equal')` and `ax.set_aspect('equal')` don't work on 3D.
def set_axes_equal_3d(ax: plt.Axes):
    """	
    https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
    """
    limits = np.array([
        ax.get_xlim3d(),
        ax.get_ylim3d(),
        ax.get_zlim3d(),
    ])
    origin = np.mean(limits, axis=1)
    radius = 0.5 * np.max(np.abs(limits[:, 1] - limits[:, 0]))
    _set_axes_radius(ax, origin, radius)
    return

# set axis limits
def _set_axes_radius(ax, origin, radius):
    x, y, z = origin
    ax.set_xlim3d([x - radius, x + radius])
    ax.set_ylim3d([y - radius, y + radius])
    ax.set_zlim3d([z - radius, z + radius])
    return

# 3d plot
def plot3D(ax, contour):
	# unpack the contour values
	xe = contour[0];   	ye = contour[1];   	nye = contour[2];
	xe2 = contour[3]; 	ye2 = contour[4];  	nye2 = contour[5];
	xbell = contour[6]; ybell = contour[7]; nybell = contour[8];
	# collect and append array values
	x = []; y = [];
	x = np.append(x, xe);  y = np.append(y, ye)	
	x = np.append(x, xe2);  y = np.append(y, ye2)	
	x = np.append(x, xbell);  y = np.append(y, ybell)	
	
	# ring thickness
	thick = 5 * (x[1] - x[0]) # 0.01
	
	# draw multiple rings to create 3d structure
	for i in range(len(y)):
		X, Y, Z = ring(y[i], thick, x[i])
		ax.plot_surface(X, Y, Z, color='g')	
		
	# set correct aspect ratio
	ax.set_box_aspect([1,1,1])
	set_axes_equal_3d(ax)
	# set view
	ax.view_init(-170, -15)
	return

def plot(title, throat_radius, angles, contour):
	# Plot 3d view
	fig = plt.figure(figsize=(12,9))
	# plot some 2d information
	ax1 = fig.add_subplot(121)
	plot_nozzle(ax1, title, throat_radius, angles, contour)
	# plot 3d view
	ax2 = fig.add_subplot(122, projection='3d')
	plot3D(ax2, contour)	
	# show
	fig.tight_layout(rect=[0, 0.03, 1, 0.95])
	plt.show()
	return

# __main method__
if __name__=="__main__":
	
	# constants
	k 	= 1.21		# ratio of specific heats
	l_percent = 80	# nozzle length percntage (60, 80, 90)
	
	# typical upper stage values
	aratio = 25 			# Ae / At	
	throat_radius = 40 		# {'radius_throat': 40, 'radius_exit': 210}		

	# rao_bell_nozzle_contour
	angles, contour = bell_nozzle(k, aratio, throat_radius, l_percent)
	# plot contour
	title = 'Bell Nozzle \n [Area Ratio = ' + str(round(aratio,1)) + ', Throat Radius = ' + str(round(throat_radius,1)) + ']' 
	plot(title, throat_radius, angles, contour)

	# --------------- Nozzle-2------------------

	# typical lower stage booster values
	aratio = 7 				# Ae / At	
	throat_radius = 800 	# 	

	# rao_bell_nozzle_contour
	angles, contour = bell_nozzle(k, aratio, throat_radius, l_percent)
	# plot contour
	title = 'Bell Nozzle \n [Area Ratio = ' + str(round(aratio,1)) + ', Throat Radius = ' + str(round(throat_radius,1)) + ']' 
	plot(title, throat_radius, angles, contour)
