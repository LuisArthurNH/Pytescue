import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

###########################################################
################# signal definition #######################
###########################################################
pi = np.pi

lim_x = (-1.1, 1.1)
lim_y = (-0.9, 0.9)

# define operator a: complex 1<120º
a = 1*np.exp(1j*2*pi/3)

# define module and phase for ABC quantities

theta_a = 0*(20)*pi/180
theta_b = 0*(60)*pi/180
theta_c = 0*(+120)*pi/180

Va = 1*np.exp(+1j*( 0      + theta_a))
Vb = 1*np.exp(-1j*(2*pi/3  + theta_b))
Vc = 1*np.exp(+1j*(2*pi/3  + theta_c))

Vabc = [Va, Vb, Vc]
col_abc = ['b','g','r']

# Define Fortescue Transform

K = 1/3  # Fortescue constant (use sqrt(3)/3 for power-invariant transform)

Vap = K*(Va +   Vb*a    +  Vc*a**2)
Vam = K*(Va +   Vb*a**2 +  Vc*a)
Va0 = K*(Va +   Vb      +  Vc)

Vpm0_a = [Vap, Vam, Va0]
Vpm0_b = [Vap*a**2, Vam*a, Va0]
Vpm0_c = [Vap*a, Vam*a**2, Va0]

Vfortescue = [Vpm0_a,Vpm0_b,Vpm0_c]

col_pm0 = ['k','aqua','blueviolet']

# Plot stuff
fig, ax = plt.subplots()
# ax.set_xlim([-1.2,1.2])
# ax.set_ylim([-1,1])
ax.set(xlim=lim_x, ylim=lim_y)
ax.set_box_aspect(1)
ax.grid(ls=':')

aux = 0 
for abc in Vabc:

    # plot abc voltage phasors
    color_a = col_abc[aux]
    if abs(abc) > 0.01:
        ax.arrow(0, 0, np.real(abc), np.imag(abc), color=color_a, length_includes_head=True, width=0.005, head_width=None, head_length=None,overhang=0.6) 
    
    aux += 1

# # plot +, -, 0 phasors
# aux = 0
# off_x = 0
# off_y = 0

# # plt.show()

# for pm0 in Vfor:
#     color_p = col_pm0[aux]
#     pm0 = Vfor[aux]
#     if abs(abc) > 0.01:
#         ax.arrow(off_x, off_y, np.real(pm0), np.imag(pm0), color=color_p, length_includes_head=True, width=0.005, head_width=None, head_length=None,overhang=0.6) 
#         off_x += np.real(pm0)
#         off_y += np.imag(pm0)
#     aux += 1


###### Slider data
# Magnitude limits
mag_ini = 1                        # initial frequency plot
mag_max = 1.1                      # max frequency on range [Hz]
mag_min = 0                        # min freq


###### theta limits
theta_ini = 0                        # initial frequency plot
theta_max = 180                      # max frequency on range [Hz]
theta_min = -180                        # min freq



axcolor = 'lightgoldenrodyellow'

width = 0.02

# phase a
axmag_a   = plt.axes([0.07, 0.1, width , 0.78], facecolor=axcolor)
axtheta_a = plt.axes([0.1, 0.1, width, 0.78], facecolor=axcolor)

# phase b
axmag_b   = plt.axes([0.15, 0.1, width , 0.78], facecolor=axcolor)
axtheta_b = plt.axes([0.18, 0.1, width, 0.78], facecolor=axcolor)

# phase c
axmag_c   = plt.axes([0.23, 0.1, width , 0.78], facecolor=axcolor)
axtheta_c = plt.axes([0.26, 0.1, width, 0.78], facecolor=axcolor)

###### create sliders

slid_mag_a = Slider(
    ax=axmag_a,
    label="Ma",
    valmin=mag_min,
    valmax=mag_max,
    valinit=mag_ini,
    valstep=0.1,
    orientation="vertical",
    color='b'
)

slid_theta_a = Slider(
    ax=axtheta_a,
    label="Aa",
    valmin=theta_min,
    valmax=theta_max,
    valinit=theta_ini,
    valstep=0.1,
    orientation="vertical",
    color='b'
)


slid_mag_b = Slider(
    ax=axmag_b,
    label="Mb",
    valmin=mag_min,
    valmax=mag_max,
    valinit=mag_ini,
    valstep=0.1,
    orientation="vertical",
    color='g'
)

slid_theta_b = Slider(
    ax=axtheta_b,
    label="Ab",
    valmin=theta_min,
    valmax=theta_max,
    valinit=theta_ini,
    valstep=0.1,
    orientation="vertical",
    color='g'
)



slid_mag_c = Slider(
    ax=axmag_c,
    label="Mc",
    valmin=mag_min,
    valmax=mag_max,
    valinit=mag_ini,
    valstep=0.1,
    orientation="vertical",
    color='r'
)

slid_theta_c = Slider(
    ax=axtheta_c,
    label="Ac",
    valmin=theta_min,
    valmax=theta_max,
    valinit=theta_ini,
    valstep=0.1,
    orientation="vertical",
    color='r'
)



############################################################################
############################### Def Function ###############################
############################################################################ 

def update(val):

    theta_a = slid_theta_a.val*pi/180
    mag_a = slid_mag_a.val

    theta_b = slid_theta_b.val*pi/180
    mag_b = slid_mag_b.val

    theta_c = slid_theta_c.val*pi/180
    mag_c = slid_mag_c.val
    
    Va = mag_a*np.exp(+1j*( 0      + theta_a))
    Vb = mag_b*np.exp(-1j*(2*pi/3  + theta_b))
    Vc = mag_c*np.exp(+1j*(2*pi/3  + theta_c))

    Vabc = [Va, Vb, Vc]

    Vap = K*(Va +   Vb*a    +  Vc*a**2)
    Vam = K*(Va +   Vb*a**2 +  Vc*a)
    Va0 = K*(Va +   Vb      +  Vc)

    Vpm0_a = [Vap, Vam, Va0]
    Vpm0_b = [Vap*a**2, Vam*a, Va0]
    Vpm0_c = [Vap*a, Vam*a**2, Va0]

    Vfortescue = [Vpm0_a,Vpm0_b,Vpm0_c]


    ax.cla()
    ax.grid(ls=':')
    ax.set(xlim=lim_x, ylim=lim_y)

    aux = 0 
    for abc in Vabc:
        # plot abc voltage phasors
        color_a = col_abc[aux]
        if abs(abc) > 0.01:
            ax.arrow(0, 0, np.real(abc), np.imag(abc), color=color_a, length_includes_head=True, width=0.005, head_width=None, head_length=None,overhang=0.6) 
        
        aux += 1



    for pm0 in Vfortescue:
        aux = 0
        off_x = 0
        off_y = 0

        for phase_pm0 in pm0:
            color_p = col_pm0[aux]
            pm0 = phase_pm0
            if abs(abc) > 0.01:
                ax.arrow(off_x, off_y, np.real(pm0), np.imag(pm0), color=color_p, length_includes_head=True, width=0.005, head_width=None, head_length=None,overhang=0.6) 
                off_x += np.real(pm0)
                off_y += np.imag(pm0)
            aux += 1


        
    fig.canvas.draw_idle()


# update graph when user change a parameter
slid_mag_a.on_changed(update)
slid_theta_a.on_changed(update)

slid_mag_b.on_changed(update)
slid_theta_b.on_changed(update)

slid_mag_c.on_changed(update)
slid_theta_c.on_changed(update)




plt.show()

a=1
