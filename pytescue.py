import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons

###########################################################
################# Define Constants ########################
###########################################################

pi = np.pi

# colors
col_abc = ['b','g','r']                           # abc phasors
col_pn0 = ['tab:brown','tab:olive','tab:cyan']    # positive, negative and zero phasors


lim_x = (-1.25, 1.25)
lim_y = (-1.25, 1.25)
# define operator a: complex 1<120º
a = 1*np.exp(1j*2*pi/3)


#################################
####### Calculate Vabc and Vpn0
#################################

def calculate(mag, theta):

    # open data for each phase
    mag_a = mag[0]
    mag_b = mag[1]
    mag_c = mag[2]

    theta_a = theta[0]
    theta_b = theta[1]
    theta_c = theta[2]

    # calculate new phasors
    Va = mag_a*np.exp(+1j*( 0      + theta_a))
    Vb = mag_b*np.exp(-1j*(2*pi/3  + theta_b))
    Vc = mag_c*np.exp(+1j*(2*pi/3  + theta_c))

    Vabc = [Va, Vb, Vc]      # agregate data

    # Define Fortescue Transform
    Vap = K*(Va +   Vb*a    +  Vc*a**2)
    Vam = K*(Va +   Vb*a**2 +  Vc*a)
    Va0 = K*(Va +   Vb      +  Vc)

    Vpn0_a = [Vap, Vam, Va0]
    Vpn0_b = [Vap*a**2, Vam*a, Va0]
    Vpn0_c = [Vap*a, Vam*a**2, Va0]

    Vfortescue = [Vpn0_a,Vpn0_b,Vpn0_c]      # agregate data

    return Vabc, Vfortescue



def plots(Vabc, Vfortescue, both=True):

    if both:
        for pn0 in Vfortescue:
            aux = 0
            off_x = 0
            off_y = 0

            for phase_pn0 in pn0:
                color_p = col_pn0[aux]
                pn0 = phase_pn0
                if abs(pn0) > 0.01:
                    ax.arrow(off_x, off_y, np.real(pn0), np.imag(pn0), color=color_p, length_includes_head=True,
                             width=0.005, head_width=None, head_length=None,overhang=0.6) 
                    off_x += np.real(pn0)
                    off_y += np.imag(pn0)
                aux += 1

    aux = 0 
    for abc in Vabc:
        # plot abc voltage phasors
        color_a = col_abc[aux]
        if abs(abc) > 0.01:
            ax.arrow(0, 0, np.real(abc), np.imag(abc), color=color_a, length_includes_head=True, width=0.005, head_width=None, head_length=None,overhang=0.6) 
        
        aux += 1

    # add text legend for ABC components
    ax.text(1, 1, "Phase A \n \n \n", size=10,
            ha="center", va="center", c=col_abc[0], bbox=dict(facecolor='white', ec='grey'))

    ax.text(0.99, 1, "Phase B", size=10,
            ha="center", va="center", c=col_abc[1])

    ax.text(0.99, 0.89, "Phase C", size=10,
            ha="center", va="center", c=col_abc[2])

    # add text legend for PN0
    ax.text(1, -1, " Positive \n \n \n", size=10,
            ha="center", va="center", c=col_pn0[0], weight='bold', bbox=dict(facecolor='white', ec='grey'))

    ax.text(0.99, -1, " Negative", size=10,
            ha="center", va="center", c=col_pn0[1], weight='bold')

    ax.text(0.99, -1.1, "Zero", size=10,
            ha="center", va="center", c=col_pn0[2], weight='bold')



############################################################
################ Initial values plot #######################
############################################################

# define module and phase for ABC quantities

mag_a = 1
mag_b = 1
mag_c = 1
mag = [mag_a, mag_b, mag_c]

theta_a = 0*(20)*pi/180
theta_b = 0*(60)*pi/180
theta_c = 0*(120)*pi/180
theta = [theta_a, theta_b, theta_c]

# Fortescue constant (use sqrt(3)/3 for power-invariant transform)
K = 1/3  

# Call function
Vabc, Vfortescue = calculate(mag, theta)

# Plot stuff
fig, ax = plt.subplots()
ax.set(xlim=lim_x, ylim=lim_y)
ax.set_box_aspect(1)
ax.grid(ls=':')


plots(Vabc, Vfortescue, both=False)

###################################################################
######################### Define Sliders ##########################
###################################################################

axcolor = 'lightgrey'
        
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
axmag_a   = plt.axes([0.02, 0.1, width , 0.78], facecolor=axcolor)
axtheta_a = plt.axes([0.05, 0.1, width, 0.78], facecolor=axcolor)

# phase b
axmag_b   = plt.axes([0.10, 0.1, width , 0.78], facecolor=axcolor)
axtheta_b = plt.axes([0.13, 0.1, width, 0.78], facecolor=axcolor)

# phase c
axmag_c   = plt.axes([0.18, 0.1, width , 0.78], facecolor=axcolor)
axtheta_c = plt.axes([0.21, 0.1, width, 0.78], facecolor=axcolor)

# reset button
resetax = plt.axes([0.0854, 0.003, 0.08, 0.04])

###### create sliders

step_mag = 0.01

# phase a
slid_mag_a = Slider(ax=axmag_a,label=r"$|f_{A}|$",valmin=mag_min,valmax=mag_max,valinit=mag_ini,valstep=step_mag,
orientation="vertical",color=col_abc[0])
slid_theta_a = Slider(ax=axtheta_a,label=r"$\theta_{A}$",valmin=theta_min,valmax=theta_max,valinit=theta_ini,valstep=0.1,
orientation="vertical",color=col_abc[0])

# phase b
slid_mag_b = Slider(ax=axmag_b,label=r"$|f_{B}|$",valmin=mag_min,valmax=mag_max,valinit=mag_ini,valstep=step_mag,
orientation="vertical",color=col_abc[1])
slid_theta_b = Slider(ax=axtheta_b,label=r"$\theta_{B}$",valmin=theta_min,valmax=theta_max,valinit=theta_ini,valstep=0.1,
orientation="vertical",color=col_abc[1])

# phase c
slid_mag_c = Slider(ax=axmag_c,label=r"$|f_{C}|$",valmin=mag_min,valmax=mag_max,valinit=mag_ini,valstep=step_mag,
orientation="vertical",color=col_abc[2])
slid_theta_c = Slider(ax=axtheta_c,label=r"$\theta_{C}$",valmin=theta_min,valmax=theta_max,valinit=theta_ini,valstep=0.1,
orientation="vertical",color=col_abc[2])

# reset button
button = Button(resetax, 'Reset', color=axcolor, hovercolor='0.975')

############################################################################
############################## Def Functions ###############################
############################################################################ 

def update(val):

    # update magnitudes and phase for ABC quantities
    theta_a = slid_theta_a.val*pi/180
    mag_a = slid_mag_a.val

    theta_b = slid_theta_b.val*pi/180
    mag_b = slid_mag_b.val

    theta_c = slid_theta_c.val*pi/180
    mag_c = slid_mag_c.val

    # update vector with new data
    theta = [theta_a, theta_b, theta_c]
    mag = [mag_a, mag_b, mag_c]
    
    # call function
    Vabc, Vfortescue = calculate(mag, theta)

    ax.cla()
    ax.grid(ls=':')
    ax.set(xlim=lim_x, ylim=lim_y)

    plots(Vabc, Vfortescue)
        
    fig.canvas.draw_idle()


# update graph when user change a parameter
slid_mag_a.on_changed(update)
slid_theta_a.on_changed(update)

slid_mag_b.on_changed(update)
slid_theta_b.on_changed(update)

slid_mag_c.on_changed(update)
slid_theta_c.on_changed(update)

# define reset button

def reset(event):
    slid_mag_a.reset()
    slid_theta_a.reset()

    slid_mag_b.reset()
    slid_theta_b.reset()

    slid_mag_c.reset()
    slid_theta_c.reset()
    
button.on_clicked(reset)



plt.show()

a=1
