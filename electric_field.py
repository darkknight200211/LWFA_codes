import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import scipy.constants as cn
from scipy.interpolate import make_interp_spline
import cmath
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FormatStrFormatter

# all in cgs units
# -----enter number density n------------------------------------------------
num_den = 1e16      # units is in cm^-3
# calculating plasma frequency
me = cn.electron_mass*(1e3)
print("me: in g ", me)
e = cn.elementary_charge*(3.0e9)
print("e: in statcoloumb ", e)
c = cn.speed_of_light * 1e2
print("c: ", c)
# wpe = np.sqrt((4.0*cn.pi*num_den)/me) * e
wpe = 5.64e4*np.sqrt(num_den)
print("wpe: ", wpe, "   ", "{:.2e}".format(wpe))
lamblas = 1e-4                  # laser wavelength in cm
print("laser wavelength: ", lamblas, "  ", "{:.2e}".format(lamblas))
k0 = 2.0*cn.pi/lamblas
# print(wpe)
w0 = (2.0*cn.pi*c)/lamblas     # this is laser frequency
print("w0: ", w0, "  ", "{:.2e}".format(w0))
t0 = (2.0*cn.pi)            # dimless
tl = 600.0/(cn.pi)        # dimless this is the finite duration of laser pulse
print("tl, T0: s ", tl/w0, t0/w0, tl/t0)
eta = np.sqrt(1.0 - (wpe/w0)**2)    # this is refrective index of plasma
# eta = 1.0

print("eta: ", eta)
vg = c/eta
print("vg: ", vg)        # speed of light in plasma
lampe = (2.0*cn.pi*vg)/wpe      # dim plasma wavelength
print("plasma wavelength:", lampe, "  ", "{:.2e}".format(lampe))
print("plasma wavelength/laser wavelength: ", lampe/lamblas)
nwaves = 3.0        # this will decide the system length
                    # how many plasma wavelength long is it
l = nwaves*lampe    # dim system lenght
print("l: ", l)
area = 1.0          # system area
nr = l*area*num_den # this is number of real particles
print("nr: ", nr)
#nc = 1000001          # this is number of computational particles
nc = 1000001
wgt = nr/(nc-1)         # this is the weight of each comp. particle
print("weight: ", wgt)
ld = ((w0/wpe)**2)*lampe
print("depletion length: ", ld, "   ", "{:.2e}".format(ld))
wgt_arr = np.zeros(nc)
wgt_arr[:] = wgt
phdel = l/(nc-1)        # this is the the distance between consecutive
                        # com. particles
com_space = np.linspace(0,l,nc)     # creating a list with comp. particles
                                    # evenly spread
# tl = 60.0/cn.pi               # this is the pulse duration
# delt = 1.0/np.sqrt(2)         
delt = 1.0
#                               # this is the polarization for the laser
                                # 1.0 means plane polarized, 1.0/sqrt(2.0)
                                # circularly polarized

def wba(t, delt, tl, lampe, lamblas):
    denom = 0.0
    arg1 = tl/(2.0*cn.pi)
    arg2 = tl/(4.0*cn.pi)
    arg3 = tl/(2.0*cn.pi + 2.0*tl)
    arg4 = tl/(2.0*cn.pi - 2.0*tl)
    arg5 = tl/(4.0*cn.pi + 2.0*tl)
    arg6 = tl/(4.0*cn.pi - 2.0*tl)
    p1 = np.sin((2.0*cn.pi*t)/tl)
    p2 = np.sin((4.0*cn.pi*t)/tl)
    p3 = np.sin((2.0*cn.pi*t + 2.0*tl*t)/tl)
    p4 = np.sin((2.0*cn.pi*t - 2.0*tl*t)/tl)
    p5 = np.sin((4.0*cn.pi*t + 2.0*tl*t)/tl)
    p6 = np.sin((4.0*cn.pi*t - 2.0*tl*t)/tl)

    denom = (((3.0/8.0)*t - arg1*p1/2.0 + arg2*p2/8.0) + (2.0*(delt**2)-1.0) *(3.0*np.sin(2.0*t)/(16.0) -(1.0/4.0)*arg3*p3 - (1.0/4.0)*arg4*p4 +(1.0/16.0)*arg5*p5 +(1.0/16.0)*arg6*p6 ))
    val = 2.0*np.sqrt((lampe)/(denom*lamblas))
    return val

# con = (3.0/16.0) + (tl**2/(4.0 * (cn.pi)**2 - 4.0 * (tl**2))) - ((2.0 * (tl**2))/(128.0 * (cn.pi**2)- 32.0 * (tl**2)))
# denom = ((3.0* (tl**2))/8.0) + (con * np.sin(2.0*tl))

v = 1.0    # this is the multiplier for deciding mag of a0.

if (delt == 1.0):
    # a0 = v * 2.0 * np.sqrt(lampe/(lamblas*denom))
    a0 = v*wba(tl, delt, tl, lampe, lamblas)
elif (delt == 1.0/np.sqrt(2)):
    a0 = v * np.sqrt((32.0*lampe)/(3.0*lamblas*tl))

#----defining how oscillator displaces when laser is passing-----------------
def xdis(t, a0, delt, tl):
    p = 0.0
    arg1 = tl/(2.0*cn.pi)
    arg2 = tl/(4.0*cn.pi)
    arg3 = tl/(2.0*cn.pi + 2.0*tl)
    arg4 = tl/(2.0*cn.pi - 2.0*tl)
    arg5 = tl/(4.0*cn.pi + 2.0*tl)
    arg6 = tl/(4.0*cn.pi - 2.0*tl)
    p1 = np.sin((2.0*cn.pi*t)/tl)
    p2 = np.sin((4.0*cn.pi*t)/tl)
    p3 = np.sin((2.0*cn.pi*t + 2.0*tl*t)/tl)
    p4 = np.sin((2.0*cn.pi*t - 2.0*tl*t)/tl)
    p5 = np.sin((4.0*cn.pi*t + 2.0*tl*t)/tl)
    p6 = np.sin((4.0*cn.pi*t - 2.0*tl*t)/tl)

    p = (a0**2 / 4.0)*(((3.0/8.0)*t - arg1*p1/2.0 + arg2*p2/8.0) + (2.0*(delt**2)-1.0) *(3.0*np.sin(2.0*t)/(16.0) -(1.0/4.0)*arg3*p3 - (1.0/4.0)*arg4*p4 +(1.0/16.0)*arg5*p5 +(1.0/16.0)*arg6*p6 ))
    val = p*(vg/w0)
    # val = p
    return val


# xp = w0*me*(c**3)*xdis(tl, a0, delt, tl)
xp = xdis(tl, a0, delt, tl)
print("xmax: ", xp, 2.0*cn.pi*xp, "k_{0}*xmax:", (2.0*cn.pi*xp)/lamblas, xp/lampe)
print(r"$\delta= $",phdel)
print("xmax/lampe: ", 2.0*cn.pi*xp/lampe)
print("xmax/phdel: ", xp/phdel)

# -------pendulum model------------------------------------------------------
k = 0               # will be used to keep track of how many times did the
                    # i-th oscillator hit xp in the time that the laser took
                    # to cross the distance from i-th to 0th oscillator
tn = 0.0
ossx = np.zeros(nc) # comp. oscillator position
phii = 0.0          # this is the phase of the ith oscillator
labx = np.zeros(nc) # location in lab after laser has passed
lab = np.zeros(nc)  # location in lab before laser has passed
for i in range(nc):
    lab[i] = -(nc-i)*phdel
diff = np.zeros(nc) # the difference between oscillation position of neighbor

tsim = 0.0*tl  # time at which snapshot of the oscillators positions is
                # needed. Note: this time is taken to be t-nos*phdel
                # it is the time passed after the laser has transversed the
                # full system.

# --------binning the comp. oscillators into intevals of size xbin-----------
xbin = lampe/50.0       # size of one bin
nbin = int(l/xbin)           # no. of bins in the system length
print("xbin, nbin: ", xbin, nbin)
ndenc = np.zeros(nbin)
ndenr = np.zeros(nbin)
binspace = np.linspace(-l, 0, nbin)

bnossx = np.zeros(nbin)
bnlabx = np.zeros(nbin)
bnlab = np.zeros(nbin)
for i in range(nbin):
    bnlab[i] = -(nbin-i)*xbin

# -----density before laser hits the plasma----------------------------------
hist, bin_edges = np.histogram(lab, nbin, (-l,0))
# print(hist)
# print(bin_edges)
print(len(hist), len(bin_edges))

mean_edges = np.zeros(nbin)
for i in range(len(bin_edges)-1):
    mean_edges[i] = (bin_edges[i]+bin_edges[i+1])/2.0

for i in range(nbin):
    ndenc[i] = hist[i]/xbin
    ndenr[i] = ndenc[i]*wgt     # multipied by weight to get real density

f_ne = "{:.2e}".format(num_den)
# plt.title(r"$n_{e}=$" + f_ne)
# plt.plot(mean_edges, ndenr)
# plt.ylim(0.0, 1.5*num_den)
# plt.xlabel(r"$l_{sys}$ (cm)")
# plt.ylabel(r"$n_{e} \ (cm^{-3})$")
# plt.show()

# -------pendulum model------------------------------------------------------
for i in range(nc):
    tn = -(lab[i]*w0)/(vg)                # dimensional
    k = int((tn*(wpe/w0))/(2.0*cn.pi))   # dimensinless k in terms of dimensional
    phii = 2.0*cn.pi*k - tn*(wpe/w0)
    ossx[i] = xp*np.cos(tsim*(wpe/w0) + (phii))  # dimenional, cos arguments dimless
    labx[i] = lab[i] + ossx[i]

# the following loop is run to compute binned quantities. To find longitudinal
# electric field we do: delta(n_{e})*(oscillator position from equilibrium),
# however, after running above loop there are as many as nc elements in
# ossx. But we need only nbin number positions. So the following loop was made.
# The quantities bnlab means binned-lab position, bnlabx is binned-lab x
# bnossx is binned ossx.

for i in range(nbin):
    tn = -(bnlab[i]*w0)/(vg)                # dimensional
    k = int((tn*(wpe/w0))/(2.0*cn.pi))   # dimensinless k in terms of dimensional
    phii = 2.0*cn.pi*k - tn*(wpe/w0)
    bnossx[i] = xp*np.cos(tsim*(wpe/w0) + (phii))  # dimenional, cos arguments dimless
    bnlabx[i] = bnlab[i] + bnossx[i]

#--------density profile tsim time after the laser has passed----------------
hist, bin_edges = np.histogram(labx, nbin, (-l,0))
# print(wgt*hist)
# print(bin_edges)

mean_edges = np.zeros(nbin)
for i in range(len(bin_edges)-1):
    mean_edges[i] = (bin_edges[i]+bin_edges[i+1])/2.0

diffn = 0.0
diffr = 0.0
ex = np.zeros(nbin)
den_diff = np.zeros(nbin)

for i in range(nbin):
    ndenc[i] = hist[i]/xbin
    ndenr[i] = wgt*ndenc[i]
    diffn = num_den-ndenc[i]
    diffr = num_den-ndenr[i]
    den_diff[i] = diffr
    ex[i] = 4.0*cn.pi*e*diffn*bnossx[i]

max_num_den = max(ndenr)
print("max num den: ", max_num_den)
print( "maximum num den diff: ", max(den_diff))
max_diff = "{:.2e}".format(max(den_diff))
# max_diff = max(den_diff)

print("nc: ", int(lampe*nwaves/xbin))

# fig, axs = plt.subplots(2, sharex=True)

# print(lab[0], mean_edges[0], labx[0])
# axs[0].plot(bnlabx, bnossx, label = r"$\delta x_{i} vs. x_{i}+\delta x_{i}$")
# axs[0].legend()
# axs[0].set_ylabel(r"$\delta x_{i} \ (cm)$", fontsize=15)
# # axs[2].plot(bnlabx, den_diff)
# axs[1].plot(mean_edges, ndenr, label="real density")
# # axs[1].plot(bnlabx, ndenr, label="real density")
# axs[1].legend()
# axs[1].set_ylim(0.0, 1.2*max_num_den)
# # axs[1].axhline(y=num_den)
# axs[1].set_ylabel(r"$n_{e} \ (cm^{-3})$", fontsize=15)
# axs[1].set_xlabel(r"$l_{sys}$ (cm)", fontsize = 15)

# axex = axs[1].twinx()

# axex.plot(bnlabx, ex, label=r"$E_{x}$", color = "red")
# # axex.axhline(y=0.0, color='red')
# axex.set_ylabel(r"$E_{x}\ statvolt\ cm^{-1}$", fontsize=15)
# axex.legend()

# # ax4 = axs[1].twinx()
# # ax4.plot(bnlabx, den_diff, color='red')

# # ax5 = axs[0].twinx()
# # ax5.plot(bnlabx, bnossx, color='red')

# # fig.suptitle(r"$t_{sim}= $"+ str(round(tsim,3))+ "      " +r"$a_{0}$ = " + str(round(a0,3)) + "     " + r"$x_{max}=$" + str(round(xp,3))+ "       " + r"$\lambda_{p}= $" + str(round(lampe/(2.0*cn.pi),3)) + "        " + r"$\frac{x_{max}}{\lambda_{p}}=$"+str(round(2.0*cn.pi*xp,3)/lampe)+"        " + str((lampe)/lamblas))
# fig.suptitle(r"$a_{0}$="+str(round(a0,4)) + "        " + r"$x_{max}=$" + str(round(xp,4))  +  "     " + r"$2.0\pi \ x_{max}=$" + str(round(2.0*cn.pi*xp,4)) + "       " + r"$\lambda_{p}=$" + str(round(lampe,4)) + "     " + r"$\tilde{x}_{max}$=" + str(round(2.0*cn.pi*xp/lamblas,4)) + "      " + r"$2\pi \ \frac{\lambda_{p}}{\lambda_{0}}=$" + str(round( lampe/lamblas ,4)) + "      " + r"$\Delta n_{e,max}=$" + str(max_diff))
# plt.show()

exmax = "{:.3e}".format(max(ex))
print("maximum longitudinal electric field: ", exmax)

fig, axs = plt.subplots()

axs.plot(mean_edges[8:]/1e-1, ndenr[8:]/1e17, label=r"$n_{e}$")
# axs.axhline(y=num_den)
axs.legend()
axs.set_ylabel(r"$\mathbf{n_{e}\times 10^{17}\ cm^{-3}}$", fontsize=13)
for ax in np.ravel(axs):

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='sci',
                        axis='both',
                        scilimits=(0, 0))  # Force scientific notation

    ax.tick_params(
        axis = 'both',
        direction = 'in',
        width = 2,
        color='black',  # Tick color
        labelcolor='k',  # Label color
        labelsize='12', #label size
        top = True,
        right = True
    )
    # Set fontweight of tick labels directly
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight('bold')
    # for label in ax.get_xlabel() + ax.get_ylabel():
    #     label.set_fontsize(15)
     # Adjust legend
    ax.legend(loc='center left',
              bbox_to_anchor=(0, 0.025),
              frameon=False,
              fontsize=13
              )
axs.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))  # For left y-axis

ax1 = axs.twinx()
ax1.plot(bnlabx/1e-1, ex/1e5, label=r"$E_{x}$", color='red')
ax1.tick_params(
        axis = 'both',
        direction = 'in',
        width = 2,
        color='black',  # Tick color
        labelcolor='k',  # Label color
        labelsize='12', #label size
        top = True,
        right = True
    )
for label in axs.get_xticklabels()+axs.get_yticklabels():
        label.set_fontweight('bold')
# ax1.axhline(y=0.0, color='red')
ax1.legend(frameon=False)
ax1.set_ylabel(r"$\mathbf{E_{x}\times 10^{5}}$", fontsize = 13)
axs.set_xlabel(r"$\mathbf{l_{sys}\times 10^{-1}}$ (cm)", fontsize=13, fontweight='bold')
ax1.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))
fig.suptitle(r"$\mathbf{E_{x,max}=}$" + str(exmax), fontweight='bold', fontsize=15)


plt.show()

#****************************************************************************

# for i in range(nbin):
#     xr = -i*xbin
#     xl = -(i+1)*xbin
#     counter = 0
#     for j in range(nc):

#         if (labx[j] <= xr and labx[j] > xl):
#             counter = counter + 1
#         else:
#             counter = counter
    
#     ndenc[i] = counter/xbin
