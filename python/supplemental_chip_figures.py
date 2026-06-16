# Imports
import os

OUTPUT_DIR = os.environ.get("OUTPUT_DIR", "/app/output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

import numpy as np
from scipy.integrate import odeint
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import gridspec
import matplotlib.ticker as ticker

# Functions
def Hematopoiesis_2CellJohn(y, t, k, l):
    ''' TKI therapy
        Two parallel lineages:
        1. Normal Development
        2. Leukemic Development
        Differences between lineages:
        1. Feedback strengths are different,
            a.gami for normal development
            b.gami_L for leukemic development
        2. Decay term associated to cells
            a. d_HSC
            b. d_T
    '''

    (p0max,eta0max,
     gam0,gam1,
     phi0,c1,c2,c3,
     dHSC,dT,
     n1,n2) = k
    (p0max_L,eta0max_L,
     gam0_L,gam1_L,phi0_L,
     alpha,alpha2,alpha3,dT_L,
     n1L,n2L,n3L) = l

    dydt = np.zeros(4)

    y = np.where(y>1e-8,y,0)

    #FBKS Normal
    p0 = p0max/(1+(gam0*(c1*y[1]+c2*y[3])**n1)/(1+(phi0*phi0_L*(y[2]+c3*y[0]))**n2))
    eta0 = eta0max/(1+gam1*(y[0]+y[2]))

    #FBKS Leukemic
    p0_L = p0max_L/(1+(gam0_L*(c1*y[1]+c2*y[3])**n1L)/(1+(phi0_L*(y[2]+c3*y[0]))**n2L))
    eta0_L = eta0max_L/(1+gam1_L*(y[0]+y[2]))
    dHSC_L = alpha*(alpha2+(alpha3*y[2])**n3L/(1+(alpha3*y[2])**n3L))

    #Normal Lineage
    dydt[0] = (2*p0-1-dHSC)*eta0*y[0]
    dydt[1] = 2*(1-p0)*eta0*y[0] - dT*y[1]
    #Luekemic Lineage
    dydt[2] = (2*p0_L-1-dHSC_L)*eta0_L*y[2]
    dydt[3] = 2*(1-p0_L)*eta0_L*y[2] - dT_L*y[3]

    return dydt
# CHIP Mutant Experiments
## Control (no CHIP mutation)
# Initialize
k1 = [0.719673,1,2.4408113329290927*1e-6,0,0,1,1,.03,0,.1,1,1]
l1 = [0.719673,1,(2.4408113329290927*1e-6)*(0.8511500334069608),0,0.0016270262087911059,0.4093047,0.04515,0.01633333333333333,.1,1,3.667,.9]
chimerisms = [1e-2,2.8740*1e-2,2.875*1e-2,5e-2,1e-1]
chimchoice = [0,2,3,4,7,10,11,12]
initial = (2*k1[0]-1)/k1[2]
t = np.linspace(0,500,500)
chims = np.zeros((2,len(chimerisms),len(t)))
S90 = np.zeros(len(chimerisms))
S90.fill(np.nan)
T90 = np.zeros(len(chimerisms))
T90.fill(np.nan)
colororder = [0,3,4,1,2]

# Simulate steady state
ytemp1 = odeint(Hematopoiesis_2CellJohn,[initial,initial,0,0], np.linspace(0,10000,10000), (k1,l1))[-1]

t = np.linspace(0,400,1000)
y = np.zeros((len(chimerisms),len(t),4))

# Simulate CML
for c in range(len(chimerisms)):
    y[c] = odeint(Hematopoiesis_2CellJohn,[ytemp1[0],ytemp1[1],1.022*ytemp1[0]*chimerisms[c],1.022*ytemp1[1]*chimerisms[c]], t, (k1,l1))

# Plot
with matplotlib.rc_context({'font.size': 20, 'lines.linewidth': 3, 'svg.fonttype': 'none'}):

    labels = ['SP', 'T', 'SP$^\\text{L}$', 'T$^\\text{L}$', 'SP$^\\text{L}$ Chimerism', 'T$^\\text{L}$ Chimerism', '90% S$^\\text{L}$ Chimerism Time', '90% T$^\\text{L}$ Chimerism Time']
    lines = ['-','--']

    f = plt.figure(figsize = (10,5), dpi = 300)
    gs = gridspec.GridSpec(1, 2)
    gs.update(hspace=0.25, wspace = 0.3)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])

    ax = [ax1,ax2]

    for c, chim in enumerate(chimerisms):
        for i in range(2):
            ax[i].semilogy(t/30,100*y[c,:,i+2]/(y[c,:,i]+y[c,:,i+2]), color = 'C{}'.format(colororder[c]), label='Initial Chimerism : {:1.4g}%'.format(100*y[c,0,i+2]/(y[c,0,i]+y[c,0,i+2])), linestyle='solid')
    for i in range(2):
        ax[i].set(ylim = [1e-1,1e2], xlim = [0,10], title=labels[i+4], xlabel='Time After Dox Removal (Months)',ylabel='Chimerism (%)')
        ax[i].yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
        ax[i].set_xlabel('Time (Months)', fontsize=20)
        ax[i].set_ylabel('Chimerism (%)', fontsize=20)
        ax[i].tick_params(axis='both', labelsize=20)
        ax[i].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{:.0f}'.format(x))))
        ax[i].vlines(t[-1]/30,1e-1,1e2,linestyle='dashed',color='k')


lgd = ax[1].legend(bbox_to_anchor=(1.04,.5), loc="center left", borderaxespad=0,fontsize=15)
plt.savefig(os.path.join(OUTPUT_DIR, "SupplementalA.png"), bbox_extra_artists=(lgd,), bbox_inches="tight")
### Treatment
k1 = [0.719673,1,2.4408113329290927*1e-6,0,0,1,1,.03,0,.1,1,1]
k2 = [0.72,1,2.4408113329290927*1e-6*(2*.72-1)/(2*0.719673-1),0,0,1,1,.03,0,.1,1,1]
l1 = [0.719673,1,(2.4408113329290927*1e-6)*(0.8511500334069608),0,0.0016270262087911059,0.4093047,0.04515,0.01633333333333333,.1,1,3.667,.9]
l2 = [0.719673,1,(2.4408113329290927*1e-6)*(0.8511500334069608),0,0.0016270262087911059,0.4093047,0.09,0.01633333333333333,.1,1,3.667,.9]
chimerisms = [1e-2,2.8740*1e-2,2.875*1e-2,5e-2,1e-1]
chimchoice = [0,2,3,4,7,10,11,12]
initial = (2*k1[0]-1)/k1[2]
t = np.linspace(0,500,500)
chims = np.zeros((2,len(chimerisms),len(t)))
S90 = np.zeros(len(chimerisms))
S90.fill(np.nan)
T90 = np.zeros(len(chimerisms))
T90.fill(np.nan)
colororder = [0,3,4,1,2]

ytemp1 = odeint(Hematopoiesis_2CellJohn,[initial,initial,0,0], np.linspace(0,10000,10000), (k1,l1))[-1]
ytemp2 = odeint(Hematopoiesis_2CellJohn,[initial,initial,0,0], t, (k2,l1))[-1]


t1 = np.linspace(0,400,1000)
t2 = np.linspace(0,1000,1001)
t = np.concatenate((t1,t1[-1]+t2))
y1 = np.zeros((len(chimerisms),len(t1),4))
y2 = np.zeros((len(chimerisms),len(t2),4))

for c in range(len(chimerisms)):
    y1[c] = odeint(Hematopoiesis_2CellJohn,[ytemp2[0]*(1-chimerisms[c]),ytemp2[1]*(1-chimerisms[c]),ytemp1[0]*chimerisms[c],ytemp1[1]*chimerisms[c]], t1, (k2,l1))
    y2[c] = odeint(Hematopoiesis_2CellJohn,y1[c,-1], t2, (k2,l2))

y = np.hstack((y1,y2))

with matplotlib.rc_context({'font.size': 20, 'lines.linewidth': 3, 'svg.fonttype': 'none'}):

    labels = ['SP', 'T', 'SP$^\\text{L}$', 'T$^\\text{L}$', 'SP$^\\text{L}$ Chimerism', 'T$^\\text{L}$ Chimerism', '90% S$^\\text{L}$ Chimerism Time', '90% T$^\\text{L}$ Chimerism Time']
    lines = ['-','--']

    f = plt.figure(figsize = (10,5), dpi = 300)
    gs = gridspec.GridSpec(1, 2)
    gs.update(hspace=0.25, wspace = 0.3)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])

    ax = [ax1,ax2]

    for c, chim in enumerate(chimerisms):
        for i in range(2):
            ax[i].semilogy(t/30,100*y[c,:,i+2]/(y[c,:,i]+y[c,:,i+2]), color = 'C{}'.format(colororder[c]), label='Initial Chimerism : {:1.4g}%'.format(100*chim), linestyle='solid')
    for i in range(2):
        ax[i].set(ylim = [1e-1,1e2], xlim = [0,int(t1[-1]/30)+30], title=labels[i+4], xlabel='Time After Dox Removal (Months)',ylabel='Chimerism (%)')
        ax[i].yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
        ax[i].set_xlabel('Time (Months)', fontsize=20)
        ax[i].set_ylabel('Chimerism (%)', fontsize=20)
        ax[i].tick_params(axis='both', labelsize=20)
        ax[i].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{:.0f}'.format(x))))
        ax[i].vlines(t1[-1]/30,1e-1,1e2,linestyle='dashed',color='k')


lgd = ax[1].legend(bbox_to_anchor=(1.04,.5), loc="center left", borderaxespad=0,fontsize=15)
plt.savefig(os.path.join(OUTPUT_DIR, "SupplementalD.png"), bbox_extra_artists=(lgd,), bbox_inches="tight")
## CHIP Mutation Affects Maximum Self-Renewal
k1 = [0.719673,1,2.4408113329290927*1e-6,0,0,1,1,.03,0,.1,1,1]
k2 = [0.72,1,2.4408113329290927*1e-6,0,0,1,1,.03,0,.1,1,1]
l1 = [0.719673,1,(2.4408113329290927*1e-6)*(0.8511500334069608),0,0.0016270262087911059,0.4093047,0.04515,0.01633333333333333,.1,1,3.667,.9]
chimerisms = [1e-2,2.8740*1e-2,2.875*1e-2,5e-2,1e-1]
chimchoice = [0,2,3,4,7,10,11,12]
initial = (2*k1[0]-1)/k1[2]
t = np.linspace(0,500,500)
chims = np.zeros((2,len(chimerisms),len(t)))
S90 = np.zeros(len(chimerisms))
S90.fill(np.nan)
T90 = np.zeros(len(chimerisms))
T90.fill(np.nan)
colororder = [0,3,4,1,2]

ytemp1 = odeint(Hematopoiesis_2CellJohn,[initial,initial,0,0], np.linspace(0,10000,10000), (k1,l1))[-1]

t = np.linspace(0,400,1000)
y = np.zeros((len(chimerisms),len(t),4))

for c in range(len(chimerisms)):
    y[c] = odeint(Hematopoiesis_2CellJohn,[ytemp1[0],ytemp1[1],1.02*ytemp1[0]*chimerisms[c],1.02*ytemp1[1]*chimerisms[c]], t1, (k2,l1))


with matplotlib.rc_context({'font.size': 20, 'lines.linewidth': 3, 'svg.fonttype': 'none'}):

    labels = ['SP', 'T', 'SP$^\\text{L}$', 'T$^\\text{L}$', 'SP$^\\text{L}$ Chimerism', 'T$^\\text{L}$ Chimerism', '90% S$^\\text{L}$ Chimerism Time', '90% T$^\\text{L}$ Chimerism Time']
    lines = ['-','--']

    f = plt.figure(figsize = (10,5), dpi = 300)
    gs = gridspec.GridSpec(1, 2)
    gs.update(hspace=0.25, wspace = 0.3)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])

    ax = [ax1,ax2]

    for c, chim in enumerate(chimerisms):
        for i in range(2):
            ax[i].semilogy(t/30,100*y[c,:,i+2]/(y[c,:,i]+y[c,:,i+2]), color = 'C{}'.format(colororder[c]), label='Initial Chimerism : {:1.4g}%'.format(100*y[c,0,i+2]/(y[c,0,i]+y[c,0,i+2])), linestyle='solid')
    for i in range(2):
        ax[i].set(ylim = [1e-1,1e2], xlim = [0,10], title=labels[i+4], xlabel='Time After Dox Removal (Months)',ylabel='Chimerism (%)')
        ax[i].yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
        ax[i].set_xlabel('Time (Months)', fontsize=20)
        ax[i].set_ylabel('Chimerism (%)', fontsize=20)
        ax[i].tick_params(axis='both', labelsize=20)
        ax[i].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{:.0f}'.format(x))))
        ax[i].vlines(t1[-1]/30,1e-1,1e2,linestyle='dashed',color='k')


lgd = ax[1].legend(bbox_to_anchor=(1.04,.5), loc="center left", borderaxespad=0,fontsize=15)
plt.savefig(os.path.join(OUTPUT_DIR, "SupplementalB.png"), bbox_extra_artists=(lgd,), bbox_inches="tight")
### Treatment
k1 = [0.719673,1,2.4408113329290927*1e-6,0,0,1,1,.03,0,.1,1,1]
k2 = [0.72,1,2.4408113329290927*1e-6,0,0,1,1,.03,0,.1,1,1]
l1 = [0.719673,1,(2.4408113329290927*1e-6)*(0.8511500334069608),0,0.0016270262087911059,0.4093047,0.04515,0.01633333333333333,.1,1,3.667,.9]
l2 = [0.719673,1,(2.4408113329290927*1e-6)*(0.8511500334069608),0,0.0016270262087911059,0.4093047,0.09,0.01633333333333333,.1,1,3.667,.9]
chimerisms = [1e-2,2.8740*1e-2,2.875*1e-2,5e-2,1e-1]
chimchoice = [0,2,3,4,7,10,11,12]
initial = (2*k1[0]-1)/k1[2]
t = np.linspace(0,500,500)
chims = np.zeros((2,len(chimerisms),len(t)))
S90 = np.zeros(len(chimerisms))
S90.fill(np.nan)
T90 = np.zeros(len(chimerisms))
T90.fill(np.nan)
colororder = [0,3,4,1,2]

ytemp1 = odeint(Hematopoiesis_2CellJohn,[initial,initial,0,0], np.linspace(0,10000,10000), (k2,l1))[-1]

t1 = np.linspace(0,400,1000)
t2 = np.linspace(0,1000,1001)
t = np.concatenate((t1,t1[-1]+t2))
y1 = np.zeros((len(chimerisms),len(t1),4))
y2 = np.zeros((len(chimerisms),len(t2),4))

for c in range(len(chimerisms)):
    y1[c] = odeint(Hematopoiesis_2CellJohn,[ytemp1[0],ytemp1[1],1.02*ytemp1[0]*chimerisms[c],1.02*ytemp1[1]*chimerisms[c]], t1, (k2,l1))
    y2[c] = odeint(Hematopoiesis_2CellJohn,y1[c,-1], t2, (k2,l2))

y = np.hstack((y1,y2))

with matplotlib.rc_context({'font.size': 20, 'lines.linewidth': 3, 'svg.fonttype': 'none'}):

    labels = ['SP', 'T', 'SP$^\\text{L}$', 'T$^\\text{L}$', 'SP$^\\text{L}$ Chimerism', 'T$^\\text{L}$ Chimerism', '90% S$^\\text{L}$ Chimerism Time', '90% T$^\\text{L}$ Chimerism Time']
    lines = ['-','--']

    f = plt.figure(figsize = (10,5), dpi = 300)
    gs = gridspec.GridSpec(1, 2)
    gs.update(hspace=0.25, wspace = 0.3)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])

    ax = [ax1,ax2]

    for c, chim in enumerate(chimerisms):
        for i in range(2):
            ax[i].semilogy(t/30,100*y[c,:,i+2]/(y[c,:,i]+y[c,:,i+2]), color = 'C{}'.format(colororder[c]), label='Initial Chimerism : {:1.4g}%'.format(100*y[c,0,i+2]/(y[c,0,i]+y[c,0,i+2])), linestyle='solid')
    for i in range(2):
        ax[i].set(ylim = [1e-1,1e2], xlim = [0,int(t1[-1]/30)+30], title=labels[i+4], xlabel='Time After Dox Removal (Months)',ylabel='Chimerism (%)')
        ax[i].yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
        ax[i].set_xlabel('Time (Months)', fontsize=20)
        ax[i].set_ylabel('Chimerism (%)', fontsize=20)
        ax[i].tick_params(axis='both', labelsize=20)
        ax[i].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{:.0f}'.format(x))))
        ax[i].vlines(t1[-1]/30,1e-1,1e2,linestyle='dashed',color='k')


lgd = ax[1].legend(bbox_to_anchor=(1.04,.5), loc="center left", borderaxespad=0,fontsize=15)
plt.savefig(os.path.join(OUTPUT_DIR, "SupplementalE.png"), bbox_extra_artists=(lgd,), bbox_inches="tight")
## CHIP Mutation Affects Self-Renewal Feedback Gain
k1 = [0.719673,1,2.4408113329290927*1e-6,0,0,1,1,.03,0,.1,1,1]
k2 = [0.719673,1,0.99*2.4408113329290927*1e-6,0,0,1,1,.03,0,.1,1,1]
l1 = [0.719673,1,(2.4408113329290927*1e-6)*(0.8511500334069608),0,0.0016270262087911059,0.4093047,0.04515,0.01633333333333333,.1,1,3.667,.9]
l2 = [0.719673,1,(2.4408113329290927*1e-6)*(0.8511500334069608),0,0.0016270262087911059,0.4093047,0.09,0.01633333333333333,.1,1,3.667,.9]
chimerisms = np.concatenate((np.linspace(1e-2,1e-1,10),[2.8740*1e-2,2.875*1e-2],np.linspace(2e-1,1e0,9)))
chimerisms = [1e-2,2.8740*1e-2,2.875*1e-2,5e-2,1e-1]
chimchoice = [0,2,3,4,7,10,11,12]
initial = (2*k1[0]-1)/k1[2]
t = np.linspace(0,500,500)
chims = np.zeros((2,len(chimerisms),len(t)))
S90 = np.zeros(len(chimerisms))
S90.fill(np.nan)
T90 = np.zeros(len(chimerisms))
T90.fill(np.nan)
colororder = [0,3,4,1,2]

ytemp1 = odeint(Hematopoiesis_2CellJohn,[initial,initial,0,0], np.linspace(0,10000,10000), (k1,l1))[-1]

t = np.linspace(0,400,1000)
y = np.zeros((len(chimerisms),len(t),4))

for c in range(len(chimerisms)):
    y[c] = odeint(Hematopoiesis_2CellJohn,[ytemp1[0],ytemp1[1],1.01*ytemp1[0]*chimerisms[c],1.01*ytemp1[1]*chimerisms[c]], t, (k2,l1))


with matplotlib.rc_context({'font.size': 20, 'lines.linewidth': 3, 'svg.fonttype': 'none'}):

    labels = ['SP', 'T', 'SP$^\\text{L}$', 'T$^\\text{L}$', 'SP$^\\text{L}$ Chimerism', 'T$^\\text{L}$ Chimerism', '90% S$^\\text{L}$ Chimerism Time', '90% T$^\\text{L}$ Chimerism Time']
    lines = ['-','--']

    f = plt.figure(figsize = (10,5), dpi = 300)
    gs = gridspec.GridSpec(1, 2)
    gs.update(hspace=0.25, wspace = 0.3)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])

    ax = [ax1,ax2]

    for c, chim in enumerate(chimerisms):
        for i in range(2):
            ax[i].semilogy(t/30,100*y[c,:,i+2]/(y[c,:,i]+y[c,:,i+2]), color = 'C{}'.format(colororder[c]), label='Initial Chimerism : {:1.4g}%'.format(100*y[c,0,i+2]/(y[c,0,i]+y[c,0,i+2])), linestyle='solid')
    for i in range(2):
        ax[i].set(ylim = [1e-1,1e2], xlim = [0,10], title=labels[i+4], xlabel='Time After Dox Removal (Months)',ylabel='Chimerism (%)')
        ax[i].yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
        ax[i].set_xlabel('Time (Months)', fontsize=20)
        ax[i].set_ylabel('Chimerism (%)', fontsize=20)
        ax[i].tick_params(axis='both', labelsize=20)
        ax[i].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{:.0f}'.format(x))))
        ax[i].vlines(t[-1]/30,1e-1,1e2,linestyle='dashed',color='k')


lgd = ax[1].legend(bbox_to_anchor=(1.04,.5), loc="center left", borderaxespad=0,fontsize=15)
plt.savefig(os.path.join(OUTPUT_DIR, "SupplementalC.png"), bbox_extra_artists=(lgd,), bbox_inches="tight")
### Treatment
k1 = [0.719673,1,2.4408113329290927*1e-6,0,0,1,1,.03,0,.1,1,1]
k2 = [0.719673,1,0.99*2.4408113329290927*1e-6,0,0,1,1,.03,0,.1,1,1]
l1 = [0.719673,1,(2.4408113329290927*1e-6)*(0.8511500334069608),0,0.0016270262087911059,0.4093047,0.04515,0.01633333333333333,.1,1,3.667,.9]
l2 = [0.719673,1,(2.4408113329290927*1e-6)*(0.8511500334069608),0,0.0016270262087911059,0.4093047,0.09,0.01633333333333333,.1,1,3.667,.9]
chimerisms = [1e-2,2.8740*1e-2,2.875*1e-2,5e-2,1e-1]
chimchoice = [0,2,3,4,7,10,11,12]
initial = (2*k1[0]-1)/k1[2]
t = np.linspace(0,500,500)
chims = np.zeros((2,len(chimerisms),len(t)))
S90 = np.zeros(len(chimerisms))
S90.fill(np.nan)
T90 = np.zeros(len(chimerisms))
T90.fill(np.nan)
colororder = [0,3,4,1,2]

ytemp1 = odeint(Hematopoiesis_2CellJohn,[initial,initial,0,0], np.linspace(0,10000,10000), (k2,l1))[-1]

t1 = np.linspace(0,400,1000)
t2 = np.linspace(0,1000,1001)
t = np.concatenate((t1,t1[-1]+t2))
y1 = np.zeros((len(chimerisms),len(t1),4))
y2 = np.zeros((len(chimerisms),len(t2),4))

for c in range(len(chimerisms)):
    y1[c] = odeint(Hematopoiesis_2CellJohn,[ytemp1[0],ytemp1[1],1.01*ytemp1[0]*chimerisms[c],1.01*ytemp1[1]*chimerisms[c]], t1, (k2,l1))
    y2[c] = odeint(Hematopoiesis_2CellJohn,y1[c,-1], t2, (k2,l2))

y = np.hstack((y1,y2))

with matplotlib.rc_context({'font.size': 20, 'lines.linewidth': 3, 'svg.fonttype': 'none'}):

    labels = ['SP', 'T', 'SP$^\\text{L}$', 'T$^\\text{L}$', 'SP$^\\text{L}$ Chimerism', 'T$^\\text{L}$ Chimerism', '90% S$^\\text{L}$ Chimerism Time', '90% T$^\\text{L}$ Chimerism Time']
    lines = ['-','--']

    f = plt.figure(figsize = (10,5), dpi = 300)
    gs = gridspec.GridSpec(1, 2)
    gs.update(hspace=0.25, wspace = 0.3)

    ax1 = plt.subplot(gs[0, 0])
    ax2 = plt.subplot(gs[0, 1])

    ax = [ax1,ax2]

    for c, chim in enumerate(chimerisms):
        for i in range(2):
            ax[i].semilogy(t/30,100*y[c,:,i+2]/(y[c,:,i]+y[c,:,i+2]), color = 'C{}'.format(colororder[c]), label='Initial Chimerism : {:1.4g}%'.format(100*y[c,0,i+2]/(y[c,0,i]+y[c,0,i+2])), linestyle='solid')
    for i in range(2):
        ax[i].set(ylim = [1e-1,1e2], xlim = [0,int(t1[-1]/30)+30], title=labels[i+4], xlabel='Time After Dox Removal (Months)',ylabel='Chimerism (%)')
        ax[i].yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'.format(int(np.maximum(-np.log10(y),0)))).format(y)))
        ax[i].set_xlabel('Time (Months)', fontsize=20)
        ax[i].set_ylabel('Chimerism (%)', fontsize=20)
        ax[i].tick_params(axis='both', labelsize=20)
        ax[i].xaxis.set_major_formatter(ticker.FuncFormatter(lambda x,pos: ('{:.0f}'.format(x))))
        ax[i].vlines(t1[-1]/30,1e-1,1e2,linestyle='dashed',color='k')


lgd = ax[1].legend(bbox_to_anchor=(1.04,.5), loc="center left", borderaxespad=0,fontsize=15)
plt.savefig(os.path.join(OUTPUT_DIR, "SupplementalF.png"), bbox_extra_artists=(lgd,), bbox_inches="tight")
