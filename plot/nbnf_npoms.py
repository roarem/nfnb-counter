import matplotlib as mpl
#mpl.use('Agg')
mpl.rc('text',usetex=True)
mpl.rcParams['legend.numpoints']=1
mpl.rcParams['font.size'] = 27
mpl.rcParams['font.weight']   = 'bold'
mpl.rcParams['text.latex.preamble']=[r'\usepackage{bm} \boldmath']
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import ROOT
import linecache
from scipy.interpolate import spline
from scipy.stats import linregress

def plotSetup(xlim=None,ylim=None):
    
    majorLocator_y = ticker.MultipleLocator(2)
    majorLocator_x = ticker.MultipleLocator(5)
    minorLocator = ticker.MultipleLocator(1)

    fig, axarr = plt.subplots(2,2)#,sharex='col',sharey='row')
    axarr = np.reshape(axarr,4)
    DPI = fig.get_dpi()
    size = 1000
    fig.set_size_inches(size/DPI,size/DPI)
    
    for ax in axarr:
        ax.set_xlim(xlim[0],xlim[1])
        ax.set_ylim(ylim[0],ylim[1])
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))

        ax.grid(which='minor',alpha=0.5)

        majorFormatter = ticker.FormatStrFormatter('%d')
        minorFormatter = ticker.FormatStrFormatter('%d')
        ax.yaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_minor_locator(minorLocator)
        ax.yaxis.set_major_locator(majorLocator_y)
        ax.xaxis.set_major_locator(majorLocator_x)
        ax.xaxis.set_major_formatter(majorFormatter)
        #ax.xaxis.set_minor_formatter(minorFormatter)

        [tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
        [tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]

        ax.xaxis.set_tick_params(which='major',length=12,width=4)
        ax.yaxis.set_tick_params(which='major',length=12,width=4)
        ax.xaxis.set_tick_params(which='minor',length=8 ,width=2)
        ax.yaxis.set_tick_params(which='minor',length=8 ,width=2)

    return fig,axarr

def nbnf_npoms(fig,axarr):
    filepath    = "/home/roar/master/qgsm_analysis_tool/ana/"
    filenames   = [\
                   '/home/roar/master/qgsm_analysis_tool/qgsmpy/src/900_4m_test.root',\
                   #filepath+'out/2212/900_4M.root',\
                   filepath+'out/2212/2760_4M.root',\
                   filepath+'out/2212/7000_4M.root',\
                   filepath+'out/2212/13000_4M.root']
    titles      = ['0.9 TeV','2.76 TeV','7 TeV','13 TeV']
    prefix      = 'ptcut'
    N_bins      = 9
    for i_ax,(name,ax) in enumerate(zip(filenames,axarr)):
        f = ROOT.TFile(name)
        npom_list = []
        for i in range(1,7):
            if i!=1:
                SnH0    .Add(f.FindObjectAny(prefix+"_NPOM_{:02d}_00".format(i)))
                SnH0_nf .Add(f.FindObjectAny(prefix+"_NPOM_NF_{:02d}_00".format(i)))
            else:
                SnH0    = f.FindObjectAny(prefix+"_NPOM_{:02d}_00".format(i))
                SnH0_nf = f.FindObjectAny(prefix+"_NPOM_NF_{:02d}_00".format(i))
        #N_bins = SnH0.GetNbinsX()
        SnH0.Divide(SnH0_nf) 
        SnH0_nbnf       = np.asarray(\
                          [SnH0.GetBinContent(k) for k in range(1,N_bins)])#N_bins)])
        SnH0_nbnf_err   = np.asarray(\
                          [SnH0.GetBinContent(k) for k in range(1,N_bins)])#N_bins)])
        SnH0_nf         = np.linspace(0,N_bins,N_bins-1)#N_bins,N_bins-1)

        SnH0_nf_new = np.linspace(SnH0_nf.min(),SnH0_nf.max(),300)
        smooth = spline(SnH0_nf,SnH0_nbnf,SnH0_nf_new)
        ax.plot(SnH0_nf_new,smooth,linestyle='-',linewidth=4,color='black',\
                label='all npoms',zorder=7)
        #jpolynome = np.polyfit(SnH0_nf[:6],SnH0_nbnf[:6],1)
        polynome = np.polyfit(SnH0_nf,SnH0_nbnf,1)
        print(name)
        print(polynome)
        f.Close()
        f = ROOT.TFile(name)
        for i in range(1,7):
            SxH0    = f.FindObjectAny(prefix+"_NPOM_{:02d}_00".format(i))
            SxH0_nf = f.FindObjectAny(prefix+"_NPOM_NF_{:02d}_00".format(i))
            SxH0.Divide(SxH0_nf)
            temp_hist = [SxH0.GetBinContent(j) for j in range(1,N_bins)]#N_bins)]
            npom_list.append(np.asarray(temp_hist))

        #nf = np.linspace(0,30,31)
        #nf = np.linspace(0,30,31)
        #nfnew = np.linspace(nf.min(),nf.max(),300)
        nfnew = np.linspace(SnH0_nf.min(),SnH0_nf.max(),300)
        for i,npom in enumerate(npom_list):
            #npom = np.trim_zeros(npom,trim='b')[:11+i]
            smooth = spline(SnH0_nf[:len(npom)],npom,nfnew)
            smooth = np.trim_zeros(smooth,trim='b')
            ax.plot(nfnew[:len(smooth)],smooth,linestyle='--',linewidth=4,\
                    color='black',label='NPOMS {}'.format(i),zorder=6)
        
        [ax.text(-1,npom[0]-0.1,str(i+1),fontsize=24) for i,npom in enumerate(npom_list)]
	ax.text(0.90, 0.05, titles[i_ax],\
	        verticalalignment='bottom', horizontalalignment='right',\
	        transform=ax.transAxes,\
	        color='black', fontsize=27)
        

        plt.setp([a.get_xticklabels() for a in [axarr[0],axarr[1]]], visible=False)
        plt.setp([a.get_yticklabels() for a in [axarr[1],axarr[3]]], visible=False)
        fig.subplots_adjust(hspace=0.01)
        fig.subplots_adjust(wspace=-0.6)
        fig.text(0.5, 0.04, '$n_F$', ha='center', va='center')
        fig.text(0.27, 0.5, '$\langle n_B(n_F)\\rangle$', ha='center', va='center', rotation='vertical')
        f.Close()

    handles, labels = axarr[0].get_legend_handles_labels()    
    leg = axarr[0].legend((handles[0],handles[1]),\
                     ('soft pomerons = 1-6',\
                      '# soft pomerons = \\{1,...,6 \\}'),\
                     loc=0,fontsize=24)

if __name__=='__main__':
    fig, axarr = plotSetup(xlim=[-2,10],ylim=[0,7])
    
    nbnf_npoms(fig,axarr)
    plt.show()

