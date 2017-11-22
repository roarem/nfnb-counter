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

def plotSetup(xlim=None,ylim=None):
    
    majorLocator = ticker.MultipleLocator(5)
    minorLocator = ticker.MultipleLocator(1)

    fig, axarr = plt.subplots(1,2)
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
        ax.yaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        #ax.xaxis.set_minor_formatter(minorFormatter)

        [tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
        [tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]

        ax.xaxis.set_tick_params(which='major',length=12,width=4)
        ax.yaxis.set_tick_params(which='major',length=12,width=4)
        ax.xaxis.set_tick_params(which='minor',length=8 ,width=2)
        ax.yaxis.set_tick_params(which='minor',length=8 ,width=2)

    return fig,axarr

def nbnf_picture():
    filepath = '/home/roar/master/qgsm_analysis_tool/ana/'
    fig, axarr = plotSetup(xlim=[-1,22],ylim=[0,12])
    xsize = int(axarr[0].get_xlim()[1])
######################### Larissas points ##################################
    #larNF,larNBNF = np.loadtxt('/home/roar/master/qgsm_analysis_tool/ana/out/larissa_anfnb_exp_soft',\
    #        unpack=True)
    #ax.plot(larNF,larNBNF,label='Larissa NBNF')
######################### Experimental data ################################
    expNF,expNBNF,expXerr,expYerr = np.loadtxt(filepath+'out/nbnf_7000_exp',unpack=True)
    axarr[0].errorbar(expNF,expNBNF,yerr=expYerr,linestyle='--',\
                marker='o',markersize=15,color='black',\
                label='ALICE 7 TeV',zorder=10)

    f = ROOT.TFile(filepath+'out/2112/7000_4M.root')
    prefix = 'ptcut'
    for i in range(25):
        for j in range(25):
            if i or j:
                tempS_nbnf.Add(f.FindObjectAny(prefix+"_NPOM_{:02d}_{:02d}".format(i,j)))
                tempS_nf.Add(f.FindObjectAny(prefix+"_NPOM_NF_{:02d}_{:02d}".format(i,j)))
            else:
                tempS_nbnf = f.FindObjectAny(prefix+"_NPOM_{:02d}_{:02d}".format(i,j))
                tempS_nf = f.FindObjectAny(prefix+"_NPOM_NF_{:02d}_{:02d}".format(i,j))

    Nbins  = tempS_nbnf.GetNbinsX()
    tempS_nbnf.Divide(tempS_nf)
    simNBNF     = np.asarray([tempS_nbnf.GetBinContent(k) for k in range(1,Nbins)])
    simNBNF_err = np.asarray([tempS_nbnf.GetBinError(k) for k in range(1,Nbins)])
    simNF       = np.linspace(0,Nbins,(Nbins+1))
    axarr[0].errorbar(simNF[:xsize],simNBNF[:xsize],yerr=simNBNF_err[:xsize],\
                linestyle='--', marker='s',markersize=15,color='black',\
                label='QGSM 7 TeV'),\

    handles = [mpl.lines.Line2D([],[],color='black',marker='o',markersize=15,\
                                linestyle='',label='ALICE 7 TeV'),\
               mpl.lines.Line2D([],[],color='black',marker='s',\
                                        markersize=15,linestyle='',\
                                        label='QGSM 7 TeV')]
    axarr[0].legend(handles=handles,loc='best')
    f.Close()
    handles = []
    for energy,marker in zip([0.9,2.76,7,13],['*','D','s','^']):
        f = ROOT.TFile(filepath+'out/2112/{}_4M.root'.format(str(int(energy*1000))))
        for i in range(25):
            for j in range(25):
                if i or j:
                    tempS_nbnf.Add(f.FindObjectAny(prefix+"_NPOM_{:02d}_{:02d}".format(i,j)))
                    tempS_nf.Add(f.FindObjectAny(prefix+"_NPOM_NF_{:02d}_{:02d}".format(i,j)))
                else:
                    tempS_nbnf = f.FindObjectAny(prefix+"_NPOM_{:02d}_{:02d}".format(i,j))
                    tempS_nf = f.FindObjectAny(prefix+"_NPOM_NF_{:02d}_{:02d}".format(i,j))

        Nbins  = tempS_nbnf.GetNbinsX()
        tempS_nbnf.Divide(tempS_nf)
        simNBNF     = np.asarray([tempS_nbnf.GetBinContent(k) for k in range(1,Nbins)])
        simNBNF_err = np.asarray([tempS_nbnf.GetBinError(k) for k in range(1,Nbins)])
        simNF       = np.linspace(0,Nbins,(Nbins+1))
        axarr[1].errorbar(simNF[:xsize],simNBNF[:xsize],yerr=simNBNF_err[:xsize],\
                    linestyle='--', marker=marker,markersize=15,color='black',\
                    label='QGSM '+str(energy)+' TeV'),\
                    #zorder=9)
        handles.append(mpl.lines.Line2D([],[],color='black',marker=marker,\
                                        markersize=15,linestyle='',\
                                        label='{} TeV'.format(str(energy))))
        f.Close()

    
    #handles, labels = ax.get_legend_handles_labels()    
    axarr[1].legend(handles=handles,loc='best')
    fontdict = {'fontsize':27,'weight':'bold'}
    axarr[0].set_ylabel(r'$\langle n_B(n_F)\rangle$',fontdict=fontdict)
    #axarr[1].set_xlabel('$n_F$',fontdict=fontdict)
    plt.setp(axarr[1].get_yticklabels(), visible=False)
    fig.subplots_adjust(wspace=0.01)
    fig.text(0.51, 0.08, '$n_F$', ha='center', va='center')
    if 0:
        simNBNF_fit = tempS_nbnf.Fit('pol1','SQN','',0,30)
        fit_a = simNBNF_fit.Parameter(0)
        fit_b = simNBNF_fit.Parameter(1)
        nf = np.linspace(0,30,31)
        ax.plot(nf,self.fit(fit_a,fit_b,nf),linestyle='-',marker='',markersize=8,color='black',\
                linewidth=3,label='simNBNF self.fit',zorder=8)


##################### General settings ################################################


    #leg = plt.legend((handles[8],handles[9],handles[0],handles[7]),#handles[10]),\
    #                 (\
    #                  'ALICE pp, $\sqrt{s}$=7 TeV',\
    #                  r'QGSM $0.3<p_T<1.5$ GeV/c, $0.2<|\eta|<0.8$',\
    #                  'NPOMH = 0 and for NPOMS = \{0,...,6\}',\
    #                  'All NPOMS and NPOMH = 0'\
    #                  ),\
    #                  #'{:.3f}+{:.3f}x self.fit of simulated data'.format(self.fit_a,self.fit_b),\
    #                  loc='best',fontsize=24)
    #leg.get_frame().set_alpha(0.0)
if __name__=='__main__':
    nbnf_picture()
    plt.show()
