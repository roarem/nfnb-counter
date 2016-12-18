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


def nbnf_w_kin(filepath):

    def fit(a,b,x):
        return a+b*x
    
    fig, axs = plt.subplots(2,2)
    axs = axs.reshape(4)
    fig.subplots_adjust(wspace=0.001,hspace=0.001)
    DPI = fig.get_dpi()
    fig.set_size_inches(1000/DPI,1000/DPI)

    m  = 'o'
    mc = 'black'
    ms = 10
    ls = ''
    lc = 'red'
    filenames = ["900_4M.root","2760_4M.root","7000_4M.root","13000_4M.root"]
    energies  = [900,2760,7000,13000]

    for e,filename in enumerate(filenames):

        f = ROOT.TFile(filepath+filename)
        for i in range(25):
            for j in range(25):
                if i or j:
                    tempS_nbnf.Add(f.FindObjectAny("ptcut_NPOM_{:02d}_{:02d}".format(i,j)))
                    tempS_nf.Add(f.FindObjectAny("ptcut_NPOM_NF_{:02d}_{:02d}".format(i,j)))
                else:
                    tempS_nbnf = f.FindObjectAny("ptcut_NPOM_{:02d}_{:02d}".format(i,j))
                    tempS_nf = f.FindObjectAny("ptcut_NPOM_NF_{:02d}_{:02d}".format(i,j))

        
        Nbins = 35#tempS_nbnf.GetNbinsX()
        tempS_nbnf.Divide(tempS_nf)
        nbnf_nf = np.asarray([tempS_nbnf.GetBinContent(i) for i in range(1,Nbins)])
        nbnf_nf[nbnf_nf==0] = np.nan
        nbnf_nf_err = np.asarray([tempS_nbnf.GetBinError(i) for i in range(1,Nbins)])
        fitting = tempS_nbnf.Fit('pol1','SQN','',0,Nbins)
        a = fitting.Parameter(0)
        b = fitting.Parameter(1)
        nbnf_nf_fit = fit(a,b,np.linspace(0,Nbins,Nbins-1))
        nf = np.linspace(0,Nbins,Nbins-1) 

        axs[e].errorbar(nf,nbnf_nf,nbnf_nf_err,
                        marker=m,linestyle=ls,markersize=ms,color=mc)
        axs[e].plot(nf,nbnf_nf_fit,color=lc)
        


        f.Close()

    [ax.grid(which='major',alpha=0.5) for ax in axs]

    for ax in axs:
        ax.set_ylim(0,25)
        ax.xaxis.set_tick_params(which='major',length=14,width=2)
        ax.yaxis.set_tick_params(which='major',length=14,width=2)
        ax.xaxis.set_tick_params(which='minor',length=8 ,width=2)
        ax.yaxis.set_tick_params(which='minor',length=8 ,width=2)
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
        ax.yaxis.set_minor_locator(ticker.MultipleLocator(1))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(5))

    #axs[0].set_xticks([])
    #axs[1].set_xticks([])
    #axs[1].set_yticks([])
    #axs[3].set_yticks([])
    [label.set_visible(False) for label in axs[0].get_xticklabels()]
    [label.set_visible(False) for label in axs[1].get_xticklabels()]
    [label.set_visible(False) for label in axs[1].get_yticklabels()]
    [label.set_visible(False) for label in axs[3].get_yticklabels()]
    [axs[2].get_yticklabels()[-i].set_visible(False) for i in range(3)]
    [axs[2].get_xticklabels()[-i].set_visible(False) for i in range(3)]

    [ax.text(0.1,0.9,str(ener),horizontalalignment='center',
             verticalalignment='center',transform=ax.transAxes) for ax,ener in zip(axs,energies)]

    plt.show()
    

if __name__=='__main__':
    path = '/home/roar/master/qgsm_analysis_tool/ana/out/1612/'
    nbnf_w_kin(path)
