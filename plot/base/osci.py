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


def osci_search(filepath):

    fig, ax = plt.subplots()
    f = ROOT.TFile(filepath)
    markers = ['','1','2','3','4','_','d','*','o']
    for i in range(1,9):
        tempS_nbnf = f.FindObjectAny("SIN_NBNF_rap_le_0{}".format(i))
        tempS_nf   = f.FindObjectAny("SIN_NF_rap_le_0{}".format(i))
        #tempS_nbnf .Add(f.FindObjectAny("DOU_NBNF_rap_le_0{}".format(i)))
        #tempS_nf   .Add(f.FindObjectAny("DOU_NF_rap_le_0{}".format(i)))

        Nbins = tempS_nbnf.GetNbinsX()
        tempS_nbnf.Divide(tempS_nf)
        nbnf_nf = np.asarray([tempS_nbnf.GetBinContent(j) for j in range(1,Nbins)])
        nbnf_nf[nbnf_nf==0] = np.nan
        nf = np.linspace(0,Nbins,Nbins-1) 
        ax.plot(nf,nbnf_nf,marker='',markersize=10,
                linestyle='-',label=r'$<{:d}$'.format(i))
        print(Nbins)
    #ax.set_xlim(0,30)
    plt.legend(loc='best')
    plt.show()


if __name__=='__main__':
    path = '/home/roar/master/qgsm_analysis_tool/ana/b/900_500k.root'
    osci_search(path)
