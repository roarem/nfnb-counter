import matplotlib as mpl
mpl.rc('text',usetex=True)
mpl.rcParams['legend.numpoints']=1
mpl.rcParams['font.size'] = 27
mpl.rcParams['font.weight']   = 'bold'
mpl.rcParams['text.latex.preamble']=[r'\usepackage{bm} \boldmath']
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import ROOT

def diffraction(filepath):

    f = ROOT.TFile(filepath+'900_4M.root')
    nf_sin = f.FindObjectAny('DOU_NF_01')
    #nb_sin = f.FindObjectAny('DOU_NB_01')
    for i in [2,3]:
        nf_sin.Add(f.FindObjectAny('DOU_NF_0{}'.format(i)))
        #nb_sin.Add(f.FindObjectAny('DOU_NB_0{}'.format(i)))
    Nbins = nf_sin.GetNbinsX()
    nf = np.asarray([nf_sin.GetBinContent(i) for i in range(1,Nbins)])

    x = np.linspace(-10,10,Nbins-1)
    plt.plot(x,nf,linestyle='',marker='o')
    plt.show()


if __name__=='__main__':
    diffraction('/home/roar/master/qgsm_analysis_tool/ana/build/')
