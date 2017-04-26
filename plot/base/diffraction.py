import matplotlib as mpl
#mpl.rc('text',usetex=True)
#mpl.rcParams['legend.numpoints']=1
#mpl.rcParams['font.size'] = 27
#mpl.rcParams['font.weight']   = 'bold'
#mpl.rcParams['text.latex.preamble']=[r'\usepackage{bm} \boldmath']
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import ROOT

def diffraction(filepath):

    #f09  = ROOT.TFile(filepath+'900_4M.root')
    #f13  = ROOT.TFile(filepath+'13000_4M.root')
    files = [ROOT.TFile(filepath+'900_4M.root'),ROOT.TFile(filepath+'13000_4M.root')]
    objects = ['DOU_NF_0{}','SIN_NF_0{}','DOU_NBNF_0{}','SIN_NBNF_0{}']

    fi = []; NF = []; NFx = []
    figs = []; axs = []; 
    for k,f in enumerate(files):
        for i,obj in enumerate(objects):

            fi.append(f.FindObjectAny(obj.format(1)))
            i = i if k<1 else i+len(objects)

            for j in range(2,len(objects)):
                fi[i].Add(f.FindObjectAny(obj.format(j)))
            
            Nbins = fi[i].GetNbinsX()

            NF.append(np.asarray([fi[i].GetBinContent(j) for j in range(1,Nbins)]))

            NF[i] = np.trim_zeros(NF[i],trim='b')

            limit = [-10,10] if i < 2 else [0,100]
            NFx.append(np.linspace(limit[0],limit[1],len(NF[i])))

            fig, ax = plt.subplots()
            axs.append(ax)
            axs[i].plot(NFx[i],NF[i],linestyle='',marker='o',label='{}'.format(i))
            axs[i].legend()

    plt.show()

if __name__=='__main__':
    diffraction('/home/roar/master/qgsm_analysis_tool/ana/build/')
