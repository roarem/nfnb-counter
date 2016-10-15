import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import ROOT
from scipy.interpolate import spline

class Plot:

    def __init__(self,root_file_path=None):
        try:
            self.f = ROOT.TFile(root_file_path)
        except:
            print("Could not find root file")
        
        #self.npoms1h0 = ROOT.gROOT.FindObject("NPOM_01_00")
        #self.npoms1h0_nf = ROOT.gROOT.FindObject("NPOM_NF_01_00")
        #self.npoms1h0.Divide(self.npoms1h0_nf)
        #self.npoms1h0.Draw()
        #print(self.npoms1h0.GetNbinsX())

    def __call__(self):
        pass

    def NPOM(self):
        majorLocator = MultipleLocator(10)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(1)

        expNF,expNBNF,expXerr,expYerr = np.loadtxt('../out/nbnf_7000_exp',unpack=True)
        simNBNF_in  = ROOT.gROOT.FindObject("ptcut_div")
        Nbins       = simNBNF_in.GetNbinsX()
        simNF       = np.linspace(0,Nbins,Nbins-1)
        simNBNF     = np.asarray([simNBNF_in.GetBinContent(i) for i in range(1,Nbins)])
        simNBNF_err = np.asarray([simNBNF_in.GetBinError(i) for i in range(1,Nbins)])
            
        NPOMList = []
        for i in range(7):
            SxH0_name = "NPOM_%02d_00"%i
            SxH0_nf_name = "NPOM_NF_%02d_00"%i
            SxH0 = ROOT.gROOT.FindObject(SxH0_name)
            SxH0_nf = ROOT.gROOT.FindObject(SxH0_nf_name)
            SxH0.Divide(SxH0_nf)
            Nbins = SxH0.GetNbinsX()
            temp_hist = [SxH0.GetBinContent(i) for i in range(1,Nbins)]
            NPOMList.append(np.asarray(temp_hist))
            
        fig, ax = plt.subplots()
        nf = np.linspace(0,30,31)
        nfnew = np.linspace(nf.min(),nf.max(),300)
        graphs = []
        for i,npom in enumerate(NPOMList):
            npom = np.trim_zeros(npom,trim='b')[:11+i]
            smooth = spline(nf[:len(npom)],npom,nfnew)
            #ax.plot(nf[:len(npom)],npom,marker='',linestyle='-',label='NPOMS {}'.format(i))
            smooth = np.trim_zeros(smooth,trim='b')
            ax.plot(nfnew[:len(smooth)],smooth,linestyle='-',linewidth=1.5,\
                    color='grey',label='NPOMS {}'.format(i))

        ax.errorbar(expNF,expNBNF,yerr=expYerr,linestyle='',\
                    marker='o',markersize=8,color='black',\
                    label='experimental')
        ax.errorbar(simNF,simNBNF,yerr=simNBNF_err,linestyle='',\
                    marker='s',markersize=8,color='black',label='simulation')

        fontsize_labels = 24
        ax.set_xlim(-1,25)
        ax.set_ylim(0,15)
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))

        [ax.text(-0.4,npom[0]-0.1,str(i),fontsize=16) for i,npom in enumerate(NPOMList)]

        ax.set_ylabel('$<n_B(n_F)>$',fontsize=fontsize_labels)
        ax.set_xlabel('$n_F$',fontsize=fontsize_labels)
        ax.grid(which='minor',alpha=0.5)
        #ax.yaxis.set_major_locator(majorLocator)
        #ax.yaxis.set_major_formatter(majorFormatter)
        ax.yaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)
        [tick.label.set_fontsize(14) for tick in ax.xaxis.get_major_ticks()]
        [tick.label.set_fontsize(14) for tick in ax.yaxis.get_major_ticks()]
        ax.xaxis.set_tick_params(which='major',size=12)
        ax.yaxis.set_tick_params(which='major',size=12)
        ax.xaxis.set_tick_params(which='minor',size=8)
        ax.yaxis.set_tick_params(which='minor',size=8)
        handles, labels = ax.get_legend_handles_labels()    

        leg = plt.legend((handles[7],handles[8],handles[0]),\
                         ('ALICE pp, $\sqrt{s}$=7 TeV',\
                          'QGSM $0.3<p_T<1.5$ GeV/c, $0.2<\eta<0.8$',\
                          'NPOMH = 0 and for NPOMS = {0,...,6}'),\
                         loc='best',fontsize=24)

        leg.get_frame().set_alpha(0.0)
        plt.show()


if __name__=="__main__":
    P = Plot(root_file_path="../out/test.root")
    P.NPOM()
