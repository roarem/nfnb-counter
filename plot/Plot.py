import matplotlib as mpl
mpl.rc('text',usetex=True)
mpl.rcParams['text.latex.preamble']=[r'\usepackage{bm} \boldmath']
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import ROOT
from scipy.interpolate import spline

class Plot:

    def __init__(self,root_file_path=None):
        self.filepath = root_file_path
        try:
            self.f = ROOT.TFile(self.filepath)
        except:
            print("Could not find root file")

    def __call__(self):
        pass

    def ReOpen(self):
        self.f.Close()
        self.f = ROOT.TFile(self.filepath)

    def Close(self):
        self.f.Close()

    def Show(self):
        plt.show()
    
    def plotSetup(self,xlim=None,ylim=None):

        majorLocator = MultipleLocator(10)
        majorFormatter = FormatStrFormatter('%d')
        minorLocator = MultipleLocator(1)
        fig, ax = plt.subplots()
        DPI = fig.get_dpi()
        size = 1000
        fig.set_size_inches(size/DPI,size/DPI)
        ax.grid(which='minor',alpha=0.5)

        ax.set_xlim(xlim[0],xlim[1])
        ax.set_ylim(ylim[0],ylim[1])
        x0,x1 = ax.get_xlim()
        y0,y1 = ax.get_ylim()
        ax.set_aspect((x1-x0)/(y1-y0))

        ax.yaxis.set_minor_locator(minorLocator)
        ax.xaxis.set_major_locator(majorLocator)
        ax.xaxis.set_major_formatter(majorFormatter)
        ax.xaxis.set_minor_locator(minorLocator)

        [tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
        [tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]

        ax.xaxis.set_tick_params(which='major',length=12,width=4)
        ax.yaxis.set_tick_params(which='major',length=12,width=4)
        ax.xaxis.set_tick_params(which='minor',length=8 ,width=2)
        ax.yaxis.set_tick_params(which='minor',length=8 ,width=2)

        return fig,ax
    
    def fit(self,a,b,x):
        return a+b*x

    def NBNFPicture(self):

        fig, ax = self.plotSetup(xlim=[-1,25],ylim=[0,15])#plt.subplots()
        
############################# Experimental data ################################
        expNF,expNBNF,expXerr,expYerr = np.loadtxt('/home/roar/master/qgsm_analysis_tool/ana/out/nbnf_7000_exp',unpack=True)
        ax.errorbar(expNF,expNBNF,yerr=expYerr,linestyle='',\
                    marker='o',markersize=10,color='black',\
                    label='experimental',zorder=10)

############################# Simulated data ###################################
        simNBNF_in  = ROOT.gROOT.FindObject("ptcut_div")
        Nbins       = simNBNF_in.GetNbinsX()
        simNF       = np.linspace(0,Nbins,Nbins-1)
        simNBNF     = np.asarray([simNBNF_in.GetBinContent(i) for i in range(1,Nbins)])
        simNBNF_err = np.asarray([simNBNF_in.GetBinError(i) for i in range(1,Nbins)])
        
        simNBNF_fit = simNBNF_in.Fit('pol1','SQN','',0,30)
        fit_a = simNBNF_fit.Parameter(0)
        fit_b = simNBNF_fit.Parameter(1)

        
        ax.errorbar(simNF,simNBNF,yerr=simNBNF_err,linestyle='',\
                    marker='s',markersize=10,color='black',label='simulation',\
                    zorder=9)

        nf = np.linspace(0,30,31)
       # ax.plot(nf,self.fit(nf),linestyle='-',marker='',markersize=8,color='black',\
       #         linewidth=3,label='simNBNF self.fit',zorder=8)

############## Lines with different number of soft pomerons ####################
        NPOMList = []
        for i in range(7):
            SxH0_name = "NPOM_%02d_00"%i
            SxH0_nf_name = "NPOM_NF_%02d_00"%i
            SxH0 = ROOT.gROOT.FindObject(SxH0_name)
            SxH0_nf = ROOT.gROOT.FindObject(SxH0_nf_name)
            SxH0.Divide(SxH0_nf)
            temp_hist = [SxH0.GetBinContent(j) for j in range(1,Nbins)]
            NPOMList.append(np.asarray(temp_hist))

        nf = np.linspace(0,30,31)
        nfnew = np.linspace(nf.min(),nf.max(),300)
        graphs = []
        for i,npom in enumerate(NPOMList):
            npom = np.trim_zeros(npom,trim='b')[:11+i]
            smooth = spline(nf[:len(npom)],npom,nfnew)
            #ax.plot(nf[:len(npom)],npom,marker='',linestyle='-',label='NPOMS {}'.format(i))
            smooth = np.trim_zeros(smooth,trim='b')
            ax.plot(nfnew[:len(smooth)],smooth,linestyle='-',linewidth=4,\
                    color='grey',label='NPOMS {}'.format(i),zorder=6)


################## All soft pomerons and 0 hard pomerons########################
        self.ReOpen()
        simS = np.zeros(239)
        S = ROOT.gROOT.FindObject("NPOM_00_00")
        S_nf = ROOT.gROOT.FindObject("NPOM_NF_00_00")
        for i in range(1,25):
            S_name = "NPOM_%02d_00"%i
            S_nf_name = "NPOM_NF_%02d_00"%i
            S.Add(ROOT.gROOT.FindObject(S_name))
            S_nf.Add(ROOT.gROOT.FindObject(S_nf_name))
        S.Divide(S_nf)
        simS_err = np.asarray([S.GetBinError(i) for i in range(1,Nbins)])
        simS = np.asarray([S.GetBinContent(i) for i in range(1,Nbins)])
        simS = np.trim_zeros(simS,trim='b')[:-3]
        simSlen = len(simS)
        simS_err = simS_err[:simSlen]
        simSNF = np.linspace(0,simSlen,simSlen)
        
        simSNFnew = np.linspace(simSNF.min(),simSNF.max(),300)
        smooth = spline(simSNF,simS,simSNFnew)
        ax.plot(simSNFnew,smooth,linestyle='--',linewidth=4,color='grey',\
                label='all npoms',zorder=7)
        #ax.errorbar(simSNF,simS,yerr=simS_err,linestyle='',\
        #            marker='*',markersize=12,color='black',label='all npoms')
        

######################### Plot settings ################################################
        [ax.text(-0.4,npom[0]-0.1,str(i),fontsize=16) for i,npom in enumerate(NPOMList)]

        fontdict = {'fontsize':27,'weight':'bold'}
        ax.set_ylabel(r'$\langle n_B(n_F)\rangle$',fontdict=fontdict)
        ax.set_xlabel('$n_F$',fontdict=fontdict)

        handles, labels = ax.get_legend_handles_labels()    
        #leg = plt.legend(handles,labels)
        leg = plt.legend((handles[8],handles[9],handles[0],handles[7]),#handles[10]),\
                         (\
                          'ALICE pp, $\sqrt{s}$=7 TeV',\
                          r'QGSM $0.3<p_T<1.5$ GeV/c, $0.2<\eta<0.8$',\
                          'NPOMH = 0 and for NPOMS = \{0,...,6\}',\
                          'All NPOMS and NPOMH = 0'\
                          ),\
                          #'{:.3f}+{:.3f}x self.fit of simulated data'.format(self.fit_a,self.fit_b),\
                          loc='best',fontsize=24)

        leg.get_frame().set_alpha(0.0)
        plt.savefig('/home/roar/master/qgsm_analysis_tool/ana/analyzed/nbnf_allnpoms_0npomh.pdf')
        #plt.show()
        self.Close()

    def var_NPOMS(self):
        
        fig, ax   = self.plotSetup(xlim=[-1,25],ylim=[0,15])
        #fig1, ax1 = self.plotSetup(xlim=[-1,25],ylim=[0,15])
        outFitCSV = [['label','a','b']]
######################## Experimental data #####################################
        expNF,expNBNF,expXerr,expYerr = np.loadtxt('/home/roar/master/qgsm_analysis_tool/ana/out/nbnf_7000_exp',unpack=True)
        ax.errorbar(expNF,expNBNF,yerr=expYerr,linestyle='',\
                    marker='o',markersize=10,color='black',\
                    label='experiment',zorder=10)
######################## Experiamental self.fit #####################################
        rulerMeasuredX = [-0.5,30]
        rulerMeasuredY = [0.64,17.57]
        cof_b = (rulerMeasuredY[1]-rulerMeasuredY[0])/(rulerMeasuredX[1]-rulerMeasuredX[0])
        cof_a = rulerMeasuredY[0] + cof_b*0.5
        label = 'experiment fit'
        #outFitCSV.append([label,cof_a,cof_b])
        outFitCSV.append([label,'{:.4}'.format(cof_a),'{:.4}'.format(cof_b)])
        ax.plot(expNF,self.fit(cof_a,cof_b,expNF),linestyle='--',linewidth=2,alpha=1,label='experiment fit')

######################### Simulated for all npom ###############################
        self.ReOpen()
        simNBNF_in  = ROOT.gROOT.FindObject("ptcut_div")
        Nbins       = simNBNF_in.GetNbinsX()
        simNF       = np.linspace(0,Nbins,Nbins-1)
        simNBNF     = np.asarray([simNBNF_in.GetBinContent(i) for i in range(1,Nbins)])
        simNBNF_err = np.asarray([simNBNF_in.GetBinError(i) for i in range(1,Nbins)])
        temp_fit = simNBNF_in.Fit('pol1','SQN','',0,30)
        cof_a = temp_fit.Parameter(0)
        cof_b = temp_fit.Parameter(1)
        label='all npoms'
        #outFitCSV.append([label,cof_a,cof_b])
        outFitCSV.append([label,'{:.4}'.format(cof_a),'{:.4}'.format(cof_b)])
        ax.plot(simNF[:30],self.fit(cof_a,cof_b,simNF[:30]),linestyle='--',linewidth=2,\
                alpha=1,label=label+' fit')
        ax.errorbar(simNF,simNBNF,yerr=simNBNF_err,linestyle='',\
                    marker='s',markersize=10,color='black',label=label,\
                    zorder=9)
         
# (sum_{npoms=0}^N sum_{npomh=0}^{all} <n_{B}(n_{F})>)/(sum_{npoms=0}^N sum_{npomh=0}^{all} <n_{F}>) #
        self.ReOpen()
        NSNH = 25
        Sstart =17 
        N = [8,7,6]
        for k,n in enumerate(N):
            for i in range(n+1):
                for j in range(NSNH):
                    if i or j:
                        tempS.Add   (ROOT.gROOT.FindObject("NPOM_{:02d}_{:02d}".format(i,j)))
                        tempS_nf.Add(ROOT.gROOT.FindObject("NPOM_NF_{:02d}_{:02d}".format(i,j)))
                    else:
                        tempS      = ROOT.gROOT.FindObject("NPOM_{:02d}_{:02d}".format(i,j))
                        tempS_nf   = ROOT.gROOT.FindObject("NPOM_NF_{:02d}_{:02d}".format(i,j))

            tempS.Divide(tempS_nf)
            Nbins = tempS.GetNbinsX()
            tempAllS_nf  = np.linspace(0,Nbins,Nbins-1)
            tempAllS     = np.asarray([tempS.GetBinContent(i) for i in range(1,Nbins)])
            tempAllS_err = np.asarray([tempS.GetBinError(i) for i in range(1,Nbins)])
            label = '$N = {}$'.format(n)
            #ax.errorbar(tempAllS_nf,tempAllS,yerr=tempAllS_err,linewidth=1,label=label)
            #ax.plot(tempAllS_nf,tempAllS,label=label)

            temp_fit = tempS.Fit('pol1','SQN','',0,30)
            cof_a = temp_fit.Parameter(0)
            cof_b = temp_fit.Parameter(1)
            outFitCSV.append([label,'{:.4}'.format(cof_a),'{:.4}'.format(cof_b)])
            ax.plot(tempAllS_nf[:30],self.fit(cof_a,cof_b,tempAllS_nf[:30]),label=label)


            self.ReOpen()

######################## Plot settings #########################################
        fontdict = {'fontsize':27,'weight':'bold'}
        ax.set_ylabel(r'$\langle n_B(n_F)\rangle$',fontdict=fontdict)
        ax.set_xlabel('$n_F$',fontdict=fontdict)
        ax.set_title('$\\sum\\limits^{N}_{NPOMS=0}\\sum\\limits^{24}_{NPOMH=0} \\langle n_{B}(n_{F})\\rangle$',fontdict=fontdict)
        handles, labels = ax.get_legend_handles_labels()
        leg = plt.legend(handles,labels,loc='best')
        leg.get_frame().set_alpha(0.0)
####################### Closing remarks ########################################
        #header = 'label,a,b'
        #print(outFitCSV)
        fmt = '%s,%s,%s'
        np.savetxt('/home/roar/master/qgsm_analysis_tool/ana/analyzed/fits.csv',outFitCSV,fmt=fmt)
        plt.savefig('/home/roar/master/qgsm_analysis_tool/ana/analyzed/nbnf_Nnpoms_allnpomh.pdf')
        #plt.show()

        self.Close()

    def fix_S_var_H(self):
        self.ReOpen()
        fig, ax   = self.plotSetup(xlim=[-1,25],ylim=[0,6])
        NS = 2
        NH = [0,1,2,3,4]

        for h in NH:
            tempS      = ROOT.gROOT.FindObject("NPOM_{:02d}_{:02d}".format(NS,h))
            tempS_nf   = ROOT.gROOT.FindObject("NPOM_NF_{:02d}_{:02d}".format(NS,h))
            tempS.Divide(tempS_nf)
            Nbins = tempS.GetNbinsX()
            tempAllS_nf  = np.linspace(0,Nbins,Nbins-1)
            tempAllS     = np.asarray([tempS.GetBinContent(i) for i in range(1,Nbins)])
            tempAllS_err = np.asarray([tempS.GetBinError(i) for i in range(1,Nbins)])
            label = '$H = {}$'.format(h)
            #ax.errorbar(tempAllS_nf,tempAllS,yerr=tempAllS_err,linewidth=1,label=label)
            ax.plot(tempAllS_nf,tempAllS,label=label)

        self.ReOpen()
        #print("NPOM_{:02d}_{:02d}".format(NS,NH[0]))
        tempS = ROOT.gROOT.FindObject("NPOM_{:02d}_{:02d}".format(NS,NH[0]))
        tempS_nf = ROOT.gROOT.FindObject("NPOM_NF_{:02d}_{:02d}".format(NS,NH[0]))
        for h in range(25):
            tempS.Add   (ROOT.gROOT.FindObject("NPOM_{:02d}_{:02d}".format(NS,h)))
            tempS_nf.Add(ROOT.gROOT.FindObject("NPOM_NF_{:02d}_{:02d}".format(NS,h)))
        tempS.Divide(tempS_nf)
        Nbins = tempS.GetNbinsX()
        tempAll_nf   = np.linspace(0,Nbins,Nbins-1)
        tempAllS     = np.asarray([tempS.GetBinContent(i) for i in range(1,Nbins)])
        tempAllS_err = np.asarray([tempS.GetBinError(i) for i in range(1,Nbins)])
        label = '$H = \sum^{{{}}}_{{0}}$'.format(h)
        #ax.errorbar(tempAllS_nf,tempAllS,yerr=tempAllS_err,linewidth=1,label=label)
        ax.plot(tempAllS_nf,tempAllS,marker='o',linestyle='-',label=label)

        fontdict = {'fontsize':27,'weight':'bold'}
        ax.set_ylabel(r'$\langle n_B(n_F)\rangle$',fontdict=fontdict)
        ax.set_xlabel('$n_F$',fontdict=fontdict)
        #ax.set_title('$\\sum\\limits^{N}_{NPOMS=0}\\sum\\limits^{24}_{NPOMH=0} \\langle n_{B}(n_{F})\\rangle$',fontdict=fontdict)
        handles, labels = ax.get_legend_handles_labels()
        leg = plt.legend(handles,labels,loc='best')
        leg.get_frame().set_alpha(0.0)
    
        plt.savefig('/home/roar/master/qgsm_analysis_tool/ana/analyzed/nbnf_fixed_s_var_h.pdf')
        #plt.show()



if __name__=="__main__":
    path ="/home/roar/master/qgsm_analysis_tool/ana/out/7TeV_4M.root" 
    P = Plot(root_file_path=path)
    #P.NBNFPicture()
    P.var_NPOMS()
    P.fix_S_var_H()
    P.Show()
