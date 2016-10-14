import matplotlib.pyplot as plt
import numpy as np
import ROOT


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
        NPOMList = []
        for i in range(7):
            SxH0_name = "NPOM_%02d_00"%i
            SxH0_nf_name = "NPOM_NF_%02d_00"%i
            SxH0 = ROOT.gROOT.FindObject(SxH0_name)
            SxH0_nf = ROOT.gROOT.FindObject(SxH0_nf_name)
            SxH0.Divide(SxH0_nf)
            Nbins = SxH0.GetNbinsX()
            temp_hist = []
            for i in range(1,Nbins):
                temp_hist.append(SxH0.GetBinContent(i))
            NPOMList.append(np.asarray(temp_hist))
            
        fig, ax = plt.subplots()
        nf = np.linspace(0,30,31)

        for npom in NPOMList:
            npom = np.trim_zeros(npom,trim='b')
            ax.plot(nf[:len(npom)],npom,marker='o',linestyle='')
        ax.set_xlim(-1,30)
        plt.show()


if __name__=="__main__":
    P = Plot(root_file_path="../out/test.root")
    P.NPOM()
