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

FILEPATH    ='/home/roar/master/qgsm_analysis_tool/ana/build/' 
F_NAME      = '13000_4M.root'

class histogram:
    def __init__(self,f,dia):

        self.th1f   = f.FindObjectAny(dia[0]) 
        self.limit  = [self.th1f.GetXaxis().GetXmin(),self.th1f.GetXaxis().GetXmax()]
        self.nb     = self.th1f.GetNbinsX()
        self.f      = f 
        self.dia    = dia

        self.adding()

    def adding(self):
        for d in self.dia[1:]:
            print(d)
            self.th1f.Add(self.f.FindObjectAny(d))
            
    def draw(self,label=''):
        NF      = np.asarray([self.th1f.GetBinContent(i) for i in range(1,self.nb)])
        #NF      = np.trim_zeros(NF,trim='b')
        NFx     = np.linspace(self.limit[0],self.limit[1],len(NF))
        fig,ax  = plt.subplots()

        ax.plot(NFx,NF,linestyle='',marker='o',label=label)
        ax.legend()

    def close(self):
        self.f.Close()

            
if __name__=='__main__':

    f           = ROOT.TFile(FILEPATH+F_NAME)
    dou_dia     = ['DOU_NF_01']#,'DOU_NF_02','DOU_NF_03']
    dou         = histogram(f,dou_dia)
    dou.draw(label='dou')
    dou.close()

    #f           = ROOT.TFile(FILEPATH+F_NAME)
    #nf_alldia   = ['DOU_NF_01','DOU_NF_02','DOU_NF_03','SIN_NF_01','SIN_NF_02','SIN_NF_03']
    #nf_all      = histogram(f,nf_alldia)
    #nf_all.draw(label='all')
    #nf_all.close()
    
    #f           = ROOT.TFile(FILEPATH+F_NAME)
    #sin_dia     = ['SIN_NF_01','SIN_NF_02','SIN_NF_03']
    #sin         = histogram(f,sin_dia)
    #sin.draw(label='sin')
    #sin.close()

    #f               = ROOT.TFile(FILEPATH+F_NAME)
    #nbnf_alldia     = ['DOU_NBNF_01','DOU_NBNF_02','DOU_NBNF_03',\
    #                   'SIN_NBNF_01','SIN_NBNF_02','SIN_NBNF_03']
    #nbnf_all        = histogram(f,nbnf_alldia)
    #nbnf_all.draw(label='nbnfall')
    #nbnf_all.close()

    #f               = ROOT.TFile(FILEPATH+F_NAME)
    #nbnf_alldia     = ['DOU_NBNF_01','DOU_NBNF_02','DOU_NBNF_03']
    #nbnf_all        = histogram(f,nbnf_alldia)
    #nbnf_all.draw(label='nbnfdou')
    #nbnf_all.close()

    #f               = ROOT.TFile(FILEPATH+F_NAME)
    #nbnf_alldia     = ['SIN_NBNF_01','SIN_NBNF_02','SIN_NBNF_03']
    #nbnf_all        = histogram(f,nbnf_alldia)
    #nbnf_all.draw(label='nbnfsin')
    #nbnf_all.close()
    
    plt.show()
