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

class histogram:
    def __init__(self,f,th1f,limit,nb,dia):
        self.limit  = limit
        self.nb     = nb
        self.th1f   = th1f 
        self.f      = f 
        self.dia    = dia

        self.adding()

    def adding(self):
        for d in self.dia:
            self.th1f.Add(self.f.FindObjectAny(d))
            
    def draw(self,label=''):
        NF      = np.asarray([self.th1f.GetBinContent(i) for i in range(1,self.nb)])
        NF      = np.trim_zeros(NF,trim='b')
        NFx     = np.linspace(self.limit[0],self.limit[1],len(NF))
        fig,ax  = plt.subplots()

        ax.plot(NFx,NF,linestyle='',marker='o',label=label)

            
if __name__=='__main__':
    filepath    = '/home/roar/master/qgsm_analysis_tool/ana/build/'
    f           = ROOT.TFile(filepath+'13000_4M.root')
    nf_limit    = [-10,10]
    nf_nb       = 81
    dou_dia     = ['DOU_NF_01','DOU_NF_02','DOU_NF_03']
    dou_th1f    = ROOT.TH1F("dou","DOU",nf_nb,nf_limit[0],nf_limit[1])
    dou         = histogram(f,dou_th1f,nf_limit,nf_nb,dou_dia)
    dou.draw(label='dou')
    #plt.show()

    f.Close()
    f           = ROOT.TFile(filepath+'13000_4M.root')
    nf_alldia   = ['DOU_NF_01','DOU_NF_02','DOU_NF_03','SIN_NF_01','SIN_NF_02','SIN_NF_03']
    nf_allth1f  = ROOT.TH1F("all","ALL",nf_nb,nf_limit[0],nf_limit[1])
    nf_all      = histogram(f,nf_allth1f,nf_limit,nf_nb,nf_alldia)
    nf_all.draw(label='all')
    
    f.Close()
    f           = ROOT.TFile(filepath+'13000_4M.root')
    
    sin_dia     = ['SIN_NF_01','SIN_NF_02','SIN_NF_03']
    sin_th1f    = ROOT.TH1F("sin","SIN",nf_nb,nf_limit[0],nf_limit[1])
    sin         = histogram(f,sin_th1f,nf_limit,nf_nb,sin_dia)
    sin.draw(label='sin')

    f.Close()
    f           = ROOT.TFile(filepath+'13000_4M.root')

    nbnf_limit      = [-0.5,600.5]
    nbnf_nb         = 600
    nbnf_alldia     = ['DOU_NBNF_01','DOU_NBNF_02','DOU_NBNF_03',\
                       'SIN_NBNF_01','SIN_NBNF_02','SIN_NBNF_03']
    nbnf_allth1f    = ROOT.TH1F("nbnfall","NBNFALL",nbnf_nb,nbnf_limit[0],nbnf_limit[1])
    nbnf_all        = histogram(f,nbnf_allth1f,nbnf_limit,nbnf_nb,nbnf_alldia)
    nbnf_all.draw(label='nbnfall')
    
    plt.show()
