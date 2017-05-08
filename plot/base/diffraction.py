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

FILEPATH    ='/home/roar/master/qgsm_analysis_tool/ana/build/' 
F_NAME      = '900_4M.root'
DIA      = {'SIN0':'1','SIN1':'6','SIN2':'10','DOU0':'11','DOU1':'21','DOU2':'31'}
ETALIM      = {'0':'eta < 4','1':'eta < 6','2':'eta < 8','3':'eta < 10'}

class histogram:
    def __init__(self,f,lim=0,NBNF='NF',SINDIA=[],DOUDIA=[]):
        SIN = 'SIN_{0}_{1}{2}' 
        DOU = 'DOU_{0}_{1}{2}'
        if 'NBNF'== NBNF:
            dianbnf = SIN.format(NBNF,lim,SINDIA[0]) if SINDIA else DOU.format(NBNF,lim,DOUDIA[0]) 
            dia     = SIN.format('NBNFNF',lim,SINDIA[0]) if SINDIA else DOU.format('NBNFNF',lim,DOUDIA[0]) 
            self.th1f       = f.FindObjectAny(dianbnf) 
            self.th1fnf     = f.FindObjectAny(dia) 
            self.nb         = self.th1f.GetNbinsX()
            self.limit      = [self.th1f.GetXaxis().GetXmin(),self.th1f.GetXaxis().GetXmax()]
            self.dianbnf    = [SIN.format('NBNF',lim,l) for l in SINDIA]+\
                              [DOU.format('NBNF',lim,l) for l in DOUDIA]
            self.dia        = [SIN.format('NBNFNF',lim,l) for l in SINDIA]+\
                              [DOU.format('NBNFNF',lim,l) for l in DOUDIA]
        else:
            dia     = SIN.format('NF',lim,SINDIA[0]) if SINDIA else DOU.format('NF',lim,DOUDIA[0]) 
            self.th1f   = f.FindObjectAny(dia) 
            self.nb     = self.th1f.GetNbinsX()
            self.limit  = [self.th1f.GetXaxis().GetXmin(),self.th1f.GetXaxis().GetXmax()]
            self.dia    = [SIN.format('NF',lim,l) for l in SINDIA]+\
                          [DOU.format('NF',lim,l) for l in DOUDIA]

        self.f      = f 
        self.NBNF   = NBNF
        self.adding()

    def adding(self):
        if 'NBNF' == self.NBNF:
            for nbnf,nf in zip(self.dianbnf[1:],self.dia[1:]):
                print(nbnf)
                self.th1f  .Add(self.f.FindObjectAny(nbnf))
                self.th1fnf.Add(self.f.FindObjectAny(nf))
            self.th1f.Divide(self.th1fnf)

        else:
            for d in self.dia[1:]:
                print(d)
                self.th1f.Add(self.f.FindObjectAny(d))
            
    def draw(self):
        labels          = ['{}'.format(DIA[l[:3]+l[-1]]) for l in self.dia]
        label           = '{} {}'.format('Diagrams ',', '.join(labels))
        title           = '{} with $\eta <$ {}'\
                .format('$<n_B(n_F)>$'if self.NBNF else '$\eta$',ETALIM[self.dia[0][-2]][-2:])
        print(title)
        NF              = np.asarray([self.th1f.GetBinContent(i) for i in range(1,self.nb+1)])
        NF              = NF[:35] if self.NBNF == 'NBNF' else NF
        self.limit[1]   = len(NF)-1 if self.NBNF=='NBNF' else self.limit[1] 
        NFx             = np.linspace(self.limit[0],self.limit[1],len(NF))

        fig,ax  = plt.subplots()

        ax.plot(NFx,NF,linestyle='-',marker='o',label=label)
        ax.set_title(title)
        ax.set_xlabel('$n_F$')
        ax.set_ylabel('$<n_B(n_F)>$') if self.NBNF=='NBNF' else ax.set_ylabel('$\eta$')
        ax.grid('on')
        ax.legend()
        filename =\
        'temp_plots/{}_eta{}_dia{}.pdf'.\
        format(self.NBNF,ETALIM[self.dia[0][-2]][-2:],DIA[self.dia[0][:3]+self.dia[0][-1]])\
                .replace(" ","")
        fig.savefig(filename)

    def close_file(self):
        self.f.Close()

            
if __name__=='__main__':

    SINL = [[0,1,2]]
    DOUL = [[0,1,2]]
    limits = [0,1,2,3]
    for NBNF in ['NBNF']:#,'NF']:
        for i in limits:
            for j,k in zip(SINL,DOUL):
                f           = ROOT.TFile(FILEPATH+F_NAME)
                lim         = i 
                SIN         = j#[0,1,2]
                DOU         = k#[0,1,2]
                #NBNF        = 'NBNF'
                hist        = histogram(f,lim,NBNF,SIN,DOU)
                hist.draw()
                hist.close_file()

    #f           = ROOT.TFile(FILEPATH+F_NAME)
    #lim         = 0 
    #SIN         = [0]#[0,1,2]
    #DOU         = []#[0,1,2]
    #NBNF        = 'NF'
    #hist        = histogram(f,lim,NBNF,SIN,DOU)
    #hist.draw()
    #hist.close_file()
    #plt.show()
