from base import General,bcorr,nbnf_corr

def Base():
    path = "/home/roar/master/qgsm_analysis_tool/ana/"
    def GeV7000(save=0,nsd=0):
        name = 'out/1612/7000_4M.root' 
        P = General.General(root_file_path=path,filename=name,save=save,nsd=nsd)
        return P
    def GeV2760(save=0,nsd=0):
        name = 'out/1612/2760_4M.root' 
        P = General.General(root_file_path=path,filename=name,save=save,nsd=nsd)
        return P
    def GeV900(save=0,nsd=0):
        name = 'out/1612/900_4M.root' 
        P = General.General(root_file_path=path,filename=name,save=save,nsd=nsd)
        return P
    def GeV13000(save=0,nsd=0):
        name = 'out/1612/13000GeV_1M.root' 
        P = General.General(root_file_path=path,filename=name,save=save,nsd=nsd)
        return P

    options = {0: GeV7000,1:GeV2760,2:GeV900,3:GeV13000}
    P = options[0](save=0,nsd=0)

    plots = {0:P.NBNFPicture,1:P.var_NPOMS,2:P.fix_S_var_H,3:P.nch_dist,4:P.nf_nb}
    plots[0]()

    P.Show()

def Bcorr():
    path = "/home/roar/master/qgsm_analysis_tool/ana/out/1912/"
    bcorr.bcorr(path)

def Nbnf_corr():
    path = '/home/roar/master/qgsm_analysis_tool/ana/out/1912/'
    nbnf_corr.nbnf(path,kin=False)

if __name__=='__main__':
    Nbnf_corr()
    #Base()
    #Bcorr()
