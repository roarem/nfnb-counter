import matplotlib as mpl
#mpl.use('Agg')
mpl.rc('text',usetex=True)
mpl.rcParams['legend.numpoints']=1
mpl.rcParams['font.size'] = 27
mpl.rcParams['font.weight']   = 'bold'
mpl.rcParams['text.latex.preamble']=[r'\usepackage{bm} \boldmath']
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def bcorr(filepath):

    exps = []; experrs = []

    x = np.asarray([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 0.0, 0.4, 0.8, 0.0, 0.4, 0.0])
    exps.append(np.asarray([0.212, 0.203, 0.193, 0.182, 0.172, 0.163, 0.159, 0.335, 0.3, 0.274,\
                         0.406, 0.368, 0.452]))
    experrs.append(np.asarray([0.008935882720806042, 0.007034912934784624, 0.00795110055275369,\
                            0.00751065909225016, 0.007930952023559342, 0.007516648189186454,\
                            0.007256031973468695, 0.008836288813749808, 0.008945389874119518,\
                            0.008628441342444185, 0.009135097153287424, 0.008845903006477066,\
                            0.009334345183246651]))

    exps.append(np.asarray([0.302, 0.294, 0.285, 0.269, 0.259, 0.253, 0.247, 0.447, 0.413, 0.386,\
                          0.525, 0.488, 0.572]))
    experrs.append(np.asarray([0.011016351483136328, 0.011029052543169788, 0.008009993757800314,\
                             0.007011419257183242, 0.010007996802557442, 0.011022250223978767,\
                             0.011004090148667448, 0.014012851244482689, 0.015005332385522154,\
                             0.01700264685276972, 0.017007351351694948, 0.021002142747824568,\
                             0.022005681084665385]))

    exps.append(np.asarray([0.366, 0.358, 0.345, 0.334, 0.327, 0.316, 0.311, 0.521,\
                         0.487, 0.463, 0.598, 0.564, 0.643]))
    experrs.append(np.asarray([0.00852877482408816, 0.007910120100226039, 0.00840535543567314,\
                            0.008514693182963202, 0.007134423592694788, 0.006007495318350236,\
                            0.00920869154657707, 0.008728115489611717, 0.011901680553602505,\
                            0.01271927670899568, 0.010117806086301516, 0.012403628501370072,\
                            0.010117806086301516]))

    #plots = np.asarray([bcorrPlotSetup() for i in range(4)])
    #figs = plots[:,0]
    #axs  = plots[:,1]
    fig,axs = bcorrPlotSetup()
    csvfiles = ['900_4M_bcorr.csv','2760_4M_bcorr.csv','7000_4M_bcorr.csv','13000_4M_bcorr.csv'] 
    simbcorr = [np.loadtxt(filepath+csvfile,skiprows=1) for csvfile in csvfiles]
    
    delta = 0; fontsize=27; markersize=10
    fontdict = {'fontsize':27,'weight':'bold'}
    #print(simbcorr[2]-[0.36672, 0.36456, 0.35950, 0.35500, 0.34924, 0.34422, 0.34115, 0.53123, 0.51858, 0.50338, 0.63815,0.62057, 0.69480])

    for k,ax in enumerate(axs):
        delta = 0
        for i,j in zip([0,7,10,12],[7,10,12,13]):
            if k<3:
                ax.errorbar(x[i:j],exps[k][i:j],experrs[k][i:j],marker='s',markersize=markersize,\
                        linestyle='',color='grey',label='ALICE')
            ax.text(-0.13,simbcorr[k][i],'0.{:d}'.format(delta),fontsize=fontsize)

            ax.plot(x[i:j],simbcorr[k][i:j],marker='o',markersize=markersize,\
                    linestyle='--',color='black',label='QGSM')

            delta +=2

    for i,ax in enumerate(axs):
        ax.text(-0.13,simbcorr[i][12]+0.09,'$\delta\eta$',fontsize=fontsize)
        handles, labels = ax.get_legend_handles_labels()
        if i<2:
            handles,labels=(handles[0],handles[-1]),(labels[0],labels[1])
        else:
            handles,labels=(handles[0],),(labels[0],)
        leg = ax.legend((handles[0],),(labels[0],),loc='best')
        leg.get_frame().set_alpha(0.0)
    #plt.show()
    plt.savefig('testy.pdf')

def bcorrPlotSetup():

    titles = ['$900\, GeV$','$2760\, GeV$','$7000\, GeV$','$13000\, GeV$']
    
    fig, axs = plt.subplots(2,2)
    axs = axs.reshape(4)
    fig.subplots_adjust(wspace=0.001,hspace=0.001)
    DPI = fig.get_dpi()
    size = 1000
    fig.set_size_inches(1500/DPI,1500/DPI)

    majorLocator = ticker.MultipleLocator(0.3)
    minorLocator = ticker.MultipleLocator(0.1)
    majorFormatter = ticker.FormatStrFormatter('%.1f')
    minorFormatter = ticker.FormatStrFormatter('%.1f')
    for ax in axs:
        ax.set_xlim(-0.2,1.3)
        ax.set_ylim(0.1,0.9)

        #x0,x1 = ax.get_xlim()
        #y0,y1 = ax.get_ylim()
        #ax.set_aspect((x1-x0)/(y1-y0))

        ax.grid(which='minor',alpha=1)

        ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.1))
        ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax.xaxis.set_major_locator(ticker.MultipleLocator(0.2))
        ax.xaxis.set_major_formatter(majorFormatter)

        [tick.label.set_fontsize(20) for tick in ax.xaxis.get_major_ticks()]
        [tick.label.set_fontsize(20) for tick in ax.yaxis.get_major_ticks()]
        #[tick.label.set_fontsize(20) for tick in ax.xaxis.get_minor_ticks()]
        #[tick.label.set_fontsize(20) for tick in ax.yaxis.get_minor_ticks()]

        ax.xaxis.set_tick_params(which='major',length=14,width=3)
        ax.yaxis.set_tick_params(which='major',length=14,width=3)
        ax.xaxis.set_tick_params(which='minor',length=8 ,width=2)
        ax.yaxis.set_tick_params(which='minor',length=8 ,width=2)

        ax.set_xlabel('$\eta$')
        ax.set_ylabel('$b_{corr}$')

    
    [label.set_visible(False) for label in axs[0].get_xticklabels()]
    [label.set_visible(False) for label in axs[1].get_xticklabels()]
    [label.set_visible(False) for label in axs[1].get_yticklabels()]
    [label.set_visible(False) for label in axs[3].get_yticklabels()]
    [axs[2].get_yticklabels()[-i].set_visible(False) for i in range(2)]
    [axs[2].get_xticklabels()[-i].set_visible(False) for i in range(3)]
    [ax.text(0.5,0.9,tit,horizontalalignment='center',
             verticalalignment='center',transform=ax.transAxes) for ax,tit in zip(axs,titles)]

    return fig, axs

if __name__=='__main__':
    path = "/home/roar/master/qgsm_analysis_tool/ana/out/1912/"
    bcorr(path)
