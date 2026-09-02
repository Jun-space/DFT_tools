# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 09:31:06 2022

@author: DengJun
"""


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pymongo
import matplotlib as mpl
mpl.use('Agg')

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'

def plot_phonon(ax,s_id):
    client=pymongo.MongoClient("mongodb://192.168.50.14:27017")
    db=client['layers']
    col=db['stacking_ionic']
    
    
    data=col.find_one({'structure_id':s_id})
    
    distances=data['phonon']['data']['band']['distances']
    for num,i in enumerate(distances):
        distances[num]=np.array(i)
    
    frequencies=data['phonon']['data']['band']['frequencies']
    for num,i in enumerate(frequencies):
        frequencies[num]=np.array(i)
    
    # labels=data['phonon']['data']['labels']
    # tmp=[]
    # for num,i in enumerate(labels[:-1]):
    #     tmp.append(labels[num][-1]+'|'+labels[num+1][0])
    # ticks=labels[0][:-1]+tmp+[labels[-1][-1]]
    ticks=data['band_opp']['ticks']
    high_sym=[]
    for i in distances:
        high_sym.append(i[0])
    high_sym.append(distances[-1][-1])
    
    freq_max=[]
    freq_min=[]
    for i in frequencies:
        freq_max.append(np.max(i))
        freq_min.append(np.min(i))
    freq_max=max(freq_max)
    freq_min=min(freq_min)
    
    ymax=np.ceil(freq_max)
    ymin=freq_min-0.2
    
    xlim=[0,high_sym[-1]]
    ylim=[ymin,ymax]
    ax.axis([xlim[0],xlim[1],ylim[0],ylim[1]])
    for i in range(len(distances)):
        for j in range(len(frequencies[0][0,:])):
            ax.plot(distances[i],frequencies[i][:,j],'-',color='b',linewidth=0.5)
            
    ax.set_ylabel('Frequency (THz)',{'family':'Arial','weight':'normal','size':8})
    for i in np.arange(1,len(high_sym)):
        ax.plot([high_sym[i],high_sym[i]],[ymin,ymax],color='grey',linewidth=0.5,zorder=0)
    
    #draw Fermi Energy
    x=np.array([0,high_sym[-1]])
    y=np.array([0,0])
    ax.plot(x,y,'k--',linewidth=0.5,zorder=0)
    
    font={'family':'Arial','weight':'normal','size':8}   
    label=ticks
    for i in range(len(label)):
        if label[i][0]== "\\":
            label[i]="$"+label[i]+"$"
    ax.set_xticks([])
    d=-0.04
    high_sym=high_sym/high_sym[-1] 
    for i in np.arange(1,len(ticks)-1):
        ax.text(high_sym[i]-0.01,d,label[i],{'family':'Arial','weight':'normal','size':7},transform=ax.transAxes)
    ax.text(high_sym[0]-0.005,d,label[0],{'family':'Arial','weight':'normal','size':7},transform=ax.transAxes)
    ax.text(high_sym[-1]-0.01,d,label[-1],{'family':'Arial','weight':'normal','size':7},transform=ax.transAxes)
    
    
    #set yticks
    ax.minorticks_on()
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
    ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)
    ax.tick_params(labelsize=8)
    labels=ax.get_xticklabels()+ax.get_yticklabels()
    [label.set_fontname('Arial') for label in labels]
    
    #set line width
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    # ax.set_yticks(np.arange(0,ylim[1],0.5))

if __name__ == '__main__':
    client=pymongo.MongoClient("mongodb://192.168.50.14:27017")
    db=client['layers']
    col=db['stacking_ionic']
    query={'flag':'keep','$or':[{'phonon.stability':'high'},{'phonon.stability':'stable'}]}
    data=col.find(query,{'structure_id':1})
    for num,i in enumerate(data):
        s_id = i['structure_id']
        fig=plt.figure()
        ax=fig.add_subplot(111)
        plot_phonon(ax,s_id)
        
        # fig=plt.gcf()
        cm=10
        inch=0.3937007874
        fig.set_size_inches(cm*inch,10*inch)
        plt.savefig(r'pdf\%s.pdf'%s_id,dpi=600,bbox_inches='tight')
        plt.close()
        print(num,s_id)