#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 14:09:16 2022

@author: dj
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

def plot_band(ax,s_id):
    client=pymongo.MongoClient("mongodb://192.168.50.14:27017")
    db=client['layers']
    col=db['stacking_ionic']
    
    data=col.find_one({'structure_id':s_id})
    band=np.array(data['band_opp']['eigval'])
    efermi=data['band_opp']['efermi']
    high_sym=np.array(data['band_opp']['high_sym'])
    ticks=data['band_opp']['ticks']
    if data['band_opp']['is_metal'] == False:
        efermi = data['band_opp']['vbm']['energy']
    band_number=int(len(band[0,:])-1)
      
    xlim=[0,high_sym[-1]]
    ylim=[-2,2]
    ax.axis([xlim[0],xlim[1],ylim[0],ylim[1]])
    band_range=[]
    for ii in np.arange(1,band_number+1):
        ax.plot(band[:,0],band[:,ii]-efermi,'-',color='b',linewidth=0.5)
        band_range.append([np.min(band[:,ii]-efermi),np.max(band[:,ii]-efermi)])
    
    # ax.plot(band[:,0],band[:,38]-efermi,'-',color='r',linewidth=0.5)
    # ax.scatter(band[90,0],band[90,38]-efermi,5,color='k',linewidth=0.5)
    
    font={'family':'Arial','weight':'normal','size':6}    
    
    ax.set_ylabel('Energy (eV)',{'family':'Arial','weight':'normal','size':8})
    
    #plot high symmetry lines
    for i in np.arange(1,len(high_sym)):
        ax.plot([high_sym[i],high_sym[i]],[-8,8],color='grey',linewidth=0.5,zorder=0)
    
    #draw Fermi Energy
    x=np.array([0,high_sym[-1]])
    y=np.array([0,0])
    ax.plot(x,y,'k--',linewidth=0.5,zorder=0)
    
    #set high symmetry label
    font={'family':'Arial','weight':'normal','size':6}   
    label=ticks
    for i in range(len(label)):
        if label[i][0]== "\\":
            label[i]="$"+label[i]+"$"
    ax.set_xticks([])
    d=-0.04
    high_sym=high_sym/high_sym[-1] 
    for i in np.arange(1,len(ticks)-1):
        ax.text(high_sym[i]-0.01,d,label[i],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
    ax.text(high_sym[0]-0.005,d,label[0],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
    ax.text(high_sym[-1]-0.01,d,label[-1],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
    
    #set yticks
    ax.minorticks_on()
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
    ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)
    ax.tick_params(labelsize=6)
    labels=ax.get_xticklabels()+ax.get_yticklabels()
    [label.set_fontname('Arial') for label in labels]
    
    #set line width
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.set_yticks(np.arange(ylim[0],ylim[1]+0.1,0.5))

def plot_band_spin(ax,s_id):
    client=pymongo.MongoClient("mongodb://192.168.50.14:27017")
    db=client['layers']
    col=db['stacking_ionic']
    
    data=col.find_one({'structure_id':s_id})
    band_up=np.array(data['band_opp']['eigval_up'])
    band_dw=np.array(data['band_opp']['eigval_dw'])
    efermi=data['band_opp']['efermi']
    high_sym=np.array(data['band_opp']['high_sym'])
    ticks=data['band_opp']['ticks']
    if data['band_opp']['is_metal'] == False:
        efermi = data['band_opp']['vbm']['energy']
    band_number=int(len(band_up[0,:])-1)
      
    xlim=[0,high_sym[-1]]
    ylim=[-2,2]
    ax.axis([xlim[0],xlim[1],ylim[0],ylim[1]])
    band_range=[]
    for ii in np.arange(1,band_number+1):
        ax.plot(band_up[:,0],band_up[:,ii]-efermi,'-',color='b',linewidth=0.5)
    for ii in np.arange(1,band_number+1):
        ax.plot(band_dw[:,0],band_dw[:,ii]-efermi,'-',color='r',linewidth=0.5)

    
    # ax.plot(band[:,0],band[:,38]-efermi,'-',color='r',linewidth=0.5)
    # ax.scatter(band[90,0],band[90,38]-efermi,5,color='k',linewidth=0.5)
    
    font={'family':'Arial','weight':'normal','size':6}    
    
    ax.set_ylabel('Energy (eV)',{'family':'Arial','weight':'normal','size':8})
    
    #plot high symmetry lines
    for i in np.arange(1,len(high_sym)):
        ax.plot([high_sym[i],high_sym[i]],[-8,8],color='grey',linewidth=0.5,zorder=0)
    
    #draw Fermi Energy
    x=np.array([0,high_sym[-1]])
    y=np.array([0,0])
    ax.plot(x,y,'k--',linewidth=0.5,zorder=0)
    
    #set high symmetry label
    font={'family':'Arial','weight':'normal','size':6}   
    label=ticks
    for i in range(len(label)):
        if label[i][0]== "\\":
            label[i]="$"+label[i]+"$"
    ax.set_xticks([])
    d=-0.04
    high_sym=high_sym/high_sym[-1] 
    for i in np.arange(1,len(ticks)-1):
        ax.text(high_sym[i]-0.01,d,label[i],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
    ax.text(high_sym[0]-0.005,d,label[0],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
    ax.text(high_sym[-1]-0.01,d,label[-1],{'family':'Arial','weight':'normal','size':6},transform=ax.transAxes)
    
    #set yticks
    ax.minorticks_on()
    ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
    ax.tick_params(axis='y',direction='in',which="major",width=0.5,length=3)
    ax.tick_params(axis='y',direction='in',which="minor",width=0.5,length=1.5)
    ax.tick_params(labelsize=6)
    labels=ax.get_xticklabels()+ax.get_yticklabels()
    [label.set_fontname('Arial') for label in labels]
    
    #set line width
    ax.spines['bottom'].set_linewidth(0.5)
    ax.spines['left'].set_linewidth(0.5)
    ax.spines['right'].set_linewidth(0.5)
    ax.spines['top'].set_linewidth(0.5)
    ax.set_yticks(np.arange(ylim[0],ylim[1]+0.1,0.5))    
    

if __name__ == '__main__':
    client=pymongo.MongoClient("mongodb://192.168.50.14:27017")
    db=client['layers']
    col=db['stacking_ionic']
    query={'flag':'keep','$or':[{'phonon.stability':'high'},{'phonon.stability':'stable'}]}
    data=col.find(query,{'structure_id':1,'band_opp.is_spin_polarized':1})
    for num,i in enumerate(data):
        s_id = i['structure_id']
        fig=plt.figure()
        ax=fig.add_subplot(111)
        if i['band_opp']['is_spin_polarized'] == True:
            plot_band_spin(ax,s_id)
        else:
            plot_band(ax,s_id)
        
        # fig=plt.gcf()
        cm=10
        inch=0.3937007874
        fig.set_size_inches(cm*inch,9*inch)
        plt.savefig(r'pdf\%s.pdf'%s_id,dpi=600,bbox_inches='tight')
        plt.close()
        print(num,s_id)