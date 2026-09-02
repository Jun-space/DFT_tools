# -*- coding: utf-8 -*-
"""
Created on Wed Aug 26 09:47:04 2026

@author: dengjun
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import rcParams
from matplotlib.patches import Polygon
from matplotlib.colors import LinearSegmentedColormap
import re


rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Arial'
rcParams['mathtext.it'] = 'Arial:italic'
rcParams['mathtext.bf'] = 'Arial:bold'
rcParams['pdf.fonttype'] = 42
rcParams['ps.fonttype'] = 42
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = 8

fig = plt.figure()
ax = fig.add_subplot(111)

font = {'family':'Arial','weight':'normal','size':8}

#################### readdata ######################
path='.'
band = np.loadtxt('%s/band.dat'%path)
efermi = np.loadtxt('%s/Fermi.dat'%path)
high_sym = np.loadtxt('%s/high_sym.dat'%path)

#read projection  Se  Cl  Rh
A_proj = np.loadtxt('%s/Nb.dat'%path)  # shape (Nk, Nb+1)
B_proj = np.loadtxt('%s/O.dat'%path)
C_proj = np.loadtxt('%s/I.dat'%path)

ticks = [line.strip() for line in open('%s/ticks.dat'%path, 'r').readlines()]
band_number = band.shape[1] - 1
num_kpoints = band.shape[0]

# A = A_proj[:, 1:]
# B = B_proj[:, 1:]
# C = C_proj[:, 1:]

total = A_proj[:,11] + B_proj[:,11] + C_proj[:,11] + 1e-10
A_frac = A_proj[:,11] / total
B_frac = B_proj[:,11] / total
C_frac = C_proj[:,11] / total

A_frac=A_frac.reshape(len(band[0,:])-1,len(band[:,0])).T
B_frac=B_frac.reshape(len(band[0,:])-1,len(band[:,0])).T
C_frac=C_frac.reshape(len(band[0,:])-1,len(band[:,0])).T


# high_sym=high_sym[:-3]
# ticks=ticks[:-2]
# ticks[-1]='A'

xlim = [0, high_sym[-1]]
ylim = [-3, 3]
ax.axis([xlim[0], xlim[1], ylim[0], ylim[1]])

# 开始画投影线
for band_idx in range(band_number):
    x = band[:, 0]
    y = band[:, band_idx + 1] - efermi

    # 计算线段颜色，取相邻两点的颜色平均值
    rgb = np.stack([A_frac[:, band_idx], B_frac[:, band_idx], C_frac[:, band_idx]], axis=1)  # shape (Nk, 3)
    rgb = np.sqrt(rgb)
       
    rgb_seg = 0.5 * (rgb[:-1] + rgb[1:])  # shape (Nk-1, 3)

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    lc = LineCollection(segments, colors=rgb_seg, linewidth=1)
    ax.add_collection(lc)

# 画费米能级线
ax.plot([0, high_sym[-1]], [0, 0], 'k--', linewidth=0.5, zorder=0)

# 画高对称点竖线
for hs in high_sym[1:]:
    ax.plot([hs, hs], [-8, 8], color='grey', linewidth=0.5, zorder=2)

# 设置高对称点标签
high_sym_norm = high_sym / high_sym[-1]
d = -0.05
for i, label in enumerate(ticks):
    if label.startswith("\\"):
        label = "$" + label + "$"
    if '_' in label:
        label=re.sub(r'_(.)', r'$_\1$', label)
    ax.text(high_sym_norm[i] - 0.01, d, label, fontdict={'family':'Arial','weight':'normal','size':8}, transform=ax.transAxes)

# 轴标签和格式
ax.set_ylabel('Energy (eV)', fontdict=font)
ax.set_xticks([])
ax.minorticks_on()
ax.yaxis.set_minor_locator(plt.MultipleLocator(0.5))
ax.tick_params(axis='y', direction='in', which="major", width=0.5, length=3)
ax.tick_params(axis='y', direction='in', which="minor", width=0.5, length=1.5)
ax.tick_params(labelsize=8)

for label in ax.get_xticklabels() + ax.get_yticklabels():
    label.set_fontname('Arial')

for spine in ax.spines.values():
    spine.set_linewidth(0.5)

ax.set_yticks(np.arange(ylim[0], ylim[1] + 0.1, 1))

def rgb_triangle_legend(ax):
    # 三角形顶点，顺时针或逆时针
    verts = np.array([
        [0.5, np.sqrt(3)/2],  # 顶点A（红色）
        [0, 0],               # 顶点B（绿色）
        [1, 0]                # 顶点C（蓝色）
    ])

    # 生成三角形内网格点，用于颜色插值
    res = 200
    x = np.linspace(0, 1, res)
    y = np.linspace(0, np.sqrt(3)/2, res)
    xv, yv = np.meshgrid(x, y)
    
    # 判断点是否在三角形内（利用重心坐标）
    def in_triangle(px, py):
        A = verts[0]
        B = verts[1]
        C = verts[2]
        
        v0 = C - A  # shape (2,)
        v1 = B - A
        
        # 对所有点计算 v2
        v2 = np.stack([px - A[0], py - A[1]], axis=-1)  # shape (200,200,2)
        
        dot00 = np.dot(v0, v0)
        dot01 = np.dot(v0, v1)
        dot11 = np.dot(v1, v1)
        
        dot02 = v2[...,0]*v0[0] + v2[...,1]*v0[1]  # shape (200,200)
        dot12 = v2[...,0]*v1[0] + v2[...,1]*v1[1]  # shape (200,200)
        
        invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
        u = (dot11 * dot02 - dot01 * dot12) * invDenom
        v = (dot00 * dot12 - dot01 * dot02) * invDenom
        
        return (u >= 0) & (v >= 0) & (u + v <= 1)


    mask = in_triangle(xv, yv)

    # 计算重心坐标比例(r,g,b)
    r = np.zeros_like(xv)
    g = np.zeros_like(xv)
    b = np.zeros_like(xv)

    # 红色顶点A重心坐标是u,v,w中的w，对应顶点A的权重
    # 简化：根据顶点位置计算颜色权重(用重心坐标反算)

    # 重心坐标函数
    def barycentric_coords(px, py):
        v0 = verts[2] - verts[0]
        v1 = verts[1] - verts[0]
        v2 = np.array([px, py]) - verts[0]

        d00 = np.dot(v0, v0)
        d01 = np.dot(v0, v1)
        d11 = np.dot(v1, v1)
        d20 = v0[0]*v2[0] + v0[1]*v2[1]
        d21 = v1[0]*v2[0] + v1[1]*v2[1]
        denom = d00 * d11 - d01 * d01
        v = (d11 * d20 - d01 * d21) / denom
        w = (d00 * d21 - d01 * d20) / denom
        u = 1 - v - w
        return u, v, w

    # 逐点计算重心坐标
    for i in range(res):
        for j in range(res):
            if mask[i,j]:
                u,v,w = barycentric_coords(xv[i,j], yv[i,j])
                r[i,j] = u
                g[i,j] = v
                b[i,j] = w

    # 只保留三角形内点
    img = np.zeros((res, res, 3))
    img[mask] = np.stack([r[mask], g[mask], b[mask]], axis=1)
    img[~mask] = 1  # 三角形外背景白色

    # 显示图片，extent对应坐标
    ax.imshow(img, origin='lower', extent=[0,1,0,np.sqrt(3)/2], interpolation='bilinear')

    # 画三角形边框
    polygon = Polygon(verts, closed=True, fill=False, edgecolor='black', linewidth=1)
    # ax.add_patch(polygon)

    # 标注顶点
    ax.text(0.5, np.sqrt(3)/2 + 0.05, 'Nb', color='red', ha='center', va='bottom', fontdict={'family': 'Arial', 'size': 8})
    ax.text(1.05, 0, 'O', color='green', ha='left', va='center', fontdict={'family': 'Arial', 'size': 8})
    ax.text(-0.05, 0, 'I', color='blue', ha='right', va='center', fontdict={'family': 'Arial', 'size': 8})
    
    ax.set_xlim(-0.1, 1.1)
    ax.set_ylim(-0.1, np.sqrt(3)/2 + 0.1)
    ax.axis('off')
    ax.set_aspect('equal')

ax_legend = fig.add_axes([0.9, 0.7, 0.2, 0.1])  # [left, bottom, width, height]

rgb_triangle_legend(ax_legend)

# ax.text(-0.1,1.03,'(a)',{'family':'Arial','weight':'normal','size':8},transform=ax.transAxes)
ax.set_title('0 GPa',font)

fig=plt.gcf()
cm=12
inch=0.3937007874
fig.set_size_inches(cm*inch,8.5*inch*1.5)
plt.savefig('band_color.png',format='png',dpi=600,bbox_inches='tight')