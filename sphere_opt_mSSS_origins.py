# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 13:21:45 2023

@author: xanmc

Optimise one sphere and two sphere origins found from "sphere_opt.py" written by Eric Larson Jan 2023
Print the two-sphere origin case for use in MATLAB creation of mSSS basis
Create visuals used for manuscript
"""

import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import mne
import nibabel as nib
import scipy
from scipy.spatial import KDTree
import vedo

mindist = 2e-3

data_path = mne.datasets.sample.data_path()
subject = 'sample'
subjects_dir = data_path / 'subjects'
fname_bem = subjects_dir / subject / 'bem' / f'{subject}-5120-5120-5120-bem.fif'
trans = data_path / 'MEG' / 'sample' / 'sample_audvis_raw-trans.fif'

sample_data_folder = mne.datasets.sample.data_path()
sample_data_raw_file = os.path.join(data_path, 'MEG', 'sample','sample_audvis_raw.fif')
raw = mne.io.read_raw_fif(sample_data_raw_file, verbose=False)
channel_pos=raw._get_channel_positions(picks='meg')
mri_head_t = mne.transforms.invert_transform(mne.read_trans(trans)) #trans is the raw_fif_trans file
assert mri_head_t['from'] == mne.io.constants.FIFF.FIFFV_COORD_MRI, mri_head_t['from']

##########################################
## specify sensor positions for visual
trans_channel_pos = np.zeros([305,3])
mri_head_t_matrix=mri_head_t['trans']
add_to=[0,-0.01,0.03]
multiply_by=mri_head_t_matrix[:3,:3]
for i in range(0,305): #plot one side hemisphere of helmet
    trans_channel_pos[i]=np.matmul(multiply_by,np.transpose(channel_pos[i]))
for i in range(0,305):
    for j in range(0,3):
        trans_channel_pos[i,j]=channel_pos[i,j]+add_to[j]   
        
################# 
# create BEM brain surfaces
bem_surf = mne.read_bem_surfaces(fname_bem) 
assert bem_surf[0]['id'] == mne.io.constants.FIFF.FIFFV_BEM_SURF_ID_HEAD
assert bem_surf[2]['id'] == mne.io.constants.FIFF.FIFFV_BEM_SURF_ID_BRAIN
scalp, _, inner_skull = bem_surf
inside_scalp = mne.surface._CheckInside(scalp, mode='pyvista')
inside_skull = mne.surface._CheckInside(inner_skull, mode='pyvista')
m3_to_cc = 100 ** 3
assert inside_scalp(inner_skull['rr']).all()
assert not inside_skull(scalp['rr']).any()
b = vedo.Mesh([inner_skull['rr'], inner_skull['tris']])
s = vedo.Mesh([scalp['rr'], scalp['tris']])
s_tree = KDTree(scalp['rr'])
brain_volume = b.volume()
print(f'Brain vedo:     {brain_volume * m3_to_cc:8.2f} cc')
brain_vol = nib.load(subjects_dir / subject / 'mri' / 'brainmask.mgz')
brain_rr = np.array(np.where(brain_vol.get_fdata())).T
brain_rr = mne.transforms.apply_trans(brain_vol.header.get_vox2ras_tkr(), brain_rr) / 1000. #apply a transformation matrix
del brain_vol
brain_rr = brain_rr[inside_skull(brain_rr)]
vox_to_m3 = 1e-9
brain_volume_vox = len(brain_rr) * vox_to_m3

def _print_q(title, got, want):
    title = f'{title}:'.ljust(15)
    print(f'{title} {got * m3_to_cc:8.2f} cc ({(want - got) / want * 100:6.2f} %)')

_print_q('Brain vox', brain_volume_vox, brain_volume_vox)

# 1. Compute a naive sphere using the center of mass of brain surf verts
naive_c = np.mean(inner_skull['rr'], axis=0)
naive_r = np.min(np.linalg.norm(inner_skull['rr'] - naive_c, axis=1))
naive_v = 4 / 3 * np.pi * naive_r ** 3
_print_q('Naive sphere', naive_v, brain_volume)
s1 = vedo.Sphere(naive_c, naive_r, res=100)
_print_q('Naive vedo', s1.volume(), brain_volume)

# 2. Now use the larger radius (to head) plus mesh arithmetic
better_r = s_tree.query(naive_c)[0] - mindist
s1 = vedo.Sphere(naive_c, better_r, res=24)
_print_q('Better vedo', s1.boolean("intersect", b).volume(), brain_volume)
v = np.sum(np.linalg.norm(brain_rr - naive_c, axis=1) <= better_r) * vox_to_m3
_print_q('Better vox', v, brain_volume_vox)

# 3. Now optimize one sphere
from scipy.optimize import fmin_cobyla #constrained optimization by linear approximation

def _cost(c):
    cs = c.reshape(-1, 3)
    rs = np.maximum(s_tree.query(cs)[0] - mindist, 0.)
    resid = brain_volume
    mask = None
    for c, r in zip(cs, rs):
        if not (r and s.contains(c)):
            continue
        m = np.linalg.norm(brain_rr - c, axis=1) <= r
        if mask is None:
            mask = m
        else:
            mask |= m
    resid = brain_volume_vox
    if mask is not None:
        resid = resid - np.sum(mask) * vox_to_m3
    return resid

def _cons(c):
    cs = c.reshape(-1, 3)
    sign = np.array([2 * s.contains(c) - 1 for c in cs], float)
    cons = sign * s_tree.query(cs)[0] - mindist
    return cons

x = naive_c
c_opt_1 = fmin_cobyla(_cost, x, _cons, rhobeg=1e-2, rhoend=1e-4)
v_opt_1 = brain_volume_vox - _cost(c_opt_1)
_print_q('COBYLA 1', v_opt_1, brain_volume_vox)

# 4. Now optimize two spheres
x = np.concatenate([c_opt_1, naive_c])
c_opt_2 = fmin_cobyla(_cost, x, _cons, rhobeg=1e-2, rhoend=1e-4)
v_opt_2 = brain_volume_vox - _cost(c_opt_2)
_print_q('COBYLA 2', v_opt_2, brain_volume_vox)

# 4. Finally, three spheres (not perfect, not global opt)
x = np.concatenate([c_opt_2, naive_c])
c_opt_3 = fmin_cobyla(_cost, x, _cons, rhobeg=1e-2, rhoend=1e-4)
v_opt_3 = brain_volume_vox - _cost(c_opt_3)
_print_q('COBYLA 3', v_opt_3, brain_volume_vox)

###########################################
# create visualization of on-scalp system
import pyvista as pv
import pyvistaqt
## make my own points
cart_points=np.zeros([9,3])
r=0.13;
rho=1;
theta=[-5/10*np.pi,-4/10*np.pi,-3/10*np.pi,-2/10*np.pi,-1/10*np.pi,0,1/10*np.pi,2/10*np.pi,3/10*np.pi,4/10*np.pi]
for i in range(0,9):
    theta_i = theta[i]
    cart_points[i,2] = r*np.sin(rho)*np.cos(theta_i) #x
    cart_points[i,1] = r*np.sin(rho)*np.sin(theta_i) #y
    cart_points[i,0] = r*np.cos(rho) #z

################
# Visualize VSH basis as spheres 
plotter = pyvistaqt.BackgroundPlotter(
      shape=(1, 2), window_size=(1200, 300), #(1,3) for three sphere case
      editor=False, menu_bar=False, toolbar=False)
plotter.background_color = 'w'
brain_mesh = pv.make_tri_mesh(inner_skull['rr'], inner_skull['tris'])
scalp_mesh = pv.make_tri_mesh(scalp['rr'], scalp['tris'])
colors = matplotlib.rcParams['axes.prop_cycle'].by_key()['color']
mesh_kwargs = dict(render=False, reset_camera=False, smooth_shading=True)
for ci, cs in enumerate((c_opt_1, c_opt_2)):
    plotter.subplot(0, ci)
    plotter.camera.position = (0., -0.5, 0)
    plotter.camera.focal_point = (0., 0., 0.)
    plotter.camera.azimuth = 90
    plotter.camera.elevation = 0
    plotter.camera.up = (0., 0., 1.)
    plotter.add_mesh(brain_mesh, opacity=0.2, color='k', **mesh_kwargs)
    plotter.add_mesh(scalp_mesh, opacity=0.1, color='tan', **mesh_kwargs)
    plotter.add_points(cart_points, opacity=0.8, color='b', point_size=40, render_points_as_spheres=True)
    for c, color in zip(cs.reshape(-1, 3), colors):
        sphere = pv.Sphere(s_tree.query(c)[0] - mindist, c)
        plotter.add_mesh(sphere, opacity=0.5, color=color, **mesh_kwargs)
plotter.show()

### Output origins of 2-origin VSH case
# These centers are used for the mSSS basis
for use in (c_opt_1, c_opt_2):
    centers = mne.transforms.apply_trans(mri_head_t, use.reshape(-1, 3))
    print(centers)
    
