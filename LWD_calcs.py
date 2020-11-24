# -*- coding: utf-8 -*-
"""
Created on Mon Nov 23 21:37:14 2020

@author: jamor
"""

import numpy as np

####################################################################################################
####################################################################################################
# FUNCTIONS
####################################################################################################
####################################################################################################
def LWD_calcs(v_1, LWD, BLDR, DEP, CONST):
    if np.ndim(LWD['Diam_in'])>1:
        Diam_in = np.squeeze(LWD['Diam_in'], axis=0)
    else:
        Diam_in = LWD['Diam_in']
    if np.ndim(LWD['Len_ft'])>1:
        Len_ft = np.squeeze(LWD['Len_ft'], axis=0)
    else:
        Len_ft = LWD['Len_ft']
    ##################################
    # Boulder component
    # volume of boulder (ft^3)
    Vol_ft3 = (BLDR['Diam_in']/12)**3
    # dry weight of boulder (lbf)
    W_dry_lbf = CONST['gam_bldr']*Vol_ft3
    # submerged weight of boulder (lbf)
    W_wet_lbf = (CONST['gam_bldr']-CONST['gam_wat'])*Vol_ft3
    # weight of boulders (lbf)
    W_bldr_lbf = BLDR['Num']*BLDR['subm']*W_wet_lbf + BLDR['Num']*(1-BLDR['subm'])*W_dry_lbf
    ##################################
    # Member component
    D_rw_in = 2*Diam_in
    L_rw_ft = 1.5*Diam_in/12
    # XS area of stem (ft^2)
    A_log_ft2 = np.pi*np.square(Diam_in/12)/4
    # XS area of root wad (ft^2)
    A_rw_ft2 = np.pi*np.square(D_rw_in/12)/4
    # total volume of member (ft^3)
    Vol_ft3 = np.multiply(Len_ft-L_rw_ft, A_log_ft2)+np.multiply(L_rw_ft, A_rw_ft2)
    # total weight of member (lbf)
    F_G_lbf = Vol_ft3*CONST['gam_wood']
    # buoyant force (lbf), weight of water displaced by member less the weight of the member
    F_B_lbf = Vol_ft3*CONST['gam_wat']-F_G_lbf
    ##################################
    # Drag Calculations for single member
    # area of log normal to flow (ft^2)
    A_drag_ft2 = np.sin(np.deg2rad(LWD['alpha']))*(1-LWD['R_embed'])*\
                 np.multiply(Len_ft-L_rw_ft,Diam_in/12)+\
                 np.multiply(L_rw_ft,D_rw_in/12)
    # drag force on log normal to flow (lbf)
    F_D_log_lbf = (v_1**2/2/CONST['g'])*A_drag_ft2*LWD['subm']*CONST['C_DP']*CONST['gam_wat']
    # drag force on rootwad normal to flow (lbf)
    F_D_rw_lbf = (v_1**2/2/CONST['g'])*A_rw_ft2*LWD['subm']*CONST['C_DU']*CONST['gam_wat']
    # maximum drag force (lbf)
    F_D_single_lbf = np.maximum(F_D_log_lbf, F_D_rw_lbf)*LWD['Num']
    # volume of soil on embedded log (ft^3)
    Vol_soil_embed_ft3 = DEP['soil']*LWD['R_embed']*1.5*np.multiply(Len_ft,Diam_in/12)
    # depth of dry soil (ft), UNUSED SO COMMENTED OUT
    # D_soil_dry_ft = max(D_soil_feet+D_jam_feet-D_water_feet, 0)
    # depth of submerged soil (ft)
    D_soil_wet_ft = min(DEP['soil'], DEP['water']-LWD['dep_ft'])
    # weight of soil on embedded log (lbf)
    W_log_soil_lbf = Vol_soil_embed_ft3*\
                     (CONST['gam_soil']+D_soil_wet_ft/DEP['soil']*(CONST['gam_soil_sat']-\
                      CONST['gam_wat']-CONST['gam_soil']))
    # friction force on log (lbf)
    F_F_single_lbf = CONST['mu_gravel']*(W_bldr_lbf+(W_log_soil_lbf-F_B_lbf)*LWD['Num'])
    ##################################
    # Summary of Results for Single Log
    # uplift factor of safety
    FS_U_single = (W_bldr_lbf+W_log_soil_lbf*LWD['Num'])/LWD['Num']/F_B_lbf
    # drag factor of safety
    FS_D_single = F_F_single_lbf/F_D_single_lbf
    ##################################
    if ((np.ndim(Diam_in)==1) and (np.ndim(Len_ft)==1)):
        if W_bldr_lbf>0:
            print('Weight of boulders: W_Boulders = '+str(np.round(W_bldr_lbf/2000,2))+' tonf')
        print('Gravity force on log: F_G = '+str(np.round(F_G_lbf/2000,2))+' tonf')
        print('Buoyant force on log: F_B = '+str(np.round(F_B_lbf/2000,2))+' tonf')
        print('Drag force on log: F_D = '+str(np.round(F_D_single_lbf/2000,2))+' tonf')
        print('Friction force on log: F_F = '+str(np.round(F_F_single_lbf/2000,2))+' tonf')
    ##################################
    return FS_U_single, FS_D_single