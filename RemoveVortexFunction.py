"""
Remove a TC vortex from the wind field to compute environmental parameters such as steering flow and wind shear.
For best results, ensure the input wind field grids are *larger* than the declared radius for calculation.

Source: Winterbottom and Chassignet 2011, https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011MS000088

Input variables:
    U -- the original 2-D or 3-D zonal wind field (m/s)
    V -- the original 2-D or 3-D meridional wind field (m/s)
    dist -- 2-D field of distances from the estimated TC center (km)
    r -- radius over which to remove the TC vortex (km)

Output variables:
    Unone -- the 2-D or 3-D zonal wind field with the TC vortex removed
    Vnone -- the 2-D or 3-D meridional wind field with the TC vortex removed
"""

def removeTCvortex(U, V, dist, r):
    import numpy as np
    
    rows,cols = np.where(dist<=r)  # remove vortex within 'r' km of TC center (WC2011 used 550 km)
    Unone = U.copy()
    Unone[:] = np.nan
    Vnone = Unone.copy()
    
    shape = U.shape
    if len(shape) == 2:  # input grid is 2-D
        umod = U.copy()  # zonal wind field for later modification
        vmod = V.copy()  # meridional wind field for later modification
        for n in range(len(rows)):
            ## First instance of appling 9-point spatial smoothing operator (WC2011, page 3)
            a = rows[n]
            b = cols[n]
            umod[a,b] = np.mean(U[a-1:a+2,b-1:b+2])
            vmod[a,b] = np.mean(V[a-1:a+2,b-1:b+2])
        for k in range(1500):  # very large number of iterations that should never be reached
            umod_old = umod.copy()
            vmod_old = vmod.copy()
            for n in range(len(rows)):
                a = rows[n]
                b = cols[n]
                ## Apply 9-point spatial smoothing operator
                umod[a,b] = np.mean(umod_old[a-1:a+2,b-1:b+2])
                vmod[a,b] = np.mean(vmod_old[a-1:a+2,b-1:b+2])
            if abs(np.var(umod_old[dist<=r])-np.var(umod[dist<=r])) < 1e-4:
                if abs(np.var(vmod_old[dist<=r])-np.var(vmod[dist<=r])) < 1e-4:
                    ## Stop iteration when zonal and meridional convergence values are reached (WC2011, page 3)
                    break
        Unone = umod.copy()
        Vnone = vmod.copy()
    elif len(shape) == 3:  # input grid is 3-D
        for l in range(shape[0]):  # iterate through available levels (THIS ASSUMES SAME CENTER AT ALL LEVELS!)
            umod = U[l,:,:].copy()
            vmod = V[l,:,:].copy()
            for n in range(len(rows)):  # do initial average
                a = rows[n]
                b = cols[n]
                ## First instance of appling 9-point spatial smoothing operator 
                umod[a,b] = np.mean(U[l,a-1:a+2,b-1:b+2])
                vmod[a,b] = np.mean(V[l,a-1:a+2,b-1:b+2])
            for k in range(1500):  # very large number of iterations that should never be reached
                umod_old = umod.copy()
                vmod_old = vmod.copy()
                for n in range(len(rows)):
                    a = rows[n]
                    b = cols[n]
                    ## Apply 9-point spatial smoothing operator (WC2011, page 3)
                    umod[a,b] = np.mean(umod_old[a-1:a+2,b-1:b+2])
                    vmod[a,b] = np.mean(vmod_old[a-1:a+2,b-1:b+2])
                if abs(np.var(umod_old[dist<=r])-np.var(umod[dist<=r])) < 1e-4:
                    if abs(np.var(vmod_old[dist<=r])-np.var(vmod[dist<=r])) < 1e-4:
                        ## Stop iteration when zonal and meridional convergence values are reached (WC2011, page 3)
                        break
            Unone[l,:,:] = umod.copy()
            Vnone[l,:,:] = vmod.copy()
    return Unone, Vnone


