import cartopy.crs as ccrs       # cartopy 0.18.0
import cartopy.feature as cfeature

import matplotlib.pyplot as plt  # matplotlib 3.3.0
import matplotlib.ticker as mticker
import numpy as np               # numpy 1.19.4
import xarray as xr              # xarray 0.16.2

def shear_map(x, savefile = None):
    tclat = x.attrs["center_lat"]
    tclon = x.attrs["center_lon"]
    if tclon > 180:
        tclon = tclon - 360

    lat_vals = x.lat.values
    lon_vals = x.lon.values

    glat, glon = [y.T for y in np.meshgrid(lat_vals, lon_vals)]
    glon[glon>180.] -= 360.

    u = x.sel(component = "u").values
    v = x.sel(component = "v").values
    mag = x.sel(component = "magnitude").values

    contour_array = np.arange(5.,40.1,1.)
    colorbar_array = np.arange(5.,40.1,5.)

    plt.ioff()
    
    ## Add lat/lon lines and ensure labels occur within the domain at reasonable intervals
    def prettylines(clon, clat):
        gl = ax.gridlines(crs=pc, draw_labels=False, linewidth=0.5, linestyle='dashed')
        gl.xlocator = mticker.FixedLocator(np.arange(-180.,1.,3.))
        gl.ylocator = mticker.FixedLocator(np.arange(0.,61.,2.))
        gl2 = ax.gridlines(crs=pc, draw_labels=True, x_inline=False, linewidth=0.)
        gl2.right_labels, gl2.top_labels, gl2.rotate_labels = [False] * 3
        xmin, xmax = [np.ceil((clon-9.)/3.)*3., np.floor((clon+9.)/3.)*3.+1.]
        ymin, ymax = [np.ceil((clat-9.)/2.)*2., np.floor((clat+9.)/2.)*2.+1.]
        gl2.xlocator = mticker.FixedLocator(np.arange(xmin,xmax,3.))
        gl2.ylocator = mticker.FixedLocator(np.arange(ymin,ymax,2.))
        gl2.xlabel_style = {'size': 12}
        gl2.ylabel_style = {'size': 12}

    proj = ccrs.LambertConformal(central_longitude=tclon, central_latitude=tclat)
    pc = ccrs.PlateCarree()  # required for cartopy to plot dimensions of latitude and longitude

    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection=proj)
    ax.set_extent([np.min(glon),np.max(glon),np.min(glat),np.max(glat)], crs=pc)  # center domain on TC latitude & longitude
    ax.coastlines('50m', linewidth=2., zorder=9)
    im = ax.contourf(glon, glat, mag, contour_array, 
                    cmap=plt.get_cmap('YlOrRd'), extend='both', transform=pc, zorder=1)  # shading for shear magnitude
    pos = ax.get_position()
    left = pos.x0 - 0.01
    bottom = pos.y1 + 0.025
    width = pos.x1 - left
    height = 0.015
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=colorbar_array, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cbar.ax.tick_params(labelsize=10, pad=0)
    cbar.set_label('200-850-hPa shear magnitude [m/s]', size=10)
    lw = 4. * (mag / 40.)  # scale factor for width of streamlines based on wind speed [m/s]
    ax.streamplot(glon, glat, u, v, density=0.6, color='k', 
                linewidth=lw, transform=pc, zorder=5)
    ax.scatter(tclon, tclat, s=80, c='w', marker='X', linewidths=0.5, edgecolors='k', transform=pc, zorder=10)
    prettylines(tclon, tclat)

    shear_x = x.attrs['avg_shear'][0]
    shear_y = x.attrs['avg_shear'][1]
    sx, sy = 5*np.array([shear_x, shear_y])/np.linalg.norm([shear_x, shear_y])
    ax.arrow(tclon, tclat, sx, sy, transform = pc, color = "green", width = 0.15)

    if 'storm_direction' in x.attrs:
        vel_x = x.attrs['storm_direction'][0]
        vel_y = x.attrs['storm_direction'][1]
        vx, vy = 5*np.array([vel_x, vel_y])/np.linalg.norm([vel_x, vel_y])
        ax.arrow(tclon, tclat, vx, vy, transform = pc, color = "blue", width = 0.15)
    
    
    if savefile is None:
        plt.show()
    else:
        plt.savefig(savefile)

def two_shade_map(x, toplot, 
                  shading = np.arange(-3.,3.,.1), ticks = np.arange(-3.1,3.1,1.), 
                  savefile = None,
                  legend_title = "Integrated Circulation"):

    tclat = x.attrs["center_lat"]
    tclon = x.attrs["center_lon"]
    if tclon > 180:
        tclon = tclon - 360

    lat_vals = x.lat.values
    lon_vals = x.lon.values

    glat, glon = [y.T for y in np.meshgrid(lat_vals, lon_vals)]
    glon[glon>180.] -= 360.

    plt.ioff()

    ## Add lat/lon lines and ensure labels occur within the domain at reasonable intervals
    def prettylines(clon, clat):
        gl = ax.gridlines(crs=pc, draw_labels=False, linewidth=0.5, linestyle='dashed')
        gl.xlocator = mticker.FixedLocator(np.arange(-180.,1.,3.))
        gl.ylocator = mticker.FixedLocator(np.arange(0.,61.,2.))
        gl2 = ax.gridlines(crs=pc, draw_labels=True, x_inline=False, linewidth=0.)
        gl2.right_labels, gl2.top_labels, gl2.rotate_labels = [False] * 3
        xmin, xmax = [np.ceil((clon-9.)/3.)*3., np.floor((clon+9.)/3.)*3.+1.]
        ymin, ymax = [np.ceil((clat-9.)/2.)*2., np.floor((clat+9.)/2.)*2.+1.]
        gl2.xlocator = mticker.FixedLocator(np.arange(xmin,xmax,3.))
        gl2.ylocator = mticker.FixedLocator(np.arange(ymin,ymax,2.))
        gl2.xlabel_style = {'size': 12}
        gl2.ylabel_style = {'size': 12}

    proj = ccrs.LambertConformal(central_longitude=tclon, central_latitude=tclat)
    pc = ccrs.PlateCarree()  # required for cartopy to plot dimensions of latitude and longitude

    fig = plt.figure(figsize=(8,8))
    ax = plt.axes(projection=proj)
    ax.set_extent([np.min(glon),np.max(glon),np.min(glat),np.max(glat)], crs=pc)  # center domain on TC latitude & longitude
    ax.coastlines('50m', linewidth=2., zorder=9)
    im = ax.contourf(glon, glat, toplot, shading, 
                    cmap=plt.get_cmap('PiYG'), extend='both', transform=pc, zorder=1)  # shading 
    pos = ax.get_position()
    left = pos.x0 - 0.01
    bottom = pos.y1 + 0.025
    width = pos.x1 - left
    height = 0.015
    cax = fig.add_axes([left, bottom, width, height])
    cbar = plt.colorbar(im, cax=cax, ticks=ticks, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cbar.ax.tick_params(labelsize=10, pad=0)
    cbar.set_label(legend_title, size=10)
    prettylines(tclon, tclat)

    if savefile is None:
        plt.show()
    else:
        plt.savefig(savefile)