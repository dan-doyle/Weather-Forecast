# Weather-Forecast
This project consists of running a WRF ([Weather Research and Forecasting](https://www.mmm.ucar.edu/weather-research-and-forecasting-model)) model using Orr2, University College Dublin’s High Performance Computing Cluster. Please see the full forecast report [here](https://github.com/dan-doyle/Weather-Forecast/blob/main/Forecast%20report.pdf).

## Brief Forecast Methodology
- Download global forecast data from the GFS global forecast model through [UCAR](https://rda.ucar.edu/index.html?hash=data_user&action=register)
- Load and use WPS, a module for WRF pre-processing ([see more information](https://github.com/wrf-model/WPS/blob/master/README))
    - Update the [namelist.wps](https://github.com/dan-doyle/Weather-Forecast/blob/main/namelist.wps) file to have correct date / time of forecast and longitude / latitude of weather station. Here we chose [London City Airport](https://rp5.ru/Weather_in_London_City_(airport)) ([more information about namelist.wps](https://www2.mmm.ucar.edu/wrf/users/namelist_best_prac_wps.html))
    - Link GFS data to the WPS directory 
    - Within WPS run three programs: 1. Geogrid.exe 2. UNGRIB.exe  3. Metgrid.exe
- Edit the [namelist.input](https://github.com/dan-doyle/Weather-Forecast/blob/main/namelist.input) file to have correct start and end date for forecast and the time interval at which to output forecast data (here we chose hourly) 
- Run the initialisation program Real.exe. Following this, submit the WRF model forecast task to the Orr2 queue ([information on commands](https://www2.mmm.ucar.edu/wrf/OnLineTutorial/index.php))
- Once finished, a wrfout.nc file is created with a forecast output. The wrfout.nc file is too large to upload to the repository however, we see exploratory analysis of the file below


## Exploratory Analysis
```python
import numpy as np
import xarray as xr

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import matplotlib.pyplot as plt
%matplotlib inline
```


```python
WRFdat = xr.open_dataset('/Users/daniel/wrfout.nc')
```




```python
# We use the Cartopy library to plot data https://scitools.org.uk/cartopy/docs/latest/
# When creating axes using plt.axes, we pass the name of the Cartopy projection: ccrs.Mercator()

fig = plt.figure(figsize=[8,6])
ax = plt.axes(projection=ccrs.Mercator())

ax.coastlines(resolution='50m', color='black',linewidth=0.5)

ax.set_extent([lons.values.min(),lons.values.max(), lats.values.min()-2, lats.values.max()])
# WRF temperature is set to Kelvin, we subtract 273.15 to set to C
cf = ax.contourf(lons, lats, t2m-273.15, transform=ccrs.PlateCarree())
# When plotting WRF data, we need to tell Cartopy that it is on its own grid. The transorm=ccrs.PlateCarree() command does this
plt.colorbar(cf, shrink=0.75)

grd = ax.gridlines(draw_labels=True, linestyle='--')
grd.xlabels_top=False
grd.ylabels_right=False
# Temperate plot (Celsius)
```


    
![png](Images/output_6_0.png)
    



```python
fig = plt.figure(figsize=[8,6])
ax = plt.axes(projection=ccrs.Mercator())

ax.coastlines(resolution='50m',color='black',linewidth=0.5)

ax.set_extent([lons.values.min(),lons.values.max(), lats.values.min()-2, lats.values.max()])
            
t_cl = ax.contour(lons,lats,t2m-273.15,np.array([-2,0,2,4,6,8]),
                 colors='red',linestyles='-',transform=ccrs.PlateCarree())

plt.clabel(t_cl,inline=1,fontsize=20,fmt='%1.0f',inline_spacing=1,colors='black')

#grd = ax.gridlines(draw_labels=True,linestlye='--')
# Weird error here, why do we have to take out linestyle?
grd = ax.gridlines(draw_labels=True)
grd.xlabels_top = False
grd.ylabels_right = False
# Here we try plotting the temperature (Celsius) with contour lines however, we see this gets messy
```


    
![png](Images/output_7_0.png)
    



```python
# Advanced plots 1

# First we define our variables
# In doing so we choose the hour in which we take our data from, we choose 0 for the first hour of the day
t2m = WRFdat.T2[0,:,:]   # 2-metre temperature
z = WRFdat.HGT[0,:,:]    # height above sea level of the ground
surfp = WRFdat.PSFC[0,:,:]  # atmospheric pressure at the ground
```


```python
surft = t2m + (6.5*z/1000)  
mslp = surfp*np.exp(9.81/(287.0*surft)*z)*0.01 + (6.7*z/ 1000)   
# The MSLP variable is the Mean Sea Level Pressure variable which is used to show
# regions of low and high pressure
```


```python
fig = plt.figure(figsize=[8,6])   
ax = plt.axes(projection=ccrs.Mercator())

ax.coastlines(resolution='50m',color='black',linewidth=0.5)  

ax.set_extent([lons.values.min(),lons.values.max(), lats.values.min()-2, lats.values.max()])

mslp_cl = ax.contour(lons,lats,mslp,np.array([990,995,1000,1005,1010,1015,1020,1025,1028,1030,1035]),
                 colors='black',linestyles='-',transform=ccrs.PlateCarree())

plt.clabel(mslp_cl,inline=1,fontsize=10,fmt='%1.0f',inline_spacing=1,colors='black')

#grd = ax.gridlines(draw_labels=True,linestlye='--')
# Weird error here, why do we have to take out linestyle?
grd = ax.gridlines(draw_labels=True)
grd.xlabels_top = False
grd.ylabels_right = False
# Plotting Contours of Mean Sea Level Pressure (MSLP) 
```


    
![png](Images/output_10_0.png)
    



```python
# Advanced plots 2

# First we define our variables
U = WRFdat.U10[0,:,:]
V = WRFdat.V10[0,:,:]

Wind_speed = np.sqrt(U**2 +V**2)

# Calculate area of max windspeed
iwsmx, jwsmx = np.unravel_index(Wind_speed.argmax(), Wind_speed.shape)
# iwsmx is the row the max value takes in the array and jwsmx is the column
# By running this code we avoid flatten the array, find the max and then put the array back together to find the index of the max value
# np.where(Wind_speed==np.amax(Wind_speed))
```


```python
fig = plt.figure(figsize=[8,6])
ax = plt.axes(projection=ccrs.Mercator())

ax.coastlines(resolution='50m', color='black',linewidth=0.5)

ax.set_extent([lons.values.min(),lons.values.max(), lats.values.min()-2, lats.values.max()])

NAMEnow = ax.contourf(lons, lats, Wind_speed, transform=ccrs.PlateCarree())  # What does PlateCarree do?

plt.plot(lons[iwsmx,:], lats[jwsmx,:], 'r--', transform=ccrs.PlateCarree())
# How does the above work because its variable about latitude too?

plt.colorbar(NAMEnow,shrink=0.75, cmap='jet')  #why doesn't jet work here?

grd = ax.gridlines(draw_labels=True, linestyle='--')
grd.xlabels_top=False
grd.ylabels_right=False

plt.title("Transect for cross-section plot")
# Plot 10-metre wind speed and transect through the location of maximum wind
```




    Text(0.5, 1.0, 'Transect for cross-section plot')




    
![png](Images/output_12_1.png)
    



```python
# Advanced plots 3

# Defining variables 
# First we note that wind speeds are calculated at staggered grids
# We must calculate wind speed variables that are 'unstaggered'
u = WRFdat.U[0,:,:]
v = WRFdat.V[0,:,:]
um = 0.5*(u[:,:,:-1] + u[:,:,1:])
vm = 0.5*(v[:,:-1,:] + v[:,1:,:])

# um and vm variavles represent the wind speed in the middle of each cell, so all values are stored in the same place now

ph = WRFdat.PH[0,:,:,:]       # These two variables proxy height
phb = WRFdat.PHB[0,:,:,:]

ApproxElev = 0.5*(phb[:-1,:,:]+ph[:-1,:,:]+phb[1:,:,:]+ph[1:,:,:])/9.81

# How do we incorporate the max wind speed values in..

#Now calculating the absolute value of the wind speeds
m = np.sqrt(um.values**2 + vm.values**2)

# new lons array which has matching dimensions
lons_new = np.tile(lons[iwsmx,:], (32,1))
```


```python
# Plot a vertical cross-section of wind speed along the transect shown in the above plot
fig = plt.figure(figsize=[8,6])
ax = plt.axes()

cf_three = ax.contourf(lons_new, ApproxElev[:,iwsmx,:], m[:,iwsmx,:])

plt.colorbar(cf_three,shrink=0.75)  

# Label axis
plt.title("Cross-section through max 10m wind speed transect")
plt.xlabel("longitude")
plt.ylabel("Elevation (m)")
```




    Text(0, 0.5, 'Elevation (m)')




    
![png](Images/output_14_1.png)
