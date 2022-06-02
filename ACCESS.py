#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 11 10:22:53 2020
Expanded grid and automatic download
Sun Mar  6 10:06:42 CET 2022

@author: ahah
"""
import sys
import xarray as xr
import numpy as np
from wrf import interplevel
from datetime import datetime,timedelta
from glob import glob
import cftime
from pyesgf.search import SearchConnection
import future_wind as fw

def read_and_interp(orog,dps,dts,dqs,du,dv,times):  

    lev = [50.,100.,200.]
    levs = slice(20.,300.)

    u = fw.combine_hemispheres(du.ua,time=times,lev=levs,for_interp=True)
    v = fw.combine_hemispheres(dv.va,time=times,lev=levs,for_interp=True)
    print("Times read",u.time[0].dt.strftime("%Y-%m-%d_%H").values,
          u.time[-1].dt.strftime("%Y-%m-%d_%H").values)

    u_mid = u.interp(lat=orog.lat,lon=orog.lon)
    v_mid = v.interp(lat=orog.lat,lon=orog.lon)
    
    wspd = np.sqrt(u_mid*u_mid + v_mid*v_mid) 
    _,z = xr.broadcast(wspd,u.lev)  ##z = a + (b-1.) * orog  # height relative to ground
    aux = interplevel(wspd, np.log(z), np.log(lev), meta=False)
    ws_z = xr.DataArray(aux, \
        coords=[z.time,lev,z.lat,z.lon], 
        dims=['time','level','lat','lon'])
    
    u_z = interplevel(u_mid, z, lev)
    v_z = interplevel(v_mid, z, lev)
    wd_z = np.rad2deg(np.arctan2(u_z,v_z)) + 180.

    tas = dts.tas
    ts = fw.combine_hemispheres(tas,time=slice(wspd.time[0],wspd.time[-1]))
    # print(ts.time[0]-timedelta(days=1),ts.time[-1]+timedelta(days=1))
    
    qas = dqs.huss
    qs = fw.combine_hemispheres(qas,\
        time=slice(
            ts.time[0]-timedelta(days=1),
            ts.time[-1]+timedelta(days=2)))
    # print(qs.time)
    # interpolate to same 6-hour time as tas
    qs = qs.interp(time=ts.time)
    # if the first hour is missing, copy from the second in array
    hour = qs.time.dt.hour
    if (hour[0] == 6):
        qs[0,:,:] = qs.isel(time=1)
        
    ps = fw.combine_hemispheres(dps.ps,time=wspd.time)
    rho = fw.air_density(ps,ts,qs)

    filename = "wspd_wdir_"+u.time[0].dt.strftime("%Y-%m-%d_%H").values+".nc"
    
    return ws_z,wd_z,rho,filename

def main():

    # server = 'https://esgf-node.llnl.gov/esg-search'
    server = 'https://esgf-data.dkrz.de/esg-search'
    conn = SearchConnection(server, distrib=True)

    # model = 'ACCESS-CM2'
    model = 'ACCESS-ESM1-5' # also for ACCESS-CM2
    calendar = 'proleptic_gregorian'
    experiment = sys.argv[1]
    if (experiment == "historical"):
        year = 1980; last_year = 2014
    else:
        year = 2015; last_year = 2050
    variant = sys.argv[2]  
    print("Retrieve data for",\
        "\n model:  ",model,"\n experiment:",experiment,\
        "\n variant:",variant)

    # Download oro file for U and V interpolation from Arakawa-C grid
    oro_file = fw.search_esgf(conn,"orog",model,experiment,variant,table='fx')
    print("Open ORO dataset",oro_file)
    ds = xr.open_dataset(oro_file)
    orog = fw.combine_hemispheres(ds.orog)

    # What filenames already exist in the directory
    filenames = "wspd_wdir_????-??-??_??.nc"
    old_files = sorted(glob(filenames))

    if not old_files:    # This is necessary for the scenario files that start at 00Z
        print("No previous files")
        date = cftime.datetime(year,1,1,6,calendar=calendar) # Files start at 06 not 00
    else:
        ff = xr.open_dataset(old_files[-1],decode_times=True,use_cftime=True)
        date = ff.time[-1] + timedelta(hours=6)
        print("Next date:",date.values)
        year = date.dt.year.values
        date = fw.datetime_to_cftime(date,calendar=calendar) # New date to start download

    last_date = cftime.datetime(last_year+1,1,1,0,calendar=calendar)
    nn = 32   # number of days in one download. About a month is a good number of 

    while (year <= last_year):

        u_file = fw.search_esgf(conn,"ua",model,experiment,variant,date=date)
        print("U dataset found",u_file)
        # First start and end dates in file from the filename
        dates = u_file.split("_")[-1].split("-")[0][0:10]
        start_time_in_file = fw.date_to_cftime(dates,calendar=calendar)
        dates = u_file.split("_")[-1].split("-")[1][0:10]
        end_time_in_file = fw.date_to_cftime(dates,calendar=calendar)
        print("Dates in file:",start_time_in_file,end_time_in_file)

        # print("Open U datasets",u_file)
        du = xr.open_dataset(u_file,decode_times=True,use_cftime=True)

        v_file = fw.search_esgf(conn,"va",model,experiment,variant,date=date)
        print("Open V dataset",v_file)
        dv = xr.open_dataset(v_file,decode_times=True,use_cftime=True)

        ps_file = fw.search_esgf(conn,"ps",model,experiment,variant,date=date)
        print("Open PS dataset",ps_file)
        dps = xr.open_dataset(ps_file,decode_times=True,use_cftime=True)

        if (model == "ACCESS-CM2"):
            ts_file = fw.search_esgf(conn,"tas",model,experiment,variant,\
                date=date,table='3hr') 
            print("Open TS dataset",ts_file)
            dts = xr.open_dataset(ts_file,decode_times=True,use_cftime=True)
        else:
            ts_file = fw.search_esgf(conn,"tas",model,experiment,variant,\
                date=date,table='6hrPlevPt') 
            print("Open TS dataset",ts_file)
            dts = xr.open_dataset(ts_file,decode_times=True,use_cftime=True)

        if (model == "ACCESS-CM2"):
            qs_file = fw.search_esgf(conn,"huss",model,experiment,variant,\
                date=date,table='3hr')
            print("Open QS dataset",qs_file)
            dqs = xr.open_dataset(qs_file,decode_times=True,use_cftime=True)
        else:

            qs_file = fw.search_esgf(conn,"huss",model,experiment,variant,\
                date=date,table='day')
            print("Open QS dataset",qs_file)
            dqs = xr.open_dataset(qs_file,decode_times=True,use_cftime=True)

        while (date <= min([end_time_in_file,last_date])):

            date_end = date + timedelta(days=nn)
            date_end = min([date_end,end_time_in_file])
            ws,wd,rho,filename = read_and_interp(
                orog,dps,dts,dqs,du,dv,slice(date,date_end))

            ds = fw.make_data_set(du,ws,wd,rho)
            ds.to_netcdf(filename,mode="w",engine="netcdf4",
                        unlimited_dims='time')
            print(filename," written to disk")

            date = date_end + timedelta(hours=6)
            print("Next date:",date)
        
        year = date.year
        print("Next date:",date.strftime("%Y-%m-%d_%H"),"year:",year)

if __name__ == "__main__":
    main()
