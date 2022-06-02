#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 16 15:30:53 2020
Expanded grid and automatic download
Mon Jan 17 16:33:54 CET 2022

@author: Andrea Hahmann, DTU Wind
"""
import sys
import xarray as xr
import numpy as np
from wrf import interplevel
from datetime import datetime,timedelta
import cftime
from glob import glob
from pandas import to_datetime

from pyesgf.search import SearchConnection
import future_wind as fw

def read_and_interp(orog,dps,dts,dqs,du,dv,times):

    lev = [50.,100.,200.]
    levs = slice(20.,300.)

    # This data is in the 3hr file
    # should be subsampled to 6hr
    times_6hr = du.time.sel(time=times)
    ts = fw.combine_hemispheres(dts.tas,time=times_6hr)
    qs = fw.combine_hemispheres(dqs.huss,time=times_6hr)
    ps = fw.combine_hemispheres(dps.ps,time=times_6hr)
    rho = fw.air_density(ps,ts,qs)

    u = fw.combine_hemispheres(du.ua,time=times,lev=levs,for_interp=True)
    v = fw.combine_hemispheres(dv.va,time=times,lev=levs,for_interp=True)
    print("Times read",u.time[0].dt.strftime("%Y-%m-%d_%H").values,
          u.time[-1].dt.strftime("%Y-%m-%d_%H").values)

    u_mid = u.interp(lat=orog.lat,lon=orog.lon)
    v_mid = v.interp(lat=orog.lat,lon=orog.lon)
    
    wspd = np.sqrt(u_mid*u_mid + v_mid*v_mid) 
    ws2,z = xr.broadcast(wspd,u.lev)  ##z = a + (b-1.) * orog  # height relative to ground
    aux = interplevel(wspd, np.log(z), np.log(lev), meta=False)
    ws_z = xr.DataArray(aux, \
        coords=[z.time,lev,z.lat,z.lon], 
        dims=['time','level','lat','lon'])

    u_z = interplevel(u_mid, z, lev)
    v_z = interplevel(v_mid, z, lev)
    wd_z = np.rad2deg(np.arctan2(u_z,v_z)) + 180.
    filename = "wspd_wdir_"+u.time[0].dt.strftime("%Y-%m-%d_%H").values+".nc"
    
    return ws_z,wd_z,rho,filename

def main():

    # Different servers might be available at different 
    # times
    # server = 'https://esgf-node.llnl.gov/esg-search'
    server = 'https://esgf-data.dkrz.de/esg-search'
    # server = 'https://esgf-index1.ceda.ac.uk/esg-search'
    conn = SearchConnection(server, distrib=True)

    # model = 'HadGEM3-GC31-MM'
    model = 'HadGEM3-GC31-LL'
    calendar = '360_day'
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
    oro_file = fw.search_esgf(
        conn,"orog",model,'hist-1950','r1i1p1f1',table='fx')
    print("Open ORO dataset",oro_file)
    ds = xr.open_dataset(oro_file)
    orog = fw.combine_hemispheres(ds.orog)

    filenames = "wspd_wdir_????-??-??_??.nc"
    old_files = sorted(glob(filenames))

    if not old_files:
        print("No previous files")
        date = cftime.datetime(year,1,1,6,calendar=calendar) # Files start at 06 not 00
    else:
        ff = xr.open_dataset(old_files[-1],decode_times=True,use_cftime=True)
        date = ff.time[-1] + timedelta(hours=6)
        print("Next date:",date.values)
        year = date.dt.year.values
        date = fw.datetime_to_cftime(date,calendar=calendar) # New date to start download

    nn = 1   # (=360*4/20); 20 files per year
    last_date = cftime.datetime(last_year+1,1,1,0,calendar=calendar)
    print("Last date:",last_date)

    while (year <= last_year):

        print(date)
        date_end = date + timedelta(days=nn)
        print(date_end)

        u_file = fw.search_esgf(
            conn,"ua",model,experiment,variant,date=date)
        print("Open U datasets",u_file)
        
        # First start and end dates in file from the filename
        dates = u_file.split("_")[-1].split("-")[0][0:10]
        start_time_in_file = fw.date_to_cftime(dates,calendar=calendar)
        dates = u_file.split("_")[-1].split("-")[1][0:10]
        end_time_in_file = fw.date_to_cftime(dates,calendar=calendar)
        print("Dates in file:",start_time_in_file.strftime(),\
            end_time_in_file.strftime())

        print("Process data from",date,"to",date_end)
        du = xr.open_dataset(u_file)

        ps_file = fw.search_esgf(
            conn,"ps",model,experiment,variant,date=date,table='3hr')
        print("Open PS dataset",ps_file)
        dps = xr.open_dataset(ps_file)

        v_file = fw.search_esgf(
            conn,"va",model,experiment,variant,date=date)
        print("Open V dataset",v_file)
        dv = xr.open_dataset(v_file)

        ts_file = fw.search_esgf(
            conn,"tas",model,experiment,variant,date=date,table='3hr')
        print("Open TS dataset",ts_file)
        dts = xr.open_dataset(ts_file)

        qs_file = fw.search_esgf(
            conn,"huss",model,experiment,variant,date=date,table='3hr')
        print("Open QS dataset",qs_file)
        dqs = xr.open_dataset(qs_file)

        while (date <= min([end_time_in_file,last_date])):

            ws,wd,rho,filename = read_and_interp(
                orog,dps,dts,dqs,du,dv,slice(date,date_end))
            ds = fw.make_data_set(du,ws,wd,rho)
            ds.to_netcdf(filename,mode="w",engine="netcdf4",
                        unlimited_dims='time')
            print(filename," written to disk")
            sys.exit()

            date = date_end + timedelta(hours=6)
            print("Next date:",date)

        year = date.year
        print("Next date:",date.strftime("%Y-%m-%d_%H"),"year:",year)

if __name__ == "__main__":
    main()