# -*- coding: utf-8 -*-
"""
Created on Tue May 20 18:16:27 2025

@author: chenq4
"""

import os
import numpy as np
import gdal
import gdalconst
import datetime
import calendar
import pandas as pd


#other Functions
def get_month_dates(year, month):
    _, last_day = calendar.monthrange(year, month)
    return [datetime.date(year, month, d) for d in [10, 20, last_day]]

def generate_10day_dates(year):
    dates = []
    for month in range(1, 13):
        for d in get_month_dates(year, month):
            dates.append(datetime.datetime(d.year, d.month, d.day))
    return dates

def generate_daily_dates(year):
    return [datetime.datetime(year, 1, 1) + datetime.timedelta(days=i) for i in range(365)]

def interpolate_10day_to_daily(ts_10day, year):
    try:
        date_10 = generate_10day_dates(year)
        date_daily = generate_daily_dates(year)
        series = pd.Series(ts_10day, index=date_10)
        dailylai = series.reindex(series.index | date_daily).interpolate(method='index').loc[date_daily]
        dailylai=dailylai.fillna(method='bfill')
        return dailylai
    except Exception:
        return np.full((365,), np.nan)


# Main Processing Functions
def compute_abLAI_block(lai_block, baseline_indices, target_index):
    layers = np.array(np.split(lai_block, 15, axis=0))
    baseline = layers[baseline_indices, :, :, :]
    baseline_mean = np.nanmean(baseline, axis=0)
    baseline_std = np.nanstd(baseline, axis=0)
    target = layers[target_index, :, :, :]
    diff = target - baseline_mean
    with np.errstate(divide='ignore'):
        ablai = np.where(baseline_std != 0., diff / baseline_std, 0)
    return ablai

def calculate_abLAI(input_tif, output_tif, target_year_index, baseline_indices):
    #read data 
    src = gdal.Open(input_tif, gdalconst.GA_ReadOnly)
    proj = src.GetProjection()
    geotrans = src.GetGeoTransform()
    
    #write output data
    Nx, Ny = 2500, 4000
    drv = gdal.GetDriverByName("GTIFF")
    dst_ds = drv.Create(output_tif, Ny, Nx, 36, gdal.GDT_Float32, options=["INTERLEAVE=BAND", "TILED=YES", "BIGTIFF=YES", "COMPRESS=LZW"])
    dst_ds.SetGeoTransform(geotrans)
    dst_ds.SetProjection(proj)
    
    #read in blocks
    xsize, ysize = src.RasterXSize, src.RasterYSize
    block_xsize, block_ysize = 1000, 1000
    for x in range(0, xsize, block_xsize):
        if x + block_xsize < xsize:
            cols = block_xsize
        else:
            cols = xsize - x
        for y in range(0, ysize, block_ysize):
            if y + block_ysize < ysize:
                rows = block_ysize
            else:
                rows = ysize - y
            lai_block = src.ReadAsArray(x, y, cols, rows)
            print(x,y)
            ablai_block = compute_abLAI_block(lai_block, baseline_indices, target_year_index)
            #save 
            for c in range(36):
                dst_ds.GetRasterBand(c + 1).WriteArray(ablai_block[c, :, :], x, y)
    dst_ds = None

def interpolate_abLAI(input_tif, output_tif, year):
    src = gdal.Open(input_tif, gdalconst.GA_ReadOnly)
    proj = src.GetProjection()
    geotrans = src.GetGeoTransform()
    
    Nx, Ny = 2500, 4000
    drv = gdal.GetDriverByName("GTIFF")
    dst_ds = drv.Create(output_tif, Ny, Nx, 365, gdal.GDT_Float32, options=["INTERLEAVE=BAND", "TILED=YES", "BIGTIFF=YES", "COMPRESS=LZW"])
    dst_ds.SetGeoTransform(geotrans)
    dst_ds.SetProjection(proj)

    xsize, ysize = src.RasterXSize, src.RasterYSize
    block_xsize, block_ysize = 1000, 1000
    for x in range(0, xsize, block_xsize):
        if x + block_xsize < xsize:
            cols = block_xsize
        else:
            cols = xsize - x
        for y in range(0, ysize, block_ysize):
            if y + block_ysize < ysize:
                rows = block_ysize
            else:
                rows = ysize - y
            ablai_block = src.ReadAsArray(x, y, cols, rows)
            print(x,y)
            #interpolation
            ablai365 = np.apply_along_axis(lambda ts: interpolate_10day_to_daily(ts, year), 0, ablai_block)
            for c in range(365):
                dst_ds.GetRasterBand(c + 1).WriteArray(ablai365[c, :, :], x, y)
    dst_ds = None

def merge_abLAI(input_2017, input_2018, output_tif):
    src_2018 = gdal.Open(input_2018, gdalconst.GA_ReadOnly)
    src_2017 = gdal.Open(input_2017, gdalconst.GA_ReadOnly)
    proj = src_2017.GetProjection()
    geotrans = src_2017.GetGeoTransform()
    
    Nx, Ny = 2500, 4000
    drv = gdal.GetDriverByName("GTIFF")
    dst_ds = drv.Create(output_tif, Ny, Nx, 730, gdal.GDT_Float32, options=["INTERLEAVE=BAND", "TILED=YES", "BIGTIFF=YES", "COMPRESS=LZW"])
    dst_ds.SetGeoTransform(geotrans)
    dst_ds.SetProjection(proj)

    xsize, ysize = src_2017.RasterXSize, src_2017.RasterYSize
    block_xsize, block_ysize = 1000, 1000
    for x in range(0, xsize, block_xsize):
        if x + block_xsize < xsize:
            cols = block_xsize
        else:
            cols = xsize - x
        for y in range(0, ysize, block_ysize):
            if y + block_ysize < ysize:
                rows = block_ysize
            else:
                rows = ysize - y
            print(x,y)
            lai_2017 = src_2017.ReadAsArray(x, y, cols, rows)
            lai_2018 = src_2018.ReadAsArray(x, y, cols, rows)
            laiall = np.append(lai_2017, lai_2018, axis=0)
            for c in range(730):
                dst_ds.GetRasterBand(c + 1).WriteArray(laiall[c, :, :], x, y)
    dst_ds = None


def cut_abLAI_with_phenology(abLAI_tif, sos_tif, eos_tif, output_tif):
    src = gdal.Open(abLAI_tif, gdalconst.GA_ReadOnly)
    sos_src = gdal.Open(sos_tif, gdalconst.GA_ReadOnly)
    eos_src = gdal.Open(eos_tif, gdalconst.GA_ReadOnly)
    proj = src.GetProjection()
    geotrans = src.GetGeoTransform()
    
    Nx, Ny = 2500, 4000
    drv = gdal.GetDriverByName("GTIFF")
    dst_ds = drv.Create(output_tif, Ny, Nx, 730, gdal.GDT_Float32, options=["INTERLEAVE=BAND", "TILED=YES", "BIGTIFF=YES", "COMPRESS=LZW"])
    dst_ds.SetGeoTransform(geotrans)
    dst_ds.SetProjection(proj)

    xsize, ysize = src.RasterXSize, src.RasterYSize
    block_xsize, block_ysize = 1000, 1000
    for x in range(0, xsize, block_xsize):
        if x + block_xsize < xsize:
            cols = block_xsize
        else:
            cols = xsize - x
        for y in range(0, ysize, block_ysize):
            if y + block_ysize < ysize:
                rows = block_ysize
            else:
                rows = ysize - y
            ablai = src.ReadAsArray(x, y, cols, rows)
            sos = sos_src.ReadAsArray(x, y, cols, rows)
            eos = eos_src.ReadAsArray(x, y, cols, rows)
            print(x, y)

            cut_lai = np.full(ablai.shape, np.nan)
            Nx, Ny = np.shape(sos)
            for ilat in np.arange(Nx):
                for ilon in np.arange(Ny):
                    print(ilat, ilon)
                    lai = pd.Series(ablai[:, ilat, ilon])
                    if np.all(np.isnan(lai)):
                        cut_lai[:, ilat, ilon] = np.nan
                    else:
                        try:
                            cut = lai[int(sos[ilat, ilon]):int(eos[ilat, ilon])]
                            df = pd.DataFrame({'s1': lai, 's2': cut})
                            cut = df['s2']
                            cut_lai[:, ilat, ilon] = cut
                        except:
                            cut_lai[:, ilat, ilon] = np.nan
            for c in range(730):
                dst_ds.GetRasterBand(c + 1).WriteArray(cut_lai[c, :, :], x, y)
    dst_ds = None






