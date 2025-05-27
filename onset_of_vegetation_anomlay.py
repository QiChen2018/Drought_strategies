# -*- coding: utf-8 -*-
"""
Created on Tue May 20 19:59:46 2025

@author: chenq4
"""

import pandas as pd
import numpy as np
import gdal
import gdalconst




def GetMaxevent(LAIevent):
    maks = max(LAIevent, key=lambda k: len(LAIevent[k]))
    return LAIevent[maks]

def _sign_grouper(anomalies):
    sign = np.sign(anomalies)
    sign[sign == 0] = -1
    runs = (sign != sign.shift(1)).astype(int).cumsum()
    runs[anomalies.isnull()] = np.nan
    runs = runs - (runs.min() - 1)
    return anomalies.groupby(runs)

def get_runs(anomalies):
    return {
        name: _sign_grouper(anomalies=anomalies).get_group(name=name)
        for name in _sign_grouper(anomalies=anomalies).indices
    }

def pool_runs(runs, pooling_method='None', show_positives=False, **kwargs):
    if pooling_method == 'None':
        runs_pooled = runs
    if show_positives:
        return runs_pooled
    else:
        return {num: run for num, run in runs_pooled.items() if run.sum() < 0}

def get_longlai(lai):
    longestlai = np.full(730, np.nan)
    if np.all(np.isnan(lai)):
        longestlai[:] = np.nan
    else:
        try:
            x0 = pd.Series([-1] * 730)
            lairuns = get_runs(anomalies=lai - x0)
            LAIevent = pool_runs(
                runs=lairuns,
                pooling_method='None',
            )
            lailong = GetMaxevent(LAIevent) - 1
            df = pd.DataFrame({'s1': lai, 's2': lailong})
            lailong = df['s2']
            longestlai[:] = lailong
        except:
            longestlai[:] = np.nan
    return longestlai

def compute_longest_lai_event(input_tif, output_tif):
    src = gdal.Open(input_tif, gdalconst.GA_ReadOnly)
    proj = src.GetProjection()
    geotrans = src.GetGeoTransform()
    
    Nx, Ny = 2500, 4000
    drv = gdal.GetDriverByName("GTIFF")
    dst_ds = drv.Create(output_tif, Ny, Nx, 730, gdal.GDT_Float32,
                        options=["INTERLEAVE=BAND", "TILED=YES", "BIGTIFF=YES", "COMPRESS=LZW"])
    dst_ds.SetGeoTransform(geotrans)
    dst_ds.SetProjection(proj)

    xsize, ysize = src.RasterXSize, src.RasterYSize
    block_xsize, block_ysize = 4000, 100
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
            print(x, y)
    
            longestlai_ = np.full(ablai.shape, np.nan)
            Nz, Nx, Ny = ablai.shape
            for ilat in np.arange(Nx):
                Inputs = []
                for ilon in np.arange(Ny):
                    print(x, y, ilat, ilon)
                    lai = pd.Series(ablai[:, ilat, ilon])
                    Inputs.append(lai)
                ret = [get_longlai(lai) for lai in Inputs]
                longestlai_[:, ilat, :] = np.transpose(np.array(ret), axes=[1, 0])
            for i in range(730):
                dst_ds.GetRasterBand(i + 1).WriteArray(longestlai_[i, :, :], x, y)
    dst_ds = None
    print('Finished computing longest LAI event.')

def compute_longest_lai_onset(input_tif, output_tif):
    src = gdal.Open(input_tif, gdalconst.GA_ReadOnly)
    proj = src.GetProjection()
    geotrans = src.GetGeoTransform()
    
    Nx, Ny = 2500, 4000
    drv = gdal.GetDriverByName("GTIFF")
    dst_ds = drv.Create(output_tif, Ny, Nx, 1, gdal.GDT_Float32,
                        options=["INTERLEAVE=BAND", "TILED=YES", "COMPRESS=LZW"])
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
            longestlai = src.ReadAsArray(x, y, cols, rows)
            print(x, y)
            
            Nz, Nx, Ny = longestlai.shape
            lonset = np.full((Nx, Ny), np.nan)
            for ilat in np.arange(Nx):
                for ilon in np.arange(Ny):
                    onellai = longestlai[:, ilat, ilon]
                    if np.all(np.isnan(onellai)):
                        lonset[ilat, ilon] = np.nan
                    elif not np.any(onellai):
                        lonset[ilat, ilon] = np.nan
                    else:
                        lonset[ilat, ilon] = np.argwhere(np.isfinite(onellai))[0][0]
            dst_ds.GetRasterBand(1).WriteArray(lonset, x, y)
    dst_ds = None
    print('Finished computing onset')


