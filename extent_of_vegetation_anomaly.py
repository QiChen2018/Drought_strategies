# -*- coding: utf-8 -*-
"""
Created on Sun May 25 15:38:54 2025

@author: chenq4
"""

import numpy as np
import gdal
import gdalconst

def get_ablaibsum(index, threshold):
    z, x, y = np.shape(index)
    laibsum_ = np.full((x, y), np.nan)

    for ilat in range(x):
        for ilon in range(y):
            a = index[:, ilat, ilon]
            laiball_ = []
            if np.all(np.isnan(a)):
                laibsum_[ilat, ilon] = np.nan
            elif len(np.where(a < threshold)[0]) == 0:
                laibsum_[ilat, ilon] = 0
            else:
                laibx = np.where(a < threshold)[0]
                for i in laibx:
                    laiball_.append(a[i])
                laibsum_[ilat, ilon] = np.nansum(laiball_)

    return laibsum_

def write_ablaibsum_output(src, dst_ds, threshold):
    xsize = src.RasterXSize
    ysize = src.RasterYSize
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
            abLAI_block = src.ReadAsArray(x, y, cols, rows)
            print(x, y)
            abLAI_block_sum = get_ablaibsum(abLAI_block, threshold)
            dst_ds.GetRasterBand(1).WriteArray(abLAI_block_sum, x, y)

    dst_ds = None

def process_ablaibsum_file(input_path, output_path, threshold):
    src = gdal.Open(input_path, gdalconst.GA_ReadOnly)
    proj = src.GetProjection()
    geotrans = src.GetGeoTransform()

    Nx, Ny = 2500, 4000
    drv = gdal.GetDriverByName("GTIFF")
    dst_ds = drv.Create(output_path, Ny, Nx, 1, gdal.GDT_Float32, options=["INTERLEAVE=BAND", "TILED=YES"])
    dst_ds.SetGeoTransform(geotrans)
    dst_ds.SetProjection(proj)

    write_ablaibsum_output(src, dst_ds, threshold)

