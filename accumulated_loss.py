# -*- coding: utf-8 -*-
"""
Created on Sun May 25 16:09:47 2025

@author: chenq4
"""

import numpy as np
import gdal
import gdalconst
from itertools import accumulate

# Function to calculate accumulation
def get_accum(yearablai):
    if np.all(np.isnan(yearablai)):
        ablai_accum = np.zeros((730))
        ablai_accum = ablai_accum * np.NaN
    else:
        yearablai[np.isnan(yearablai)] = 0.0
        yearablai[np.where(yearablai > -1)] = 0.0
        ablai_accum = list(accumulate(yearablai))
    return ablai_accum

# Function to process abLAI in blocks and write output
def write_accumulated_output(src, dst_ds):
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
            ablai_accum_730 = np.apply_along_axis(get_accum, 0, abLAI_block)
            for c in range(730):
                print(c + 1)
                dst_ds.GetRasterBand(c + 1).WriteArray(ablai_accum_730[c, :, :], x, y)

    dst_ds = None

# Main driver
def process_abLAI_accumulation(input_path, output_path):
    src = gdal.Open(input_path, gdalconst.GA_ReadOnly)
    proj = src.GetProjection()
    geotrans = src.GetGeoTransform()

    Nx, Ny = 2500, 4000
    drv = gdal.GetDriverByName("GTIFF")
    dst_ds = drv.Create(output_path, Ny, Nx, 730, gdal.GDT_Float32,
                        options=["INTERLEAVE=BAND", "TILED=YES", "BIGTIFF=YES", "COMPRESS=LZW"])
    dst_ds.SetGeoTransform(geotrans)
    dst_ds.SetProjection(proj)

    write_accumulated_output(src, dst_ds)


