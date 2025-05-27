# -*- coding: utf-8 -*-
"""
Created on Sun May 25 16:22:13 2025

@author: chenq4
"""
import sys
import os
try:
    base_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    base_dir = os.getcwd()
sys.path.append(base_dir)

from calculate_abLAI import (calculate_abLAI, interpolate_abLAI, merge_abLAI, cut_abLAI_with_phenology)
from onset_of_vegetation_anomlay import (compute_longest_lai_event, compute_longest_lai_onset)
from extent_of_vegetation_anomaly import process_ablaibsum_file
from accumulated_loss import process_abLAI_accumulation



def main():
    data_dir = os.path.join(base_dir, "data")
    output_dir = os.path.join(base_dir, "results")
    os.makedirs(output_dir, exist_ok=True)
    
    # calculate_abLAI
    lai_input = os.path.join(data_dir, "LAI_0418_001_eu.tiff")
    sos = os.path.join(data_dir, "phenology_sos_001.tiff")
    eos = os.path.join(data_dir, "phenology_eos_001.tiff")
    
    ablai_2018 = os.path.join(output_dir, "abLAIblock_EU_001_10day.tiff")
    ablai_2017 = os.path.join(output_dir, "abLAIblock_EU_001_10day_2017.tiff")
    ablai_2018_365 = os.path.join(output_dir, "abLAI365block_EU_001.tiff")
    ablai_2017_365 = os.path.join(output_dir, "abLAI365block_EU_001_2017.tiff")
    merged = os.path.join(output_dir, "abLAIblock_EU_001_365_1718.tiff")
    final_cut = os.path.join(output_dir, "abLAIblock_EU_001_730_cutsoseos.tiff")


    calculate_abLAI(lai_input, ablai_2018, target_year_index=14, baseline_indices=list(range(14)))
    calculate_abLAI(lai_input, ablai_2017, target_year_index=13, baseline_indices=list(range(14)))
    interpolate_abLAI(ablai_2018, ablai_2018_365, 2018)
    interpolate_abLAI(ablai_2017, ablai_2017_365, 2017)
    merge_abLAI(ablai_2017_365, ablai_2018_365, merged)
    cut_abLAI_with_phenology(merged, sos, eos, final_cut)

    # onset
    longest_event_output = os.path.join(output_dir, "longestLAIevent_001_730_phenology.tiff")
    onset_output = os.path.join(output_dir, "longestLAI_onset_001_phenology.tiff")

    compute_longest_lai_event(input_tif=final_cut, output_tif=longest_event_output)
    compute_longest_lai_onset(input_tif=longest_event_output, output_tif=onset_output)

    # extent
    extent_output = os.path.join(output_dir, "abLAIb1sum_EU_001_cutphenology.tiff")
    process_ablaibsum_file(input_path=final_cut, output_path=extent_output, threshold=-1)

    # accumulation
    accum_output = os.path.join(output_dir, "abLAIb1_EU_001_cutphenology_accum.tiff")
    process_abLAI_accumulation(input_path=final_cut, output_path=accum_output)

    print("All processing completed")


if __name__ == "__main__":
    main()
