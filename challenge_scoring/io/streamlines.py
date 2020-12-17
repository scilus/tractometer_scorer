#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from dipy.io.streamline import load_tractogram, save_tractogram
from dipy.io.stateful_tractogram import (
    Origin, Space, StatefulTractogram)
from dipy.io.vtk import load_vtk_streamlines


def format_needs_orientation(tract_fname):
    _, extension = os.path.splitext(tract_fname)
    return extension in ['.vtk', '.fib']


def get_tractogram_in_voxel_space(
    tract_fname,
    ref_anat_fname,
    tracts_attribs={'orientation': 'unknown'},
    origin=Origin.NIFTI
):
    if tracts_attribs['orientation'] != 'unknown':
        to_lps = tracts_attribs['orientation'] == 'LPS'
        streamlines = load_vtk_streamlines(tract_fname, to_lps)
        sft = StatefulTractogram(streamlines, ref_anat_fname, Space.RASMM)
    else:
        sft = load_tractogram(
            tract_fname, ref_anat_fname,
            bbox_valid_check=False, trk_header_check=False)
    sft.to_vox()
    sft.to_origin(origin)
    return sft


def save_tracts_from_voxel_space(
    tract_fname, ref_anat_fname, tracts,
    data_per_streamline=None, data_per_point=None
):
    sft = StatefulTractogram(
        tracts, ref_anat_fname, Space.VOX,
        data_per_streamline=data_per_streamline, data_per_point=data_per_point)
    save_tractogram(sft, tract_fname, bbox_valid_check=False)


def save_valid_connections(extracted_vb_info, tractogram,
                           segmented_out_dir, basename, ref_anat_fname,
                           out_tract_type, save_vbs=False, save_full_vc=False):

    if not save_vbs and not save_full_vc:
        return

    full_vcs_idx = []

    for bundle_name, bundle_info in extracted_vb_info.items():
        if bundle_info['nb_streamlines'] > 0:
            out_fname = os.path.join(segmented_out_dir, basename +
                                     '_VB_{}.{}'.format(
                                         bundle_name, out_tract_type))

            idx = list(bundle_info['streamlines_indices'])

            if save_full_vc:
                full_vcs_idx.extend(idx)

            if save_vbs:
                vc_strl = tractogram.streamlines[idx]
                vc_dps = tractogram.data_per_streamline[idx]
                vc_dpp = tractogram.data_per_point[idx]
                save_tracts_from_voxel_space(
                    out_fname, ref_anat_fname,
                    vc_strl, vc_dps, vc_dpp)

    if save_full_vc and len(full_vcs_idx):
        out_name = os.path.join(
            segmented_out_dir, basename + '_VC.{}'.format(out_tract_type))
        full_vcs = tractogram.streamlines[full_vcs_idx]
        full_vc_dps = tractogram.data_per_streamline[full_vcs_idx]
        full_vc_dpp = tractogram.data_per_point[full_vcs_idx]
        save_tracts_from_voxel_space(out_name, ref_anat_fname,
                                     full_vcs, full_vc_dps, full_vc_dpp)


def save_invalid_connections(ib_info, id_invalids, tractogram, ic_clusters,
                             out_segmented_dir, base_name,
                             ref_anat_fname, out_tract_type,
                             save_full_ic=False, save_ibs=False):

    # ib_info is a dictionary containing all the pairs of ROIs that were
    # assigned to some IB. The value of each element is a list containing the
    # clusters indices of clusters that were assigned to that ROI pair.
    if not save_full_ic and not save_ibs:
        return

    full_ic_idx = []
    for k, v in ib_info.items():
        idx = []
        for c_idx in v:
            idx.extend(ic_clusters[c_idx].indices)

        # Get idx from invalid streamline ids
        invalid_idx = [id_invalids[i] for i in idx]

        # Get actual invalid streamlines
        out_strl = tractogram.streamlines[invalid_idx]
        out_dps = tractogram.data_per_streamline[invalid_idx]
        out_dpp = tractogram.data_per_point[invalid_idx]

        if save_ibs:
            out_fname = os.path.join(out_segmented_dir,
                                     base_name +
                                     '_IB_{0}_{1}.{2}'.format(
                                         k[0], k[1], out_tract_type))

            save_tracts_from_voxel_space(out_fname, ref_anat_fname, out_strl,
                                         out_dps, out_dpp)

        if save_full_ic:
            full_ic_idx.extend(invalid_idx)

    if save_full_ic and len(full_ic_idx):
        out_name = os.path.join(out_segmented_dir,
                                base_name + '_IC.{}'.format(out_tract_type))
        full_ic_strl = tractogram.streamlines[full_ic_idx]
        full_ic_dps = tractogram.data_per_streamline[full_ic_idx]
        full_ic_dpp = tractogram.data_per_point[full_ic_idx]
        save_tracts_from_voxel_space(out_name, ref_anat_fname, full_ic_strl,
                                     full_ic_dps, full_ic_dpp)
