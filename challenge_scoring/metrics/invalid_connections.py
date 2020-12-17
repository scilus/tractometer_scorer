#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import division

from collections import Counter
import logging
import os
import random

import numpy as np

from dipy.segment.clustering import QuickBundles
from scipy.spatial.distance import cdist

from challenge_scoring.io.streamlines import save_invalid_connections
from challenge_scoring.utils.filenames import get_root_image_name


def find_closest_distance_points_to_region(points, roi_volume):
    roi_coords = roi_volume
    dists = cdist(points, roi_coords, 'euclidean')

    return np.min(dists, axis=1)


def find_closest_region(points, rois):
    # points is a 2D array, each row being 1 point.
    # We expect 2 rows (start and end points)
    # rois is a list of (region name, region_data)

    min_global_dists = [100000, 100000]
    closest_region_names = ["", ""]

    for region_name, region_data in rois:
        closest_dists = find_closest_distance_points_to_region(points,
                                                               region_data)

        if closest_dists[0] < min_global_dists[0]:
            min_global_dists[0] = closest_dists[0]
            closest_region_names[0] = region_name
        if closest_dists[1] < min_global_dists[1]:
            min_global_dists[1] = closest_dists[1]
            closest_region_names[1] = region_name

    return (closest_region_names, min_global_dists)


def get_closest_roi_pairs_for_bundle(streamlines, rois):
    # Needs to be 2D for cdist
    start_point = np.reshape(streamlines[0][0], (-1, 3))

    closest_rois_pairs = []

    for s in streamlines:
        endpoints = np.vstack([s[0], s[-1]])
        endpoints_dists = cdist(start_point, endpoints).flatten()

        # Make sure we all start from the same "orientation" for streamlines,
        # to try to get the same region as the first region
        if endpoints_dists[0] > endpoints_dists[1]:
            endpoints = np.vstack([s[-1], s[0]])

        closest_region_names, min_dists = find_closest_region(endpoints, rois)
        closest_rois_pairs.append(tuple(closest_region_names))

    occurences = Counter(closest_rois_pairs)

    return occurences.most_common(1)[0][0]


def get_closest_roi_pairs_for_all_streamlines(streamlines, rois):
    """
    Find the closest pair of ROIs from the endpoints of each provided
    streamline.

    :param streamlines: list of streamlines to assign rois to
    :param rois: list of pairs of roi "names" and data
    :return: list of pairs of the closest regions for each bundle head and tail
    """

    # Needs to be 2D for cdist
    start_point = np.reshape(streamlines[0][0], (-1, 3))

    closest_rois_pairs = []

    for s in streamlines:
        endpoints = np.vstack([s[0], s[-1]])
        endpoints_dists = cdist(start_point, endpoints).flatten()

        # Make sure we all start from the same "orientation" for streamlines,
        # to try to get the same region as the first region
        if endpoints_dists[0] > endpoints_dists[1]:
            endpoints = np.vstack([s[-1], s[0]])

        closest_region_names, min_dists = find_closest_region(endpoints, rois)
        closest_rois_pairs.append(tuple(closest_region_names))

    return closest_rois_pairs


def group_and_assign_ibs(tractogram, candidate_ids, ROIs,
                         save_ibs, save_full_ic,
                         out_segmented_dir, base_name,
                         ref_anat_fname, out_tract_type):
    ic_counts = 0
    ib_pairs = {}

    rejected_indices = []

    # Start by clustering all the remaining potentiel IC using QB.

    # Fix seed to always generate the same output
    # Shuffle to try to reduce the ordering dependency for QB
    random.seed(0.2)
    random.shuffle(candidate_ids)

    # TODO threshold on distance as arg for other datasets
    qb = QuickBundles(threshold=20., metric='MDF_12points')
    candidate_streamlines = tractogram.streamlines[candidate_ids]
    clusters = qb.cluster(candidate_streamlines)

    logging.debug("Found {} potential IB clusters".format(len(clusters)))

    # Prefetch information about the bundles endpoints regions of interest.
    # Is used in the get_closest_roi_pairs... function.
    rois_info = []
    for roi in ROIs:
        rois_info.append((
            get_root_image_name(os.path.basename(roi.get_filename())),
            np.array(np.where(roi.get_data())).T))

    all_ics_closest_pairs = get_closest_roi_pairs_for_all_streamlines(
        candidate_streamlines, rois_info)

    for c_idx, c in enumerate(clusters):
        closest_for_cluster = [all_ics_closest_pairs[i]
                               for i in c.indices]

        # Clusters containing only a single streamlines are rejected.
        if len(c.indices) > 1:
            ic_counts += len(c.indices)
            occurences = Counter(closest_for_cluster)

            # TODO could be changed in future to allow an equality
            most_frequent = occurences.most_common(1)[0][0]

            val = ib_pairs.get(most_frequent)
            if val is None:
                # Check if flipped pair exists
                val1 = ib_pairs.get((most_frequent[1], most_frequent[0]))
                if val1 is not None:
                    val1.append(c_idx)
                else:
                    ib_pairs[most_frequent] = [c_idx]
            else:
                val.append(c_idx)
        else:
            rejected_indices.append(c.indices[0])

    if save_ibs or save_full_ic:
        save_invalid_connections(ib_pairs, candidate_ids, tractogram,
                                 clusters, out_segmented_dir,
                                 base_name,
                                 ref_anat_fname,
                                 out_tract_type,
                                 save_full_ic=save_full_ic,
                                 save_ibs=save_ibs)

    return rejected_indices, ic_counts, len(ib_pairs.keys())
