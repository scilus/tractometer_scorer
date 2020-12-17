#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import json
import os


def load_attribs(attribs_filepath):
    with open(attribs_filepath, 'rb') as attribs_file:
        attribs = json.load(attribs_file)
    return attribs


def get_attribs_for_file(attribs_filepath, filename):
    attribs = load_attribs(attribs_filepath)
    return attribs.get(os.path.basename(filename))


def save_attribs(attribs_filepath, attribs):
    with open(attribs_filepath, 'wb') as attribs_file:
        json.dump(attribs, attribs_file)


def merge_attribs(orig_attribs, additional_attribs, overwrite=False):
    for fname, attr in additional_attribs.iteritems():
        orig_value = orig_attribs.get(fname, None)

        if orig_value is None:
            # Item was not in dictionary
            orig_attribs[fname] = attr
        else:
            for new_attr, new_val in attr.iteritems():
                orig_attrib_val = orig_value.get(new_attr, None)

                # If we didn't find the attribute, or the attribute
                # was found and overwrite is set.
                if orig_attrib_val is None or overwrite:
                    orig_value[new_attr] = new_val
                else:
                    raise ValueError('Couldn\'t overwrite, flag is false.')

    return orig_attribs


# Considers that all files are children of the root dir.
# Attrib_func must accept a filename as a param
def compute_attrib_files(root_dir, attrib_func, attrib_name):
    attribs = {}

    fnames = glob.glob(os.path.join(root_dir, '*'))

    for fname in fnames:
        val = attrib_func(fname)
        attribs[os.path.basename(fname)] = {attrib_name: val}

    return attribs
