#!/usr/bin/env python

import pyfits
import subprocess
import os
import h5py
import healpy
import requests
import numpy as np
import argparse
import tempfile
import StringIO
import shutil


parser = argparse.ArgumentParser(description='Download FEBECop beams and produce HDF file.')
parser.add_argument('input_file', help='Text file with lines of "FREQLABEL IPIX"')
parser.add_argument('output_prefix', help='Prefix of output HDF file; "FREQLABEL.h5" will be appended')
parser.add_argument('--threshold', default=1e-300, type=float,
                    help='Truncate beam at this level, relative to maximum in beam map')
parser.add_argument('--release', default='PR2', help='Data release label')
parser.add_argument('--nside', default=1024, type=int, help='nside in use by portal...')
args = parser.parse_args()


def process_beam(args, freq_label, ipix):
    theta, phi = healpy.pix2ang(args.nside, ipix)
    theta = theta * 360. / (2. * np.pi) - 90
    phi = phi * 360. / (2. * np.pi)
    
    try:
        tmp_dir = tempfile.mkdtemp()

        # Get JSON data pointing to gzipped FITS file
        r = requests.get('http://pla.esac.esa.int/pla-sl/data-action?'
                        'ROI_LON={phi:.4f}&ROI_LAT={theta:.4f}&FREQUENCY=30&RELEASE={release}&ProductType=BEAM'.format(
                        theta=theta, phi=phi, release=args.release))
        uri = r.json()[0]

        # Check that we agree with portal on ipix..
        print ipix, theta, phi, uri
        ##assert ('_%s.fits.gz' % ipix) in uri

        # Download it to disk
        r = requests.get(uri)
        with open('%s/beam.fits.gz' % tmp_dir, 'w') as f:
            f.write(r.content)
        # Unzip and read map into memory
        subprocess.check_call(['gunzip', '%s/beam.fits.gz' % tmp_dir])
        primary, extension = pyfits.open('%s/beam.fits' % tmp_dir)
        map = extension.data.field(0)
    finally:
        shutil.rmtree(tmp_dir)

    nside = int(np.sqrt(map.shape[0] // 12))
    assert nside == args.nside  # or else user needs to pass --nside
    threshold = args.threshold * map.max()
    indices = (map < threshold).nonzero()[0].astype(np.int32)
    values = map[indices].astype(np.float32)
        
    # Convert to sparse. We re-open the HDF file every time in order to be more
    # robust against crashes (which would leave the file corrupt if opened)
    with h5py.File('%s%s.h5' % (args.output_prefix, freq_label), 'a') as f:
        group_name = '/%d' % ipix
        f.create_dataset(group_name + '/indices', data=indices)
        f.create_dataset(group_name + '/values', data=values)
        grp = f[group_name]
        grp.attrs['nside'] = nside


with open(args.input_file) as f:
    for line in f.readlines():
        cols = [x.strip() for x in line.split()]
        freq_label, ipix = cols[:2]
        ipix = int(ipix)
        process_beam(args, freq_label, ipix)
