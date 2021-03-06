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


FREQLABEL_TO_NSIDE = {
  '30': 1024,
  '43': 1024
}


parser = argparse.ArgumentParser(description='Download FEBECop beams and produce HDF file.')
parser.add_argument('input_file', help='Text file with lines of "FREQLABEL IPIX"')
parser.add_argument('output_prefix', help='Prefix of output HDF file; "FREQLABEL.h5" will be appended')
parser.add_argument('--threshold', default=1e-300, type=float,
                    help='Truncate beam at this level, relative to maximum in beam map')
parser.add_argument('--release', default='PR2', help='Data release label')
args = parser.parse_args()


def process_beam(args, freq_label, ipix):
    nside = FREQLABEL_TO_NSIDE[freq_label]
    theta, phi = healpy.pix2ang(nside, ipix)
    theta = 90 - theta * 180. / np.pi
    phi = phi * 180. / np.pi    
    try:
        tmp_dir = tempfile.mkdtemp()

        # Get JSON data pointing to gzipped FITS file
        r = requests.get('http://pla.esac.esa.int/pla-sl/data-action?'
                        'ROI_LON={phi:.4f}&ROI_LAT={theta:.4f}&FREQUENCY={freq_label}&RELEASE={release}&ProductType=BEAM'.format(
                        theta=theta, phi=phi, release=args.release, freq_label=freq_label))
        uri = r.json()[0]
        print 'Downloading: %s' % uri

        # Check that we agree with portal on ipix..
        assert ('_%s.fits.gz' % ipix) in uri

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

    assert 12 * nside**2 == map.shape[0]  # otherwise we are wrong about data hosted in portal
    threshold = args.threshold * map.max()
    indices = (np.abs(map) > threshold).nonzero()[0].astype(np.int32)
    values = map[indices].astype(np.float32)
            
    # Convert to sparse. We re-open the HDF file every time in order to be more
    # robust against crashes (which would leave the file corrupt if opened)
    with h5py.File('%s%s.h5' % (args.output_prefix, freq_label), 'a') as f:
        group_name = '/%d' % ipix
        f.create_dataset(group_name + '/indices', data=indices)
        f.create_dataset(group_name + '/values', data=values)
        grp = f['/']
        grp.attrs['nside'] = nside


with open(args.input_file) as f:
    lines = f.readlines()
    n = len(lines)
    for i, line in enumerate(lines):
        cols = [x.strip() for x in line.split()]
        if cols == [] or cols[0].startswith('#'):
            continue
        freq_label, ipix = cols[:2]
        ipix = int(ipix)
        print '%.1f %%' % ((i * 100.) / n)
        process_beam(args, freq_label, ipix)
