"""
This module receives and processes gravitational wave alerts.

This module provides the AlertReceiver class which receives alerts
and processes them into a format that can be used by other modules in the
GECKO project.

Example:
    To receive and process alerts, create an instance of AlertReceiver
    class and call its start() method.

    alert_receiver = AlertReceiver()
    alert_receiver.start()

Attributes:
    ALERT_SERVER (str): The URL of the server from which to receive alerts.
"""
#%%
from base64 import b64decode
from io import BytesIO
import json
from pprint import pprint

import tempfile
from astropy.table import Table
import astropy_healpix as ah
from astropy import units as u
from gcn_kafka import Consumer
import numpy as np
from astropy.io import fits
from ligo.skymap.io import write_sky_map

import os
#	Append the path to the GECKO project root directory
import sys
sys.path.append('..')
from config.config import GCN_KAFKA_CONFIG
#%%
def AlertReceiver(record):
    record = json.loads(record)

    # Only respond to mock events. Real events have GraceDB IDs like
    # S1234567, mock events have GraceDB IDs like M1234567.
    # NOTE NOTE NOTE replace the conditional below with this commented out
    # conditional to only parse real events.
    # if record['superevent_id'][0] != 'S':
    #    return
    if record['superevent_id'][0] != 'M':
        return

    if record['alert_type'] == 'RETRACTION':
        print(record['superevent_id'], 'was retracted')
        return

    # Respond only to 'CBC' events. Change 'CBC' to 'Burst' to respond to
    # only unmodeled burst events.
    if record['event']['group'] != 'CBC':
        return
    print(f"="*60)
    print(f"EVENT: {record['superevent_id']}-{record['alert_type']}")
    print(f"="*60)
    # Parse sky map
    # skymap_str = record.get('event', {}).pop('skymap')
    skymap_str = record['event']['skymap']
    if skymap_str:
        # Decode, parse skymap, and print most probable sky location
        skymap_bytes = b64decode(skymap_str)
        skymap = Table.read(BytesIO(skymap_bytes))

        level, ipix = ah.uniq_to_level_ipix(
            skymap[np.argmax(skymap['PROBDENSITY'])]['UNIQ']
        )
        ra, dec = ah.healpix_to_lonlat(ipix, ah.level_to_nside(level),
                                       order='nested')
        print(f'Most probable sky location (RA, Dec) = ({ra.deg:.3f}, {dec.deg:.3f})')

        # Print some information from FITS header
        print(f'Distance = {skymap.meta["DISTMEAN"]:.1f} +/- {skymap.meta["DISTSTD"]:.1f} Mpc')

        #   Calculate the area of the skymap
        skymap.sort('PROBDENSITY', reverse=True)
        level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
        pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
        prob = pixel_area * skymap['PROBDENSITY']
        cumprob = np.cumsum(prob)
        i = cumprob.searchsorted(0.9)
        area_90 = pixel_area[:i].sum()

        #   Put additional information to the record
        record['ramax'] = ra.deg
        record['decmax'] = dec.deg
        record['area_90'] = area_90.to_value(u.deg**2)
        record['distmean'] = skymap.meta['DISTMEAN']
        record['diststd'] = skymap.meta['DISTSTD']

    print(f"-"*60)
    # Print remaining fields
    print('Record:')
    pprint(record)


    # 새로 생성할 디렉토리 이름
    new_directory = f"../output/{record['superevent_id']}_{record['alert_type']}"
    # 디렉토리가 존재하지 않으면 생성
    if not os.path.exists(new_directory):
        os.makedirs(new_directory)
    # 딕셔너리를 JSON 문자열로 변환
    record_str = json.dumps(record)
    # JSON 문자열을 파일로 저장
    with open(f'{new_directory}/record.json', 'w') as file:
        file.write(record_str)

    write_skymap_to_fits(skymap, path_output=f"{new_directory}/skymap.fits")


    return record, skymap

def read_skymap_bytes_to_table(skymap_str):
    skymap_bytes = b64decode(skymap_str)
    skymap = Table.read(BytesIO(skymap_bytes))
    return skymap

def write_skymap_to_fits(skymap, path_output):
    # 아래 예시에 대한 부분; skymap array의 경우도 원하는 컬럼만 활용할 수도 있음.
    with tempfile.NamedTemporaryFile(suffix='.fits') as f:
        write_sky_map(f.name, skymap, nest=True)
        for card in fits.getheader(f.name, 1).cards:
            print(str(card).rstrip())

    # 실제 파일로 저장하는 부분
    write_sky_map(path_output, skymap, nest=True)
#%%
#	Online Test
# consumer = Consumer(client_id=GCN_KAFKA_CONFIG['client_id'], client_secret=GCN_KAFKA_CONFIG['client_secret'],)
# consumer.subscribe(['igwn.gwalert'])

# while True:
#     for message in consumer.consume():
#         record, skymap = AlertReceiver(message.value())
#%%
#	Offline Test
# with open('../data/MS181101ab-preliminary.json', 'r') as f:
#     record = f.read()

# AlertReceiver(record)