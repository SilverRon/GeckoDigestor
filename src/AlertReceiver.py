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
from config.config import *
#%%
def AlertReceiver(record):
    record = json.loads(record)

    # Only respond to mock events. Real events have GraceDB IDs like
    # S1234567, mock events have GraceDB IDs like M1234567.
    # NOTE NOTE NOTE replace the conditional below with this commented out
    # conditional to only parse real events.
    # if record['superevent_id'][0] != 'S':
    #    return
    # if record['superevent_id'][0] != 'M':
    #     return record, None

    if record['alert_type'] == 'RETRACTION':
        print(record['superevent_id'], 'was retracted')
        return record, None

    # Respond only to 'CBC' events. Change 'CBC' to 'Burst' to respond to
    # only unmodeled burst events.
    # if record['event']['group'] != 'CBC':
    #     return
    print(f"="*60)
    print(f"EVENT: {record['superevent_id']}-{record['alert_type']}")
    print(f"="*60)
    # Parse sky map
    # skymap_str = record.get('event', {}).pop('skymap')
    # if (skymap_str) & (skymap_str != None):
    # if (record['event']['skymap'] != None):
    # if record is not None and 'event' in record and 'skymap' in record['event']:
    # if record is not None and 'event' in record:


    #   Burst-like event
    # In the context of gravitational waves,
    # a signal candidate that is detected without a template and without prior knowledge of the waveform.
    # Examples of potential sources of gravitational-wave bursts include
    # 1. high mass BBH mergers, 2. core-collapse supernovae, and 3. cosmic string cusps.

    if record is not None and 'event' in record and record['alert_type'] != 'RETRACTION':
        try:
            if ('skymap' in record['event']):
                #   24.05.01
                """
                skymap_str = record['event']['skymap']
                # Decode, parse skymap, and print most probable sky location
                skymap_bytes = b64decode(skymap_str)
                skymap = Table.read(BytesIO(skymap_bytes))

                level, ipix = ah.uniq_to_level_ipix(
                    skymap[np.argmax(skymap['PROBDENSITY'])]['UNIQ']
                )
                ra, dec = ah.healpix_to_lonlat(ipix, ah.level_to_nside(level),
                                            order='nested')

                #   Calculate the area of the skymap
                skymap.sort('PROBDENSITY', reverse=True)
                level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
                pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
                prob = pixel_area * skymap['PROBDENSITY']
                cumprob = np.cumsum(prob)
                i = cumprob.searchsorted(0.9)
                area_90 = pixel_area[:i].sum()
                #
                if (record['event']['group'] == 'CBC'):
                    pass
                elif (record['event']['group'] == 'Burst'):
                    #   Put some distance info
                    skymap.meta['DISTMEAN'] = 1000
                    skymap.meta['DISTSTD'] = 100
                #
                print(f'Most probable sky location (RA, Dec) = ({ra.deg:.3f}, {dec.deg:.3f})')
                print(f'Distance = {skymap.meta["DISTMEAN"]:.1f} +/- {skymap.meta["DISTSTD"]:.1f} Mpc')
                """

                skymap_str = record['event']['skymap']
                # Decode, parse skymap, and print most probable sky location
                skymap_bytes = b64decode(skymap_str)
                skymap = Table.read(BytesIO(skymap_bytes))

                level, ipix = ah.uniq_to_level_ipix(
                    skymap[np.argmax(skymap['PROBDENSITY'])]['UNIQ']
                )
                ra, dec = ah.healpix_to_lonlat(ipix, ah.level_to_nside(level),
                                            order='nested')

                #   Calculate the area of the skymap
                skymap.sort('PROBDENSITY', reverse=True)
                level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
                pixel_area = ah.nside_to_pixel_area(ah.level_to_nside(level))
                prob = pixel_area * skymap['PROBDENSITY']
                cumprob = np.cumsum(prob)
                i = cumprob.searchsorted(0.9)
                area_90 = pixel_area[:i].sum()
                #
                if (record['event']['group'] == 'CBC'):
                    skymap.meta.setdefault('DISTMEAN', 1000)  # Default mean distance
                    skymap.meta.setdefault('DISTSTD', 100)   # Default standard deviation
                    pass
                elif (record['event']['group'] == 'Burst'):
                    # Set default distance info if not available
                    skymap.meta.setdefault('DISTMEAN', 1000)  # Default mean distance
                    skymap.meta.setdefault('DISTSTD', 100)   # Default standard deviation

                # Output the most probable sky location and distance information safely
                print(f'Most probable sky location (RA, Dec) = ({ra.deg:.3f}, {dec.deg:.3f})')
                distmean = skymap.meta.get('DISTMEAN', 'Unknown')
                diststd = skymap.meta.get('DISTSTD', 'Unknown')
                if distmean == 'Unknown' or diststd == 'Unknown':
                    print('Distance information is not available.')
                else:
                    print(f'Distance = {distmean:.1f} +/- {diststd:.1f} Mpc')


                if record['event']['group'] == 'Burst':
                    print(f"*Burst event: Distance is brutly set to run the following code")
                print(f"Area = {area_90.to(u.deg**2).value:.1f} deg2")

                #   Put additional information to the record
                record['ramax'] = ra.deg
                record['decmax'] = dec.deg
                record['area_90'] = area_90.to_value(u.deg**2)
                record['distmean'] = skymap.meta['DISTMEAN']
                record['diststd'] = skymap.meta['DISTSTD']
        except json.JSONDecodeError as e:
            print(f"[ERROR] JSON decoding failed: {e}")
            return 


        print(f"-"*60)
        # Print remaining fields
        # print('Record:')
        # pprint(record)


        # # 새로 생성할 디렉토리 이름
        # path_output = f"../output/{record['superevent_id']}_{record['alert_type']}"
        # # 디렉토리가 존재하지 않으면 생성
        # if not os.path.exists(path_output):
        #     os.makedirs(path_output)
        # # 딕셔너리를 JSON 문자열로 변환
        # record_str = json.dumps(record)
        # # JSON 문자열을 파일로 저장
        # with open(f'{path_output}/record.json', 'w') as file:
        #     file.write(record_str)

        # write_skymap_to_fits(skymap, path_output=f"{path_output}/skymap.fits")


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
#	Offline Test
# with open('../data/MS181101ab-preliminary.json', 'r') as f:
#     record = f.read()
# record, skymap = AlertReceiver(record)

#   Burst signal
# with open('../data/S230528ay-preliminary.json,1', 'r') as f:
#     record = f.read()
# record, skymap = AlertReceiver(record)


 #%%
#	Online Test
config = {
	'group.id': '',
	'max.poll.interval.ms': 1000000,
	# 'auto.offset.reset': 'earliest',
    # 'auto.offset.reset': 'latest',
	}

consumer = Consumer(
	config=config,
	client_id=GCN_KAFKA_CONFIG['client_id'],
	client_secret=GCN_KAFKA_CONFIG['client_secret'],
	)
consumer.subscribe(['igwn.gwalert'])

print("Waiting for a GW alert...")
# while False:
#   Before 240904
"""if __name__ == "__main__":
    while True:
        for message in consumer.consume(timeout=30*24*60*60):
            record, skymap = AlertReceiver(message.value())
            #	Superevent?
            if "S" == record['superevent_id'][0:1]:
                print("This is a superevent")
                slack = True
            elif "MS" == record['superevent_id'][0:2]:
                slack = False
                print("This is not a superevent")
            else:
                slack = True
                print("I don't know what this is")

            # 새로 생성할 디렉토리 이름
            path_output = f"{path_out}/{record['superevent_id']}_{record['alert_type']}"
            # 디렉토리가 존재하지 않으면 생성
            if not os.path.exists(path_output):
                os.makedirs(path_output)
            # 딕셔너리를 JSON 문자열로 변환
            record_str = json.dumps(record)
            # JSON 문자열을 파일로 저장
            with open(f'{path_output}/record.json', 'w') as file:
                file.write(record_str)
            if "S" == record['superevent_id'][0:1]:
                eventlogtbl = Table.read(f"{path_out}/event.log", format='ascii.fixed_width')
                if record['alert_type'] == 'RETRACTION':
                    most_probable_event = "None"
                    eventlogtbl.add_row(
                        [record['superevent_id'], record['alert_type'], most_probable_event, 0.0, 0.0, 0.0, 0.0, 0.0, f"{record['superevent_id']}_{record['alert_type']}"]
                    )

                else:
                    if record['event']['group'] == 'CBC':
                        most_probable_event = max(record['event']['classification'], key=record['event']['classification'].get)
                    else:
                        most_probable_event = 'Burst'

                    try:
                        eventlogtbl.add_row(
                            [record['superevent_id'], record['alert_type'], most_probable_event, record['ramax'], record['decmax'], record['area_90'], record['distmean'], record['diststd'], f"{record['superevent_id']}_{record['alert_type']}"]
                        )
                    except:
                        eventlogtbl.add_row(
                            [record['superevent_id'], record['alert_type'], most_probable_event, 0.0, 0.0, 0.0, 0.0, 0.0, f"{record['superevent_id']}_{record['alert_type']}"]
                        )

                    for key in eventlogtbl.keys():
                        if key in ['ramax', 'decmax', 'area_90', 'distmean', 'diststd']:
                            eventlogtbl[key].format = ".3f"
                        
                # eventlogtbl.write(f"{path_out}/event.log", format='ascii.fixed_width', overwrite=True)"""

import time
import logging



if __name__ == "__main__":
    while True:
        try:
            for message in consumer.consume(timeout=30*24*60*60):
                record, skymap = AlertReceiver(message.value())
                #	Superevent?
                if "S" == record['superevent_id'][0:1]:
                    print("This is a superevent")
                    slack = True
                elif "MS" == record['superevent_id'][0:2]:
                    slack = False
                    print("This is not a superevent")
                else:
                    slack = True
                    print("I don't know what this is")

                # 새로 생성할 디렉토리 이름
                path_output = f"{path_out}/{record['superevent_id']}_{record['alert_type']}"
                # 디렉토리가 존재하지 않으면 생성
                if not os.path.exists(path_output):
                    os.makedirs(path_output)
                # 딕셔너리를 JSON 문자열로 변환
                record_str = json.dumps(record)
                # JSON 문자열을 파일로 저장
                with open(f'{path_output}/record.json', 'w') as file:
                    file.write(record_str)
                if "S" == record['superevent_id'][0:1]:
                    eventlogtbl = Table.read(f"{path_out}/event.log", format='ascii.fixed_width')
                    if record['alert_type'] == 'RETRACTION':
                        most_probable_event = "None"
                        eventlogtbl.add_row(
                            [record['superevent_id'], record['alert_type'], most_probable_event, 0.0, 0.0, 0.0, 0.0, 0.0, f"{record['superevent_id']}_{record['alert_type']}"]
                        )

                    else:
                        if record['event']['group'] == 'CBC':
                            most_probable_event = max(record['event']['classification'], key=record['event']['classification'].get)
                        else:
                            most_probable_event = 'Burst'

                        try:
                            eventlogtbl.add_row(
                                [record['superevent_id'], record['alert_type'], most_probable_event, record['ramax'], record['decmax'], record['area_90'], record['distmean'], record['diststd'], f"{record['superevent_id']}_{record['alert_type']}"]
                            )
                        except:
                            eventlogtbl.add_row(
                                [record['superevent_id'], record['alert_type'], most_probable_event, 0.0, 0.0, 0.0, 0.0, 0.0, f"{record['superevent_id']}_{record['alert_type']}"]
                            )

                        for key in eventlogtbl.keys():
                            if key in ['ramax', 'decmax', 'area_90', 'distmean', 'diststd']:
                                eventlogtbl[key].format = ".3f"
                            
                    # eventlogtbl.write(f"{path_out}/event.log", format='ascii.fixed_width', overwrite=True)"""
        except Exception as e:
            logging.error(f"Error while consuming messages: {e}")
            print("Encountered an error while trying to consume Kafka messages.")
            # 일정 시간 대기 후 재시도
            time.sleep(10)  # 10초 대기
