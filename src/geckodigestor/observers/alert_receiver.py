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
import logging
import os
import sys
from typing import Any, Dict, Optional
import json
import traceback
from astropy_healpix import healpy as ah
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from ..core.digestor import GeckoDigestor
from ..config.config_manager import config
from ..utils.logger import logger
import tempfile
from base64 import b64decode
from io import BytesIO
from gcn_kafka import Consumer
import numpy as np
from astropy.io import fits
from ligo.skymap.io import write_sky_map

#%%
def AlertReceiver(record):
    """
    Process a gravitational wave alert record
    
    Args:
        record: JSON string containing the alert data
        
    Returns:
        tuple: (processed_record, skymap)
        processed_record: Processed alert record with additional metadata
        skymap: Astropy Table containing the skymap data
    """
    try:
        logger.info("Starting alert processing")
        record = json.loads(record)
        
        # Process alert based on type
        if record['alert_type'] == 'RETRACTION':
            logger.info(f"Retraction received for {record['superevent_id']}")
            return record, None

        # Check if this is a real event (starts with 'S')
        if not record['superevent_id'].startswith('S'):
            logger.warning(f"Skipping non-real event {record['superevent_id']}")
            return record, None

        # Log event details
        logger.info(f"Processing event {record['superevent_id']} - {record['alert_type']}")
        logger.debug(f"Event details: {json.dumps(record, indent=2)}")

        # Process skymap if present
        if record is not None and 'event' in record and 'skymap' in record['event']:
            try:
                event_id = record['superevent_id']
                alert_type = record['alert_type']
                event_version = record['event_version']
                
                logger.info(f"Processing skymap for {event_id} - {alert_type} (version {event_version})")
                skymap_str = record['event']['skymap']
                skymap_bytes = b64decode(skymap_str)
                skymap = Table.read(BytesIO(skymap_bytes))
                logger.debug(f"Skymap loaded successfully with {len(skymap)} pixels")
                
                # Get most probable sky location
                level, ipix = ah.uniq_to_level_ipix(
                    skymap[np.argmax(skymap['PROBDENSITY'])]['UNIQ']
                )
                ra, dec = ah.healpix_to_lonlat(ipix, ah.level_to_nside(level), order='nested')
                logger.info(f"Most probable location: RA={ra.deg:.4f}, Dec={dec.deg:.4f}")
                
                # Calculate 90% confidence area
                skymap.sort('PROBDENSITY', reverse=True)
                prob = skymap['PROBDENSITY'] * ah.nside_to_pixel_area(ah.level_to_nside(ah.uniq_to_level_ipix(skymap['UNIQ'])))
                cumprob = np.cumsum(prob)
                area_90 = np.sum(prob[cumprob <= 0.9])
                logger.info(f"90% confidence area: {area_90:.2f} deg²")
                
                # Add event parameters to record
                record['event']['ra'] = ra.deg
                record['event']['dec'] = dec.deg
                record['event']['area_90'] = area_90
                
                logger.info(f"Skymap processing completed for {event_id} - {alert_type} (version {event_version})")
                
                # Handle event type
                if record['event']['group'] == 'CBC':
                    pass
                elif record['event']['group'] == 'Burst':
                    # Set default distance info for burst events
                    skymap.meta['DISTMEAN'] = 1000
                    skymap.meta['DISTSTD'] = 100

                # Log sky location and distance
                logger.info(f'Most probable sky location (RA, Dec) = ({ra.deg:.3f}, {dec.deg:.3f})')
                logger.info(f'Distance = {skymap.meta["DISTMEAN"]:.1f} +/- {skymap.meta["DISTSTD"]:.1f} Mpc')
                
                return record, skymap
            except Exception as e:
                logger.error(f"Error processing skymap for {event_id} - {alert_type} (version {event_version}): {str(e)}")
                logger.error(f"Error details: {traceback.format_exc()}")
                return record, None

        return record, None

    except Exception as e:
        logger.error(f"Error processing alert: {str(e)}")
        logger.error(f"Error details: {traceback.format_exc()}")
        return record, None

    # JSON decoding error handling
    except json.JSONDecodeError as e:
        logger.error(f"JSON decoding failed: {e}")
        return record, None
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
consumer = Consumer(
    config={
        'group.id': 'geckodigestor',
        'max.poll.interval.ms': 1000000,
        'auto.offset.reset': 'latest'
    },
    client_id='geckodigestor',
    client_secret=None  # Set this in your config if needed
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
