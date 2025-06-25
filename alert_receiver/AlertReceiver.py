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
import matplotlib.pyplot as plt
import matplotlib as mpl
from IPython.core.interactiveshell import InteractiveShell
InteractiveShell.ast_node_interactivity = "last_expr"

mpl.rcParams["axes.titlesize"] = 14
mpl.rcParams["axes.labelsize"] = 20
plt.rcParams['savefig.dpi'] = 500
plt.rc('font', family='serif')

import os
#	Append the path to the GECKO project root directory
import sys
sys.path.append('..')


#%

def create_mock_eventlog(numeric_keys):
    """
    Create a mock event log table with predefined schema and mock data.
    
    Args:
        numeric_keys (list): List of numeric column names
    
    Returns:
        astropy.table.Table: Mock event log table
    """
    # Define required columns
    required_columns = ['superevent_id', 'alert_type', 'most_probable_event'] + numeric_keys + ['output_dir', 'processed']
    
    # Create new table with complete schema
    eventlogtbl = Table(
        names=required_columns,
        dtype=['S20', 'S20', 'S20'] + ['f8'] * len(numeric_keys) + ['S50', 'B']
    )
    
    # Add mock data
    mock_data = [
        'MS19941026GP',        # superevent_id
        'MOCK',                # alert_type
        'GPAEK',               # most_probable_event
        -99.0,                 # ramax
        -99.0,                 # decmax
        -99.0,                 # area_90
        -99.0,                 # distmean
        -99.0,                 # diststd
        'MS19941026GP_MOCK_0', # output_dir
        True                   # processed
    ]
    eventlogtbl.add_row(mock_data)
    print("Create mock data to event log table")
    
    return eventlogtbl

#%

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
                    skymap.meta.setdefault('DISTMEAN', TEMP_DIST_MEAN)  # Default mean distance
                    skymap.meta.setdefault('DISTSTD', TEMP_DIST_STD)   # Default standard deviation
                    pass
                elif (record['event']['group'] == 'Burst'):
                    # Set default distance info if not available
                    skymap.meta.setdefault('DISTMEAN', TEMP_DIST_MEAN)  # Default mean distance
                    skymap.meta.setdefault('DISTSTD', TEMP_DIST_STD)   # Default standard deviation

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

def safe_get(rec, key, default=np.nan, cast=float):
    """
    Safely get a value from a dictionary with casting.
    
    Args:
        rec (dict): Input dictionary
        key (str): Key to retrieve
        default: Default value if key not found
        cast (type): Type to cast the value to
        
    Returns:
        value: Casted value or default if failed
    """
    try:
        return cast(rec.get(key, default))
    except Exception:
        return default

#%%
print("=====             Start AlertReceiver.py           =====\n")
print("===== GECKO prowls the cosmic savannah for GW prey =====\n")
print("""                                        
                       .. .             
                      .@%=@=.             |
                     .=*@@@%=           --+--
                     .+=#@@@%  .......    |
                        =@@= .*%@@@@@%*.
     .-=#+@@*.   .+@@@@@@@@@+*@@@@@@@@@@
    =*#@@@@@@#. -@@@@@@@@@@@@@@@@@@@@@@=
    .+@@@**@@@++@@@@@@@@@@@@@@@@@@@@%+. 
    .#--@+=@@@@@@@@@@@@@@@@@@- +%#-.    
       .. .%@@@@@@#=*%@%+@@%.           
     ..:+@@@@@@@@@@#-.  =@@#%..         
   .*@@@@@@@%+:#@@@@@=  -%@@@@#.        
 .#@@@#+-..     +=#@+     -@@@*         
:@@*.          .=#@@#+.                 
%@-            .*%%%-:.                 
@%                :*.                   
%@.                                     
:%@:.                                   
 .+@@@@%:                               
   .....                                
                                        """)

print(f"--------------------------------------------------------\n")

#%%
import yaml

with open('../conf/config.yaml', 'r', encoding='utf-8') as f:
    config = yaml.safe_load(f)


# API Configuration
api_config = config['API_CONFIG']
# GCN Kafka Configuration
GCN_KAFKA_CONFIG = api_config['GCN_KAFKA']
client_id = GCN_KAFKA_CONFIG['client_id']
client_secret = GCN_KAFKA_CONFIG['client_secret']

# Slack Configuration
SLACK_API_CONFIG = api_config['SLACK']
slack_token = SLACK_API_CONFIG['OAuth_Token']
slack_channel = SLACK_API_CONFIG['channel']

# GW Configuration
gw_config = config['GW_SETTINGS']
TEMP_DIST_MEAN = gw_config['DISTANCE']['mean']
TEMP_DIST_STD = gw_config['DISTANCE']['std']

# Observation Settings
obs_config = config['OBSERVATION_SETTINGS']
# Path Configuration
path_out = obs_config['PATHS']['OUTPUT']
path_eventlog = obs_config['PATHS']['EVENT_LOG']

# Data Paths
PATH_TO_GALAXY_CATALOG = obs_config['PATHS']['CATALOG']
PATH_TO_FACILITY = obs_config['PATHS']['FACILITY']

# Tiling Limits
NUMBER_GALAXY_LIMIT = obs_config['LIMITS']['GALAXIES']
NUMBER_TILING_LIMIT = obs_config['LIMITS']['TILES']

# Telescope Settings
telescopes4tile = obs_config['TELESCOPES']['enabled']
EXPTIME_LIMIT_SMNET = obs_config['LIMITS']['EXPOSURE']['SMNET']
EXPTIME_LIMIT_KMTNET = obs_config['LIMITS']['EXPOSURE']['KMTNET']

# Time Settings
time_config = config['TIME_SETTINGS']
TIME_TO_SLEEP = time_config['SLEEP']

# Visualization Settings
vis_config = config['VISUALIZATION']
plot_config = vis_config['PLOT']
plt.rcParams['savefig.dpi'] = int(plot_config['DPI'])
mpl.rcParams["axes.titlesize"] = int(plot_config['TITLE_SIZE'])
mpl.rcParams["axes.labelsize"] = int(plot_config['LABEL_SIZE'])

# Debugging Mode
debug_mode = config['DEBUG_MODE']

print(f"GCN client_id    : {bool(client_id)}")
print(f"Slack OAuth token: {bool(slack_token)}")
print(f"Output path      : {path_out}")
print(f"Event log path   : {path_eventlog}")

# Check and initialize event log table
numeric_keys = ['ramax', 'decmax', 'area_90', 'distmean', 'diststd']

# Ensure output directory exists
os.makedirs(os.path.dirname(path_eventlog), exist_ok=True)

# Try reading existing event log, otherwise create a new one
try:
    eventlogtbl = Table.read(path_eventlog, format='ascii.fixed_width')
    # Check if all required columns exist
    required_columns = ['superevent_id', 'alert_type', 'most_probable_event'] + numeric_keys + ['output_dir', 'processed']
    missing_columns = [col for col in required_columns if col not in eventlogtbl.colnames]
    if missing_columns:
        print(f"Warning: Missing columns in event log: {missing_columns}")
        eventlogtbl = create_mock_eventlog(numeric_keys)
except (FileNotFoundError, OSError):
    eventlogtbl = create_mock_eventlog(numeric_keys)

# Write initial table to ensure it exists
if not os.path.exists(path_eventlog):
    eventlogtbl.write(path_eventlog, format='ascii.fixed_width', overwrite=True)

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
    'broker.address.family': 'v4',
	}

consumer = Consumer(
	config=config,
	client_id=GCN_KAFKA_CONFIG['client_id'],
	client_secret=GCN_KAFKA_CONFIG['client_secret'],
	)
consumer.subscribe(['igwn.gwalert'])



import time
import logging


if __name__ == "__main__":
    while True:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}] IDLE (wait {TIME_TO_SLEEP} sec)", end='\r')
        try:
            for message in consumer.consume(timeout=30*24*60*60):
                record, skymap = AlertReceiver(message.value())
                #	Only Superevent will be distributed to Slack
                if "S" == record['superevent_id'][0:1]:
                    slack = True
                elif "MS" == record['superevent_id'][0:2]:
                    slack = False
                else:
                    slack = True
                    print("I don't know what this is")

                #   New directory
                nn = 0 # Unique number for overlapped name
                path_output = f"{path_out}/{record['superevent_id']}_{record['alert_type']}_{nn}"
                while os.path.exists(path_output):
                    nn += 1
                    path_output = f"{path_out}/{record['superevent_id']}_{record['alert_type']}_{nn}"

                #   Generate directory
                print(f"Generating directory: {path_output}")
                os.makedirs(path_output)
                print(f"Dumping record.json")
                record_str = json.dumps(record)
                #   JSON string to file
                with open(f'{path_output}/record.json', 'w') as file:
                    file.write(record_str)
                
                #   Write skymap to fits
                # Handle superevent entries
                if record['superevent_id'].startswith('S'):
                    # Prepare common variables
                    superevent_id = record.get('superevent_id', '')
                    alert_type    = record.get('alert_type', '')

                    # Determine most_probable_event based on alert_type and event group
                    if alert_type == 'RETRACTION':
                        most_probable_event = None
                    else:
                        group = record.get('event', {}).get('group', '')
                        if group == 'CBC':
                            most_probable_event = max(
                                record['event']['classification'],
                                key=record['event']['classification'].get,
                                default=None
                            )
                        else:
                            most_probable_event = 'Burst'

                    # Safely extract numeric values
                    values = [safe_get(record, k) for k in numeric_keys]

                    # Construct output directory name
                    # output_dir = f"{superevent_id}_{alert_type}"
                    output_dir = os.path.basename(path_output)

                    # Add a single row with all prepared values
                    eventlogtbl.add_row([superevent_id, alert_type, most_probable_event] + values + [output_dir, False])

                    # Format numeric columns for readability
                    for k in numeric_keys:
                        eventlogtbl[k].format = ".3f"

                    # Write the updated table back to the file
                    eventlogtbl.write(path_eventlog, format='ascii.fixed_width', overwrite=True)
                            
            time.sleep(TIME_TO_SLEEP)
            
        except Exception as e:
            logging.error(f"Error while consuming messages: {e}")
            print("Encountered an error while trying to consume Kafka messages.")
            time.sleep(TIME_TO_SLEEP)
