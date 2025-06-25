# GeckoDigestor
GECKO project's code for digesting the alert and managing observation to detect kilonovae faster using gravitational wave alerts and rational methods for ranking observations across GECKO facilities worldwide.

# Contact & Author
Gregory S.H. Paek (백승학)
- Email: gregorypaek94___at___gmail.com
- Affiliation: Institute for Astronomy (IfA)

## Version log
- 2025.06.25: version 1.0
- 2023.09.01: version 0.1

# 1. What is GECKO?
Gratitational-wave Electromagnetic-wave Counterpart in Korea Observatory (GECKO; [M. Im et al. (2019)](http://yokohamagrb2019.wikidot.com/proceedings)) project aims to find kilonovae (KNe), the optical/NIR counterparts of gravitational waves (GW), using a network of 1-2m class telescopes worldwide. GECKO has demonstrated its capability in following up GW events through observations of the binary neutron star-black hole merger candidate S230518h ([Paek et al., 2025](https://iopscience.iop.org/article/10.3847/1538-4357/adaf99)) and the binary neutron star merger GW190425 ([Paek et al., 2024](https://iopscience.iop.org/article/10.3847/1538-4357/ad0238)). These studies showcase GECKO's wide-field follow-up capabilities and its role in constraining properties of KNe from NSBH mergers.

## 1.1. Facilities
Both projects share the same facilities. They consist of more than 10 telescopes described in below:

|Facility|Location|Note|
|:---:|:---:|:---:|
|SAO        |Korea       |-|
|LOAO       |USA         |-|
|LSGT       |Austrailia  |-|
|KCT        |Austrailia  |Tiling|
|RASA36     |Chile       |Tiling|
|CBNUO      |Korea       |Tiling|
|KMTNet |Austrailia, South Africa, Chile  |Tiling (SSO, SAAO, CTIO)|
|UKIRT      |Hawaii      |Tiling (Temp for ToO)|

# 2. Requirements
We recommend the following requirements:
- Python >= 3.11.3

**We highly recommend using the latest version of Python (>= 3.11) to maximize computing speed and compatibility.**

- `numpy >= 1.23.5`
- `scipy >= 1.10.1`
- `matplotlib >= 3.7.1`
- `astropy >= 5.1`
- `astroscrappy >= 1.1.0`
- `ccdproc >= 2.4.0`
- `requests >= 2.28.1`: For Slack API integration
- `ligo.skymap`: For sky map processing
- `gcn_kafka`: For GCN alert reception
- `yaml`: For configuration management
- `healpy`: For HEALPix operations
- `astropy_healpix`: For HEALPix operations
- `matplotlib`: For visualization
- `IPython`: For interactive shell features
- `matplotlib`: For visualization
- `mpl_toolkits`: For advanced plotting

# 3. Installation
```
$ git clone https://github.com/SilverRon/gppy
```

# 4. Structure and Usage
`GeckoDigestor` consists of two main components:

## 4.1. `AlertReceiver.py`
- **Purpose**: Receives and processes gravitational wave alerts
- **Main Features**:
  - Receives alerts through GCN Kafka in real time
  - Processes alerts into a format usable by `GeckoDigestor.py`
- **Usage**: 
  ```bash
  cd alert_receiver
  python AlertReceiver.py
  ```

## 4.2. `GeckoDigestor.py`
- **Purpose**: Analyzes gravitational wave events and generates observation plans
- **Main Features**:
  - Processes gravitational wave sky maps in real time in coperation with the `AlertReceiver.py`
  - Process the skymap locatlization area
  - Select and Prioritize host galaxy candidates targets and tiles based on $\rm P_{3D}$ and $\rm P_{3D} * M_{\odot}$ (stellar mass)
  - Generates sorted tiles and host galaxy candidates
  - Integrates Slack notifications
- **Usage**:
  ```bash
  cd gecko_digestor
  python GeckoDigestor.py
  ```
- **Key Components**:
  - Integration of various astronomical libraries (astropy, healpy)
  - Structured configuration management
  - Observation facility constraint handling
  - Tiling and galaxy limits

Both modules work together in the GECKO project's workflow:
1. `AlertReceiver` receives gravitational wave alerts
2. `GeckoDigestor` analyzes the alerts and generates observation plans
3. The plans are used by GECKO's observation network
  eventlogtbl = Table.read(path_eventlog, format='ascii.fixed_width')
  
- **Key Components**:
  - Integration of various astronomical libraries (astropy, healpy)
  - Structured configuration management
  - Observation facility constraint handling
  - Tiling and galaxy limits

Both modules work together in the GECKO project's workflow:
1. `AlertReceiver` receives gravitational wave alerts
2. `GeckoDigestor` analyzes the alerts and generates observation plans
3. The plans are used by GECKO's observation network


# 5. Outputs
`GeckoDigestor` generates comprehensive outputs for each gravitational wave event, organized in a dedicated directory named after the event ID. The outputs include:

## 5.1. Catalog Files
- `HostGalaxyCatalog_90.csv`: Catalog of potential host galaxies within the 90% confidence region
- `HostGalaxyCatalog_90.reg`: Region file for visualizing host galaxies in astronomical software
- `SkyGridCatalog_[TELESCOPE]_S230518h_90%.csv`: Observation grid catalogs for each telescope (7DT, CBNUO, KCT, KMTNet, RASA36, UKIRT)

## 5.2. Visualization Files
- `skymap.png`: Probability distribution map of the event
- `skymap_no_rot.png`: Unrotated version of the sky map
- `skymap_rot90.png`: Rotated version of the sky map
- `tiling_[TELESCOPE].png`: Tiling pattern for each telescope
- `tiling_[TELESCOPE]_cumulative_dist.png`: Cumulative distribution of tiling points
- `tiling_[TELESCOPE]_zoom_[RA]_[Dec].png`: Zoomed-in view of tiling around the most probable position
- `cumulative_p3d_HostGalaxy.png`: 3D probability distribution of host galaxies
- `lc.png`: Light curve analysis results

## 5.3. Data Files
- `record.json`: Detailed processing record and metadata
- `summary.txt`: Summary of event information including:
  - Event ID and trigger time
  - Classification and distance estimates
  - Sky position coordinates
  - Sky area coverage statistics

Example output directory structure:
```
S230518h_PRELIMINARY/
├── HostGalaxyCatalog_90.*
├── SkyGridCatalog_[TELESCOPE]_S230518h_90%.csv
├── skymap*.png
├── tiling_[TELESCOPE]*.png
├── cumulative_p3d_HostGalaxy.png
├── lc.png
├── record.json
└── summary.txt
```

These outputs provide a complete set of information for follow-up observations of gravitational wave events, including target selection, observation planning, and analysis results. 

# 6. Configuration

`GeckoDigestor` uses a comprehensive configuration file (`config.yaml`) to manage all system settings. The configuration is organized into several sections:

## 6.1. API Configuration
- **GCN Kafka Settings**: Authentication for receiving gravitational wave alerts
- **Slack Integration**: Settings for sending notifications to Slack

## 6.2. Gravitational Wave Settings
- **Distance Settings**: Mean and standard deviation for distance estimates
- **Confidence Intervals**: Thresholds for high (90%) and medium (50%) confidence regions

## 6.3. Observation Settings
- **Paths**: Configuration for data directories (catalogs, facilities, skymaps, output)
- **Limits**: Maximum number of galaxies and tiles
- **Telescope Settings**: Enabled telescopes and exposure limits

## 6.4. Time Settings
- **Processing**: Timeout and check interval settings
- **Sleep**: Waiting time between operations

## 6.5. Visualization Settings
- **Plot Settings**: Graph resolution and text sizes
- **Slack Messages**: Image inclusion and message length limits

## 6.6. Event Processing Settings
- **Timeouts**: Maximum processing time
- **Probability Thresholds**: Criteria for NS and burst signals
- **Distance Limits**: Maximum mean and standard deviation

## 6.7. Debug Mode
- **Enabled/Disabled**: Controls debug logging and output

Example configuration file:
```yaml
# config.yaml
API_CONFIG:
  GCN_KAFKA:
    client_id: "your_client_id"
    client_secret: "your_client_secret"
  SLACK:
    OAuth_Token: "your_slack_token"
    channel: "your_channel"

GW_SETTINGS:
  DISTANCE:
    mean: 100  # Mpc
    std: 10    # Mpc
  CONFIDENCE_INTERVALS:
    HIGH: 0.9
    MEDIUM: 0.5

OBSERVATION_SETTINGS:
  PATHS:
    CATALOG: "../data/GLADE+_240919.fits" # Galaxy catalog (GLADE+)
    FACILITY: "../data/gecko.facilities.tsv"
    OUTPUT: "/data/gpaek/GeckoDigestor/output"

  LIMITS:
    GALAXIES: 500
    TILES: 500

  TELESCOPES:
    enabled: ['7DT', 'KMTNet', 'RASA36', 'KCT', 'CBNUO', 'UKIRT']

TIME_SETTINGS:
  SLEEP: 5  # seconds
  PROCESSING:
    TIMEOUT: 3600
    CHECK_INTERVAL: 60

DEBUG_MODE: true
```

# 7. License and Copyright
TBD

# 7. Contact
- Gregory S.H. Paek (gregorypaek94___at___gmail.com) @Seoul National University

# 8. Acknowledgments
<!-- Thanks to our GECKO team members, Prof. Myungshin Im, Dr. Changsu Choi, Dr. Seo-won Chang, Dr. Gu Lim, and Dr. Sophia Kim.
Especially, special thanks to Dr. Changsu Choi, who made techincal foundations in the beggining of the IMSNG project, inspired and motivated me to develop this pipeline. -->
TBD
