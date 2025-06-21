# GeckoDigestor

A comprehensive package for processing gravitational wave alerts and generating observation plans. GeckoDigestor receives gravitational wave alerts, processes sky localization data, crossmatches with galaxy catalogs, and generates pointing plans for follow-up observations.

## Features

- Real-time processing of gravitational wave alerts from GCN/Kafka
- Sky localization analysis and visualization
- Galaxy catalog crossmatching
- Generation of all possible pointing plans
- Event version tracking for duplicate alerts
- Comprehensive logging and monitoring
- Slack integration for notifications
- Configurable output formats and directories

## Installation

### 1. Using Conda (Recommended)

```bash
# Create a new conda environment
conda create -n geckodigestor python=3.8
conda activate geckodigestor

# Install required packages
conda install -c conda-forge numpy>=1.21.0 pandas>=1.3.0 astropy>=5.0.0 ligo.skymap>=0.6.0
conda install -c conda-forge pydantic>=2.0.0 slack_sdk>=3.20.0 matplotlib>=3.4.0

# Install the package in development mode
pip install -e .
```

### 2. Using pip (Alternative)

```bash
# Create a virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install required packages
pip install numpy>=1.21.0 pandas>=1.3.0 astropy>=5.0.0 ligo.skymap>=0.6.0
pip install pydantic>=2.0.0 slack_sdk>=3.20.0 matplotlib>=3.4.0

# Install the package in development mode
pip install -e .
```

## Configuration

The package uses a JSON configuration file located at `config/config.json`. The configuration includes:

```json
{
    "paths": {
        "data_dir": "data",
        "output_dir": "output",
        "log_dir": "logs",
        "catalog_dir": "catalogs",
        "skymap_dir": "skymaps"
    },
    "galaxy_catalog": {
        "file": "GLADE_2.4.txt",
        "version": "2.4",
        "columns": {
            "ra": "RAJ2000",
            "dec": "DEJ2000",
            "z": "z"
        }
    },
    "slack": {
        "enabled": true,
        "channel": "#gw-alerts",
        "token": "${SLACK_API_TOKEN}"
    },
    "logging": {
        "level": "INFO",
        "max_size": 10485760,
        "backup_count": 5
    }
}
```

## Usage

### 1. Running the Alert Receiver

The AlertReceiver component listens for gravitational wave alerts from the GCN/Kafka topic and processes them:

```bash
# Run the AlertReceiver
python -m geckodigestor.observers.alert_receiver

# Example with custom configuration
python -m geckodigestor.observers.alert_receiver --config /path/to/config.json
```

The AlertReceiver will:
1. Connect to the Kafka topic `gcn.classic.voevent`
2. Listen for gravitational wave alerts
3. Process each alert in real-time
4. Generate logs in the configured log directory
5. Send Slack notifications if enabled

### 2. Processing Events Manually

You can also process events manually using the GeckoDigestor:

```python
from geckodigestor.core.digestor import GeckoDigestor
from pathlib import Path

# Initialize the digestor
config_path = Path("config/config.json")
digestor = GeckoDigestor(config_path)

# Process a specific event directory
event_dir = Path("/path/to/event/S230528ay_INITIAL_0")
digestor.process_event(event_dir)
```

### 3. Running with Docker (Optional)

If you have Docker installed, you can run GeckoDigestor in a container:

```bash
# Build the Docker image
docker build -t geckodigestor .

# Run the container
docker run -v /path/to/config:/app/config \
          -v /path/to/data:/app/data \
          -v /path/to/output:/app/output \
          geckodigestor
```

### 4. Running with Environment Variables

You can override configuration values using environment variables:

```bash
# Set Slack token
export SLACK_API_TOKEN="your_token_here"

# Set custom log level
export GECKODIGESTOR_LOG_LEVEL="DEBUG"

# Run the AlertReceiver
python -m geckodigestor.observers.alert_receiver
```

### 5. Running in Development Mode

For development and testing:

```bash
# Install in development mode
pip install -e .

# Run with debug logging
python -m geckodigestor.observers.alert_receiver --debug

# Run with test configuration
python -m geckodigestor.observers.alert_receiver --config config/test_config.json
```

## Directory Structure

```
geckodigestor/
├── config/              # Configuration files
│   └── config.json     # Main configuration
├── data/               # Input data directory
├── output/             # Processed event outputs
│   ├── S230528ay_INITIAL_0/  # Event directories
│   │   ├── pointings/
│   │   │   └── all_pointings.txt
│   │   ├── sky_localization.png
│   │   └── pointing_summary.png
│   └── ...
├── logs/               # Log files
│   ├── alert_receiver.log
│   ├── geckodigestor.log
│   ├── event_S230528ay_INITIAL_0.log
│   └── ...
└── catalogs/          # Galaxy catalogs
    └── GLADE_2.4.txt
```

## Event Versioning

The system automatically handles duplicate events by:

1. Creating unique event identifiers with version numbers
2. Storing each version in separate directories
3. Maintaining separate logs for each version
4. Including version information in Slack notifications

Example event versions:
- `S230528ay_INITIAL_0` (first INITIAL alert)
- `S230528ay_INITIAL_1` (second INITIAL alert)
- `S230528ay_UPDATE_0` (first UPDATE alert)

## Logging

The system maintains multiple levels of logging:

1. General logs in `logs/geckodigestor.log`
2. Alert receiver logs in `logs/alert_receiver.log`
3. Event-specific logs in `logs/event_{event_id}.log`

Log levels:
- INFO: General processing information
- DEBUG: Detailed data and configuration
- WARNING: Non-critical issues
- ERROR: Processing errors with full context

## Error Handling

The system includes comprehensive error handling with:

1. Detailed error messages
2. Full traceback logging
3. Event-specific error tracking
4. Slack notifications for critical errors

## Monitoring

The system provides monitoring through:

1. Real-time logging
2. Slack notifications
3. Event progress tracking
4. Performance metrics

## Contributing

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## License

This project is licensed under the MIT License - see the LICENSE file for details.

# Author
1. Gregory S.H. Paek (백승학)
2. Hyunho Choi (최현호)

## Version log
- 2023.09.01: version 0.1

# INDEX
- What is GECKO?
- TBD
- ...

# 1. What is GECKO?
Gratitational-wave Electromagnetic-wave Counterpart in Korea Observatory (GECKO; [M. Im et al. (2019)](http://yokohamagrb2019.wikidot.com/proceedings)) project is aiming to find kilonovae (KNe), the optical/NIR counterpart of GW, with the network of 1-2m class telescopes in the world.

## 1.1. Facilties
Both projects share the same facilities. They consist of more than 10 telescopes described in below:

|Facility|Location|Description|
|:---:|:---:|:---:|
|SAO        |Korea       |TBD|
|LOAO       |USA         |TBD|
|LSGT       |Austrailia  |TBD|
|KCT        |Austrailia  |TBD|
|RASA36     |Chile       |TBD|
|CBNUO      |Korea       |TBD|
|KMTNet_SSO |Austrailia  |TBD|
|KMTNet_SAAO|South Africa|TBD|
|KMTNet_CTIO|Chile       |TBD|

# 2. Requirements
We recommand the following requirements about the version:
- This pipeline requires the version of Python == 3.11.3

**We highly recommend to use the state-of-the-art `Python >= 3.11` to maximize the computing speed.**

## 2.1. Python Library
- `numpy == 1.23.5`
- `scipy == 1.10.1`
- `matplotlib == 3.7.1`
- `astropy == 5.1`
- `astroscrappy == 1.1.0`
- `ccdproc == 2.4.0`
- `requests == 2.28.1`: To utilize slack API
- *`alipy`
<!-- - `PyRAF >= X` -->

## 2.2. External Software
- TBD

# 3. Installation
```
$ git clone https://github.com/SilverRon/gppy
```

# 4. Structure and Usage
`GeckoDigestor` consists of # parts.

## 4.1. `AlertReciever.py`
- ...

# 5. Features 
## 5.1. Structure
...

## 5.2. Output
```
TBD
```

# 6. Future Update
- TBD

# 7. License and Copyright
TBD

# 8. Contact
- Gregory S.H. Paek (gregorypaek94___at___gmail.com) @Seoul National University
- Hyunho Choi (...___at___gmail.com) @Seoul National University

# 9. Reference
1. GECKO

# 10. Acknowledgments
<!-- Thanks to our GECKO team members, Prof. Myungshin Im, Dr. Changsu Choi, Dr. Seo-won Chang, Dr. Gu Lim, and Dr. Sophia Kim.
Especially, special thanks to Dr. Changsu Choi, who made techincal foundations in the beggining of the IMSNG project, inspired and motivated me to develop this pipeline. -->
TBD
