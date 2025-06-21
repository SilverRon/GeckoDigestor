"""
Core module for processing gravitational wave alerts and generating observation plans
"""

import json
import logging
from typing import Optional, Dict, Any
from pathlib import Path
import traceback

import astropy_healpix as ah
import numpy as np
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time
from ligo.skymap.io import read_sky_map
from ligo.skymap.postprocess import crossmatch

from ..config.config_manager import config
from ..utils.logger import logger
from ..utils.slack_utils import send_slack_notification
from astropy.visualization import make_lupton_rgb
import matplotlib.pyplot as plt

# Configure logging
log_dir = Path(config.log_file).parent
log_dir.mkdir(parents=True, exist_ok=True)

def setup_event_logger(event_id: str) -> logging.Logger:
    """Set up event-specific logger"""
    event_logger = logging.getLogger(f'geckodigestor.event.{event_id}')
    event_logger.setLevel(logging.DEBUG)
    
    # Create event-specific file handler
    event_fh = logging.FileHandler(log_dir / f'event_{event_id}.log')
    event_fh.setLevel(logging.DEBUG)
    
    # Create formatter
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    event_fh.setFormatter(formatter)
    
    event_logger.addHandler(event_fh)
    return event_logger

class GeckoDigestor:
    """Main class for processing gravitational wave alerts and generating observation plans"""
    
    def __init__(self):
        """Initialize the GeckoDigestor"""
        self.logger = logger
        self.glade_catalog = self._load_glade_catalog()
        self.slack_config = config.slack
        
        # Create required directories
        for path in config.paths.values():
            Path(path).mkdir(parents=True, exist_ok=True)
        
        # Log initialization
        self.logger.info("GeckoDigestor initialized successfully")
        self.logger.debug(f"Using catalog: {self.glade_catalog.meta['catalog_name']} v{self.glade_catalog.meta['catalog_version']}")
        self.logger.debug(f"Slack notifications: {'enabled' if self.slack_config and self.slack_config['enabled'] else 'disabled'}")
        
    def _create_sky_localization_plot(self, skymap: Table, event_id: str) -> str:
        """Create sky localization plot and save as PNG
        
        Args:
            skymap: Skymap data
            event_id: Event ID for filename
            
        Returns:
            Path to the saved PNG file
        """
        try:
            # Create plot
            fig, ax = plt.subplots(figsize=(10, 8))
            
            # Plot skymap
            level, ipix = ah.uniq_to_level_ipix(skymap['UNIQ'])
            ra, dec = ah.healpix_to_lonlat(ipix, ah.level_to_nside(level), order='nested')
            
            # Plot probability density
            ax.scatter(ra.deg, dec.deg, c=skymap['PROBDENSITY'],
                      cmap='viridis', s=100)
            
            # Add grid and labels
            ax.grid(True, linestyle='--', alpha=0.7)
            ax.set_xlabel('RA (deg)')
            ax.set_ylabel('Dec (deg)')
            ax.set_title(f'Sky Localization for {event_id}')
            
            # Save plot
            output_dir = Path(config.paths.output_dir) / event_id
            output_dir.mkdir(parents=True, exist_ok=True)
            
            plot_path = str(output_dir / 'sky_localization.png')
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            
            return plot_path
            
        except Exception as e:
            self.logger.error(f"Error creating sky localization plot: {str(e)}")
            self.logger.error(f"Error details: {traceback.format_exc()}")
            raise

    def _create_detailed_plots(self, event_id: str, crossmatch_results: Table) -> list:
        """Create detailed analysis plots
        
        Args:
            event_id: Event ID
            crossmatch_results: Table of crossmatch results
            
        Returns:
            List of paths to created plot files
        """
        try:
            plots = []
            output_dir = Path(config.paths.output_dir) / event_id
            
            # Create probability vs. number of galaxies plot
            fig, ax = plt.subplots(figsize=(10, 8))
            ax.plot(crossmatch_results['confidence'],
                   crossmatch_results['number_of_galaxies'],
                   marker='o')
            ax.set_xlabel('Confidence')
            ax.set_ylabel('Number of Galaxies')
            ax.set_title(f'Galaxy Confidence Distribution - {event_id}')
            
            plot_path = str(output_dir / 'galaxy_confidence.png')
            plt.savefig(plot_path, dpi=300)
            plt.close()
            plots.append(plot_path)
            
            return plots
            
        except Exception as e:
            self.logger.error(f"Error creating detailed plots: {str(e)}")
            return []

    def _send_slack_notification(self, event_id: str, skymap_path: str, 
                               crossmatch_results: Table) -> None:
        """Send Slack notification with event information
        
        Args:
            event_id: Event ID
            skymap_path: Path to sky localization plot
            crossmatch_results: Table of crossmatch results
        """
        try:
            # Create sky localization plot
            sky_plot_path = self._create_sky_localization_plot(skymap_path, event_id)
            
            # Send initial thread with basic information and sky plot
            thread_result = send_slack_notification(
                channel=self.slack_config.channel,
                text=f"New GW Event: {event_id}\n"
                    f"Sky localization area processed successfully",
                file_path=sky_plot_path
            )
            
            # Create detailed plots
            detailed_plots = self._create_detailed_plots(event_id, crossmatch_results)
            
            # Send detailed information as reply
            detailed_text = (
                f"Detailed Analysis Results:\n"
                f"- Number of candidate galaxies: {len(crossmatch_results)}\n"
                f"- Confidence levels: {', '.join(map(str, config.observations.confidence_limits.values()))}\n"
                f"- Facilities: {', '.join(f['name'] for f in config.facilities)}"
            )
            
            send_slack_notification(
                channel=self.slack_config.channel,
                text=detailed_text,
                file_path=detailed_plots[0] if detailed_plots else None,
                comment_text=None,
                thread_ts=thread_result['thread_ts']
            )
            
        except Exception as e:
            self.logger.error(f"Error sending Slack notification: {str(e)}")
            raise

    def process_event(self, event_dir: Path) -> None:
        """Process a single gravitational wave event
        
        Args:
            event_dir: Directory containing event data
        """
        try:
            self.logger.info(f"Processing event in {event_dir}")
            
            # Load event data
            record_path = event_dir / 'record.json'
            skymap_path = event_dir / 'skymap.fits'
            
            if not record_path.exists() or not skymap_path.exists():
                self.logger.error("Missing required files")
                return

            with open(record_path) as f:
                record = json.load(f)
            
            # Process skymap
            skymap = read_sky_map(str(skymap_path), moc=True)
            
            # Crossmatch with GLADE catalog
            crossmatch_results = self._crossmatch_glade(skymap)
            
            # Generate observation plans
            self._generate_observation_plans(record, crossmatch_results)
            
            # Send Slack notification
            self._send_slack_notification(record['superevent_id'], skymap, crossmatch_results)
            
        except Exception as e:
            self.logger.error(f"Error processing event: {str(e)}")
            raise

    def _load_glade_catalog(self) -> Table:
        """Load and prepare the galaxy catalog"""
        try:
            # Get catalog configuration
            catalog_config = config.galaxy_catalog
            catalog_path = Path(config.catalog_dir) / catalog_config['file']
            
            # Load catalog
            cat = Table.read(str(catalog_path))
            
            # Create SkyCoord column using configured column names
            cat['coordinates'] = SkyCoord(
                cat[catalog_config['columns']['ra']]*u.deg,
                cat[catalog_config['columns']['dec']]*u.deg
            )
            
            # Add catalog metadata
            cat.meta['catalog_name'] = catalog_config['name']
            cat.meta['catalog_version'] = catalog_config['version']
            
            return cat
        except Exception as e:
            self.logger.error(f"Error loading galaxy catalog: {str(e)}")
            raise

    def _load_facilities(self) -> Table:
        """Load and prepare the facilities catalog"""
        try:
            facilities_path = Path(config.data_dir) / 'gecko.facilities.tsv'
            gcktbl = ascii.read(str(facilities_path))
            
            # Calculate depth for different exposure times
            exposure_times = [3, 4, 10, 24, 30, 60]  # in minutes
            for exp_time in exposure_times:
                col_name = f'depth_{exp_time}min'
                gcktbl[col_name] = self._scale_depth(
                    gcktbl['depth'], 
                    gcktbl['exptime'], 
                    exp_time
                )
            
            return gcktbl
        except Exception as e:
            self.logger.error(f"Error loading facilities: {str(e)}")
            raise

    @staticmethod
    def _scale_depth(depth: float, exptime: float, target_exptime: float) -> float:
        """Scale depth based on exposure time"""
        return depth * np.sqrt(target_exptime / exptime)

    def process_event(self, event_dir: Path) -> None:
        """Process a single gravitational wave event
        
        Args:
            event_dir: Directory containing event data
        """
        try:
            self.logger.info(f"Starting event processing for {event_dir}")
            
            # Load event data
            record_path = event_dir / 'record.json'
            skymap_path = event_dir / 'skymap.fits'
            
            if not record_path.exists() or not skymap_path.exists():
                self.logger.error(f"Missing required files in {event_dir}")
                return

            with open(record_path) as f:
                record = json.load(f)
            
            # Extract event information
            event_id = record['superevent_id']
            alert_type = record['alert_type']
            event_version = record.get('event_version', 0)
            
            # Set up event-specific logger
            event_key = f"{event_id}_{alert_type}_{event_version}"
            event_logger = setup_event_logger(event_key)
            event_logger.info(f"Starting processing for event {event_id} - {alert_type} (version {event_version})")
            event_logger.debug(f"Event details: {json.dumps(record, indent=2)}")
            
            # Process skymap
            event_logger.info("Reading skymap")
            skymap = read_sky_map(str(skymap_path), moc=True)
            event_logger.info(f"Skymap loaded successfully with {len(skymap)} pixels")
            
            # Crossmatch with GLADE catalog
            event_logger.info("Performing galaxy crossmatch")
            crossmatch_results = self._crossmatch_glade(skymap)
            event_logger.info(f"Found {len(crossmatch_results)} potential targets")
            
            # Generate all possible pointings
            self._generate_pointings(record, crossmatch_results)
            
        except Exception as e:
            error_msg = f"Error processing event {record.get('superevent_id', 'unknown')}: {str(e)}"
            self.logger.error(error_msg)
            self.logger.error(f"Error details: {traceback.format_exc()}")
            if 'event_logger' in locals():
                event_logger.error(error_msg)
                event_logger.error(f"Error details: {traceback.format_exc()}")
            raise

    def _crossmatch_glade(self, skymap: Table) -> Table:
        """Crossmatch GLADE catalog with skymap
        
        Args:
            skymap: Probability density map
            
        Returns:
            Table of crossmatch results
        """
        try:
            # Get most probable sky location
            level, ipix = ah.uniq_to_level_ipix(
                skymap[np.argmax(skymap['PROBDENSITY'])]['UNIQ']
            )
            ra, dec = ah.healpix_to_lonlat(ipix, ah.level_to_nside(level), order='nested')
            
            # Crossmatch with GLADE catalog
            crossmatch_results = crossmatch(
                self.glade_catalog['coordinates'],
                ra.deg, dec.deg,
                radius=0.1 * u.deg
            )
            
            return crossmatch_results
        except Exception as e:
            self.logger.error(f"Error in crossmatch: {str(e)}")
            raise

    def _generate_pointings(self, record: Dict[str, Any], results: Table) -> None:
        """Generate all possible pointings based on crossmatch results
        
        Args:
            record: Event record
            results: Crossmatch results
        """
        try:
            event_id = record['superevent_id']
            alert_type = record['alert_type']
            event_version = record.get('event_version', 0)
            event_key = f"{event_id}_{alert_type}_{event_version}"
            
            event_logger = logging.getLogger(f'geckodigestor.event.{event_key}')
            
            # Create output directory
            event_dir = Path(config.paths.output_dir) / event_key
            pointing_dir = event_dir / 'pointings'
            pointing_dir.mkdir(parents=True, exist_ok=True)
            
            # Save pointings with metadata
            pointing_file = pointing_dir / 'all_pointings.txt'
            results.meta['catalog_name'] = self.glade_catalog.meta['catalog_name']
            results.meta['catalog_version'] = self.glade_catalog.meta['catalog_version']
            results.meta['event_id'] = event_id
            results.meta['ra'] = record['event']['ra']
            results.meta['dec'] = record['event']['dec']
            results.meta['area_90'] = record['event']['area_90']
            
            event_logger.info(f"Saving pointings to {pointing_file}")
            results.write(pointing_file, format='ascii', overwrite=True)
            event_logger.info(f"Saved {len(results)} pointings")
            
            # Create pointing summary plot
            event_logger.info("Creating pointing summary plot")
            plot_path = self._create_pointing_summary_plot(event_id, results)
            event_logger.info(f"Plot saved to {plot_path}")
            
            # Send Slack notification with plots
            if self.slack_config and self.slack_config['enabled']:
                event_logger.info("Sending Slack notification")
                send_slack_notification(
                    channel=self.slack_config['channel'],
                    message=f"New event processed: {event_id} - {alert_type} (version {event_version})\n"
                           f"Number of pointings: {len(results)}\n"
                           f"Catalog used: {self.glade_catalog.meta['catalog_name']} v{self.glade_catalog.meta['catalog_version']}\n"
                           f"90% confidence area: {record['event']['area_90']:.2f} deg²",
                    files=[plot_path]
                )
            
            event_logger.info("Pointing generation completed successfully")
            
        except Exception as e:
            error_msg = f"Error generating pointings for {record.get('superevent_id', 'unknown')}: {str(e)}"
            self.logger.error(error_msg)
            self.logger.error(f"Error details: {traceback.format_exc()}")
            if 'event_logger' in locals():
                event_logger.error(error_msg)
                event_logger.error(f"Error details: {traceback.format_exc()}")
            raise

    def _create_pointing_summary_plot(self, event_id: str, results: Table) -> str:
        """Create pointing summary plot showing all possible pointings
        
        Args:
            event_id: Event ID
            results: Crossmatch results
            
        Returns:
            Path to the saved plot file
        """
        try:
            fig, ax = plt.subplots(figsize=(12, 8))
            
            # Plot all pointings
            ax.scatter(results['ra'], results['dec'],
                      label='Pointings', alpha=0.7)
            
            ax.set_xlabel('RA (deg)')
            ax.set_ylabel('Dec (deg)')
            ax.set_title(f'All Possible Pointings - {event_id}')
            ax.legend()
            
            output_dir = Path(config.paths.output_dir) / event_id
            plot_path = str(output_dir / 'all_pointings.png')
            plt.savefig(plot_path, dpi=300)
            plt.close()
            
            return plot_path
            
        except Exception as e:
            self.logger.error(f"Error creating pointing summary plot: {str(e)}")
            return None


