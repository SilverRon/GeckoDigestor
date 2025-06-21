"""
Observation scheduler module for GECKO Digestor

This module provides the ObservationScheduler class which handles scheduling of
observations based on gravitational wave alerts.
"""

from typing import Dict, Optional
from ..config.config_manager import config
from ..core.digestor import GeckoDigestor
import logging

logger = logging.getLogger(__name__)

class ObservationScheduler:
    """Class for scheduling observations based on gravitational wave alerts."""
    
    def __init__(self):
        """Initialize the observation scheduler."""
        self.digestor = GeckoDigestor()
        self.facilities = self._load_facilities()
        
    def _load_facilities(self) -> Dict:
        """Load observation facilities from configuration."""
        facilities = {}
        try:
            with open(config.facilities_file, 'r') as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        facility = line.strip().split('\t')
                        facilities[facility[0]] = {
                            'name': facility[0],
                            'location': facility[1],
                            'type': facility[2]
                        }
        except Exception as e:
            logger.error(f"Error loading facilities: {e}")
        return facilities
    
    def schedule_observations(self, alert: Dict) -> Optional[Dict]:
        """
        Schedule observations based on alert data.
        
        Args:
            alert: Dictionary containing alert data
            
        Returns:
            Dictionary containing observation schedule or None if scheduling failed
        """
        try:
            # Process alert data
            processed_alert = self.digestor.process_alert(alert)
            
            # Create observation schedule
            schedule = {
                'alert_id': processed_alert['superevent_id'],
                'alert_type': processed_alert['alert_type'],
                'facilities': [],
                'timestamp': processed_alert['time_created']
            }
            
            # Add facilities to schedule
            for facility in self.facilities.values():
                schedule['facilities'].append({
                    'name': facility['name'],
                    'status': 'pending',
                    'scheduled_time': None
                })
            
            return schedule
            
        except Exception as e:
            logger.error(f"Error scheduling observations: {e}")
            return None
