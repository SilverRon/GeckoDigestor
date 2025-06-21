"""
Manual processor for processing specific GW events with custom configurations
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, Any

from ..config.config_manager import config_manager
from ..core.gecko_digestor import GeckoDigestor
from ..utils.logger import logger

class ManualProcessor:
    """Class for manually processing GW events"""
    
    def __init__(self):
        """Initialize the ManualProcessor"""
        self.logger = logger
        self.digestor = GeckoDigestor()
        
    def process_event(self, event_id: str, config_overrides: Dict[str, Any] = None) -> None:
        """Process a specific event with optional configuration overrides
        
        Args:
            event_id: ID of the event to process
            config_overrides: Dictionary of configuration overrides
        """
        try:
            # Apply configuration overrides if provided
            if config_overrides:
                self._apply_config_overrides(config_overrides)
                
            # Get event directory
            event_dir = Path(config_manager.config.paths.data_dir) / event_id
            
            if not event_dir.exists():
                self.logger.error(f"Event directory not found: {event_dir}")
                return
                
            self.logger.info(f"Processing event: {event_id}")
            self.digestor.process_event(event_dir)
            
        except Exception as e:
            self.logger.error(f"Error processing event {event_id}: {str(e)}")
            raise
    
    def _apply_config_overrides(self, overrides: Dict[str, Any]) -> None:
        """Apply configuration overrides
        
        Args:
            overrides: Dictionary of configuration overrides
        """
        try:
            # Update configuration with overrides
            for section, values in overrides.items():
                if hasattr(config_manager.config, section):
                    config_section = getattr(config_manager.config, section)
                    for key, value in values.items():
                        if hasattr(config_section, key):
                            setattr(config_section, key, value)
            
            # Save updated configuration
            config_manager.save_config()
            
        except Exception as e:
            self.logger.error(f"Error applying configuration overrides: {str(e)}")
            raise
    
    def process_multiple_events(self, event_ids: list, config_overrides: Dict[str, Any] = None) -> None:
        """Process multiple events
        
        Args:
            event_ids: List of event IDs to process
            config_overrides: Dictionary of configuration overrides
        """
        for event_id in event_ids:
            self.process_event(event_id, config_overrides)

def main():
    """Main entry point for manual processing"""
    parser = argparse.ArgumentParser(
        description='Manual processor for GW events')
    
    # Positional arguments
    parser.add_argument('event_ids', nargs='+',
                        help='Event IDs to process (e.g., S230528ay)')
    
    # Optional arguments
    parser.add_argument('--config', type=str,
                        help='Path to JSON file containing configuration overrides')
    parser.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        default='INFO', help='Logging level')
    
    args = parser.parse_args()
    
    # Set log level
    logger.setLevel(args.log_level)
    
    # Load configuration overrides if provided
    config_overrides = None
    if args.config:
        try:
            with open(args.config) as f:
                config_overrides = json.load(f)
        except Exception as e:
            logger.error(f"Error loading configuration file: {str(e)}")
            return
    
    # Create and run processor
    processor = ManualProcessor()
    processor.process_multiple_events(args.event_ids, config_overrides)

if __name__ == '__main__':
    main()
