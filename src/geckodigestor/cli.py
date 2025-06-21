import argparse
import logging
from typing import Optional

from .core import GeckoDigestor
from .observers import AlertReceiver
from .utils.logger import logger

def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Gecko Digestor - Astronomical Observation Management System'
    )
    
    parser.add_argument(
        '--mode',
        choices=['digestor', 'alert-receiver', 'test'],
        default='digestor',
        help='Operation mode'
    )
    
    parser.add_argument(
        '--config',
        type=str,
        help='Path to configuration file'
    )
    
    parser.add_argument(
        '--log-level',
        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
        default='INFO',
        help='Logging level'
    )
    
    return parser.parse_args()

def main():
    """Main entry point for the CLI"""
    try:
        args = parse_arguments()
        logger.setLevel(args.log_level)
        
        if args.mode == 'digestor':
            digestor = GeckoDigestor()
            digestor.run()
        elif args.mode == 'alert-receiver':
            receiver = AlertReceiver()
            receiver.start()
        elif args.mode == 'test':
            # Add test mode logic here
            pass
            
    except Exception as e:
        logger.error(f"Error in main execution: {str(e)}")
        raise

if __name__ == '__main__':
    main()
