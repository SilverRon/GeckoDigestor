import logging
from pathlib import Path
from typing import Optional

from ..config.config_manager import GeckoConfig
config = GeckoConfig()

def setup_logger(name: str, log_file: Optional[str] = None) -> logging.Logger:
    """Set up a logger with file and console handlers"""
    
    # Create logger
    logger = logging.getLogger(name)
    logger.setLevel(config.log_level)
    
    # Create console handler
    console_handler = logging.StreamHandler()
    console_handler.setLevel(config.log_level)
    
    # Create file handler if log_file is specified
    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(config.log_level)
    
    # Create formatter and add it to the handlers
    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    console_handler.setFormatter(formatter)
    if log_file:
        file_handler.setFormatter(formatter)
    
    # Add handlers to logger
    logger.addHandler(console_handler)
    if log_file:
        logger.addHandler(file_handler)
    
    return logger

# Create a root logger for the package
logger = setup_logger('geckodigestor', config.log_file)
