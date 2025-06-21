"""
Configuration models using Pydantic for validation and type safety
"""

from typing import Dict, List, Optional
from pydantic import BaseModel, Field
from pathlib import Path

class AlertReceiverConfig(BaseModel):
    """Configuration for AlertReceiver"""
    kafka_topic: str = Field(default="gcn.classic.voevent", description="Kafka topic to consume alerts from")
    consumer_timeout: int = Field(default=30*24*60*60, description="Consumer timeout in seconds")
    alert_types: List[str] = Field(default=["PRELIMINARY", "INITIAL", "UPDATE"], description="Types of alerts to process")
    
    class Config:
        title = "AlertReceiver Configuration"

class FacilityConfig(BaseModel):
    """Configuration for observation facilities"""
    name: str = Field(description="Facility name")
    exposure_times: Dict[str, float] = Field(description="Exposure times in minutes")
    filters: List[str] = Field(description="Available filters")
    depth: float = Field(description="Base depth")
    exptime: float = Field(description="Base exposure time in minutes")
    
    class Config:
        title = "Facility Configuration"

class ObservationConfig(BaseModel):
    """Configuration for observations"""
    confidence_limits: Dict[str, float] = Field(
        default={"90": 0.9, "95": 0.95, "99": 0.99},
        description="Confidence limits for observation planning"
    )
    max_exposure: float = Field(default=3600, description="Maximum exposure time in seconds")
    min_separation: float = Field(default=0.1, description="Minimum separation in degrees")
    
    class Config:
        title = "Observation Configuration"

class PathConfig(BaseModel):
    """Configuration for file paths"""
    base_dir: Path = Field(description="Base directory for the project")
    data_dir: Path = Field(description="Directory for data files")
    output_dir: Path = Field(description="Directory for output files")
    log_dir: Path = Field(description="Directory for log files")
    skymap_dir: Path = Field(description="Directory for skymap files")
    catalog_dir: Path = Field(description="Directory for catalog files")
    
    class Config:
        title = "Path Configuration"

class SlackConfig(BaseModel):
    """Configuration for Slack notifications"""
    channel: str = Field(default="#gecko--alert", description="Slack channel for notifications")
    token_env_var: str = Field(default="SLACK_API_TOKEN", description="Environment variable for Slack API token")
    
    class Config:
        title = "Slack Configuration"

class ProcessingConfig(BaseModel):
    """Configuration for processing parameters"""
    max_retries: int = Field(default=3, description="Maximum number of retries")
    retry_delay: int = Field(default=60, description="Delay between retries in seconds")
    max_processing_time: int = Field(default=3600, description="Maximum processing time in seconds")
    
    class Config:
        title = "Processing Configuration"

class Config(BaseModel):
    """Main configuration model"""
    alert_receiver: AlertReceiverConfig = Field(default_factory=AlertReceiverConfig)
    facilities: List[FacilityConfig] = Field(default_factory=list)
    observations: ObservationConfig = Field(default_factory=ObservationConfig)
    paths: PathConfig = Field(default_factory=PathConfig)
    processing: ProcessingConfig = Field(default_factory=ProcessingConfig)
    slack: SlackConfig = Field(default_factory=SlackConfig)
    
    class Config:
        title = "GeckoDigestor Configuration"
        validate_assignment = True
