"""
Utility functions for Slack notifications and file uploads
"""

import logging
import os
from typing import Optional, Dict, Any
from pathlib import Path

import requests
from slack_sdk import WebClient
from slack_sdk.errors import SlackApiError

from ..config.config_manager import GeckoConfig
config = GeckoConfig()

logger = logging.getLogger('geckodigestor.slack')

class SlackNotifier:
    """Class for handling Slack notifications and file uploads"""
    
    def __init__(self, token: Optional[str] = None):
        """Initialize the SlackNotifier
        
        Args:
            token: Slack API token (optional, will use config if not provided)
        """
        self.token = token or os.environ.get('SLACK_API_TOKEN')
        if not self.token:
            raise ValueError("Slack API token must be provided or set in environment")
            
        self.client = WebClient(token=self.token)
        
    def post_thread_with_file(
        self,
        channel: str,
        text: str,
        file_path: str,
        comment_text: Optional[str] = None,
        comment_file_path: Optional[str] = None
    ) -> Dict[str, Any]:
        """Post a message as a new thread with a file attachment, and optionally add a comment
        
        Args:
            channel: Channel ID or name
            text: Main thread text content
            file_path: Path to PNG file to upload
            comment_text: Optional text for comment under thread
            comment_file_path: Optional path to file for comment
            
        Returns:
            Dictionary containing thread timestamp and message/file results
        """
        try:
            # Upload file first
            file_upload_result = self._upload_file(file_path)
            
            # Post initial thread message with file
            thread_result = self._post_thread_message(
                channel,
                text,
                file_upload_result['file']['id']
            )
            
            # If comment is provided, post it under the thread
            if comment_text or comment_file_path:
                comment_result = self._post_thread_comment(
                    channel,
                    thread_result['ts'],
                    comment_text,
                    comment_file_path
                )
            else:
                comment_result = None
                
            return {
                'thread_ts': thread_result['ts'],
                'thread_message': thread_result,
                'file_upload': file_upload_result,
                'comment': comment_result
            }
            
        except SlackApiError as e:
            logger.error(f"Slack API error: {e.response['error']}")
            raise
        except Exception as e:
            logger.error(f"Error posting to Slack: {str(e)}")
            raise
    
    def _upload_file(self, file_path: str) -> Dict[str, Any]:
        """Upload a file to Slack
        
        Args:
            file_path: Path to file to upload
            
        Returns:
            Slack file upload response
        """
        try:
            file_path = str(Path(file_path).resolve())
            
            with open(file_path, 'rb') as f:
                response = self.client.files_upload(
                    file=f,
                    filename=os.path.basename(file_path)
                )
                
            return response.data
            
        except Exception as e:
            logger.error(f"Error uploading file: {str(e)}")
            raise
    
    def _post_thread_message(
        self,
        channel: str,
        text: str,
        file_id: str
    ) -> Dict[str, Any]:
        """Post a message as a new thread
        
        Args:
            channel: Channel ID or name
            text: Message text
            file_id: ID of the uploaded file
            
        Returns:
            Slack message response
        """
        try:
            response = self.client.chat_postMessage(
                channel=channel,
                text=text,
                files=[file_id]
            )
            return response.data
            
        except Exception as e:
            logger.error(f"Error posting thread message: {str(e)}")
            raise
    
    def _post_thread_comment(
        self,
        channel: str,
        thread_ts: str,
        text: Optional[str] = None,
        file_path: Optional[str] = None
    ) -> Dict[str, Any]:
        """Post a comment under a thread
        
        Args:
            channel: Channel ID or name
            thread_ts: Thread timestamp
            text: Comment text
            file_path: Path to file to attach
            
        Returns:
            Slack message response
        """
        try:
            # Prepare message and file
            kwargs = {
                'channel': channel,
                'thread_ts': thread_ts,
                'text': text if text else "Follow-up comment"
            }
            
            if file_path:
                file_upload = self._upload_file(file_path)
                kwargs['files'] = [file_upload['file']['id']]
                
            response = self.client.chat_postMessage(**kwargs)
            return response.data
            
        except Exception as e:
            logger.error(f"Error posting thread comment: {str(e)}")
            raise

def send_slack_notification(
    channel: str,
    text: str,
    file_path: str,
    comment_text: Optional[str] = None,
    comment_file_path: Optional[str] = None
) -> Dict[str, Any]:
    """Send a notification to Slack with thread and file support
    
    Args:
        channel: Channel ID or name
        text: Main thread text content
        file_path: Path to PNG file to upload
        comment_text: Optional text for comment under thread
        comment_file_path: Optional path to file for comment
        
    Returns:
        Dictionary containing thread timestamp and message/file results
    """
    try:
        notifier = SlackNotifier()
        return notifier.post_thread_with_file(
            channel,
            text,
            file_path,
            comment_text,
            comment_file_path
        )
        
    except Exception as e:
        logger.error(f"Error sending Slack notification: {str(e)}")
        raise
