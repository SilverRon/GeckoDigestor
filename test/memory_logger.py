import os
import psutil
import time
import json
from datetime import datetime

class MemoryLogger:
    def __init__(self, log_dir='logs'):
        self.log_dir = log_dir
        os.makedirs(log_dir, exist_ok=True)
        self.log_file = os.path.join(log_dir, f'memory_log_{datetime.now().strftime("%Y%m%d_%H%M%S")}.jsonl')
        
    def log_memory(self, section_name):
        """Log current memory usage with timestamp"""
        process = psutil.Process()
        mem_info = process.memory_info()
        
        log_entry = {
            'timestamp': datetime.now().isoformat(),
            'section': section_name,
            'memory_usage': {
                'rss': mem_info.rss / (1024 * 1024),  # Convert to MB
                'vms': mem_info.vms / (1024 * 1024),  # Convert to MB
                'shared': mem_info.shared / (1024 * 1024),  # Convert to MB
                'rss_percent': psutil.virtual_memory().percent
            }
        }
        
        with open(self.log_file, 'a') as f:
            f.write(json.dumps(log_entry) + '\n')
        
        return log_entry
    
    def get_memory_usage(self):
        """Get current memory usage without logging"""
        process = psutil.Process()
        mem_info = process.memory_info()
        return {
            'rss': mem_info.rss / (1024 * 1024),  # Convert to MB
            'vms': mem_info.vms / (1024 * 1024),  # Convert to MB
            'shared': mem_info.shared / (1024 * 1024),  # Convert to MB
            'rss_percent': psutil.virtual_memory().percent
        }
    
    def print_memory_usage(self, section_name):
        """Print current memory usage"""
        usage = self.get_memory_usage()
        print(f"\nMemory Usage in {section_name}:")
        print(f"RSS: {usage['rss']:.2f} MB")
        print(f"VMS: {usage['vms']:.2f} MB")
        print(f"Shared: {usage['shared']:.2f} MB")
        print(f"RSS Percentage: {usage['rss_percent']}%")
        return usage
