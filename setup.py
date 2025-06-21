from setuptools import setup, find_packages

setup(
    name="geckodigestor",
    version="0.1.0",
    description="Gecko Digestor for gravitational wave event processing and observation planning",
    author="Your Name",
    author_email="your.email@example.com",
    packages=find_packages('src'),
    package_dir={'': 'src'},
    install_requires=[
        'numpy>=1.21.0',
        'pandas>=1.3.0',
        'astropy>=5.0.0',
        'ligo.skymap>=0.6.0',
        'pydantic>=2.0.0',
        'slack_sdk>=3.20.0',
        'matplotlib>=3.4.0',
        'gcn-kafka>=0.3.3'
    ],
    entry_points={
        'console_scripts': [
            'alert-receiver=geckodigestor.observers.alert_receiver:main',
            'geckodigestor=geckodigestor.core.digestor:main',
            'geckodigestor-process=geckodigestor.cli.manual_processor:main'
        ],
    },
    python_requires='>=3.8',
    package_data={
        'geckodigestor': [
            'config/*',
            'templates/*',
            'catalogs/*'
        ]
    },
    include_package_data=True
)
