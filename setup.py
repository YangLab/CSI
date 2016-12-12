from setuptools import setup, find_packages

setup(
    name = "CSScore",
    version = "0.1",
    packages = find_packages('src'),
    package_dir = {'':'src'},
    license = 'MIT',
    install_requires=[
        'pysam>=0.8.4pre'
    ],
    entry_points={
        'console_scripts': [
            'CSScore = core.cs_xm:main'
        ]
    }
)
