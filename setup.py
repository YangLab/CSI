from setuptools import setup, find_packages

setup(
    name = "CSScore",
    version = "1.0",
    packages = find_packages('src'),
    package_dir = {'':'src'},
    license = 'MIT',
    install_requires=[
        'pysam>=0.8.4pre'
    ],
    entry_points={
        'console_scripts': [
            'CSI = core.cs_xm:main'
        ]
    }
)
