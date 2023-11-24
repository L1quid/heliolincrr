import setuptools

setuptools.setup(
    name='heliolincrr',
    version='0.1',
    author='Ben Engebreth',
    author_email='ben.engebreth@gmail.com',
    description='Solar system object linking using heliocentric projection and postition vector clustering at two reference epochs.',
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url='https://github.com/bengebre/heliolincrr',
    project_urls = {'Issues': 'https://github.com/bengebre/heliolincrr/issues'},
    packages=['heliolincrr'],
    install_requires=[
	'numpy>=1.20.1',
	'pandas>=1.2.4',
	'astropy>=4.2.1',
	'scikit-learn>=0.24.1',
	'poliastro @ git+https://github.com/poliastro/poliastro.git@55e96432b27301c5dffb4ef6b4f383d970c6e9c0',
	'tqdm>=4.59.0',
    	'cupy>=12.2.0']
)
