from setuptools import setup, find_packages

with open("requirements.txt", "r") as f:
    requirements = f.read().splitlines()

setup(
	name="ZTF-Check",
	version="0.1.4",
	author='Joseph Guidry',
	author_email="jaguidry@bu.edu",
	description="ZTF Check! The quick way to check out your target's FOV in ZTF.",
    url='https://github.com/astrojoeg/ZTF-Check',
    install_requires=requirements,
	packages=find_packages(),
	entry_points ={
            'console_scripts': [
                'ztfcheck = ztfcheck.ztfcheck:main'
            ]
        },
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Astronomy"
    ]
)
