from setuptools import setup

setup(
    name='snpit',
    version='1.0.0',
    author='Samuel Lipworth',
    packages=['snpit'],
    package_data={'snpit': ['../lib/*']},
    install_requires=[
        "PyVCF >= 0.6.8",
        "biopython >= 1.70"
    ],
    scripts=["bin/snpit-run.py"],
    license='unknown',
    long_description=open('README.md').read(),
    include_package_data=True,
    zip_safe=False
)
