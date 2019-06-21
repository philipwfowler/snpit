from setuptools import setup

setup(
    name="snpit",
    version="1.0.0",
    author="Samuel Lipworth",
    packages=["snpit"],
    package_data={"snpit": ["../lib/*"]},
    install_requires=["pysam", "numpy"],
    test_requirements=["pytest"],
    scripts=["bin/snpit-run.py"],
    license="unknown",
    long_description=open("README.md").read(),
    zip_safe=False,
)
