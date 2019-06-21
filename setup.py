from distutils.core import setup
from pathlib import Path

long_description = Path(__file__).parent.joinpath("README.md").read_text()


short_description = (
    "Whole genome SNP based identification of members of the "
    "Mycobacterium tuberculosis complex."
)

package_data = list(Path(__file__).parent.joinpath("lib").rglob("*"))
package_data.append("*.md")

setup(
    name="snpit",
    version="1.0.0",
    author="Samuel Lipworth",
    packages=["snpit"],
    package_data={"snpit": package_data},
    include_package_data=True,
    install_requires=["pysam"],
    test_requirements=["pytest"],
    license="MIT",
    short_description=short_description,
    long_description=long_description,
    url="https://github.com/philipwfowler/snpit",
    entry_points={"console_scripts": ["snpit = snpit.__main__:main"]},
)
