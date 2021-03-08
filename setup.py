import setuptools
import os

setuptools.setup(

    name="tbhgvs",
    version="0.0.1",
    packages=["tbhgvs"],
    license="MIT",
    long_description="Utilities to convert between hgvs format and VCF genome coordinates",
    scripts= ["scripts/%s" % x for x in os.listdir("scripts")],
)
