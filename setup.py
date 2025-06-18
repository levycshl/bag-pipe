from setuptools import setup, find_packages

setup(
    name="bag-pipe",
    version="0.1.0",
    packages=find_packages(where = "src"),
    python_requires=">=3.9",
    package_dir={"": "src"},
    install_requires=[
        "read-break @ git+https://github.com/levycshl/read-break.git",
    ],
    entry_points={
            "console_scripts": [
                "bag_pipe=bag_pipe.cli:main",
            ],
    },
    author="Dan Levy",
    author_email="levy@cshl.edu",
    description="Pipeline for converting BAG FASTQ files to read tables",
)
