import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="linnaeo",
    version="0.2.0",
    author="Aaron Wolfe",
    author_email="wolfe.aarond@gmail.com",
    description="A tool for creating and categorizing protein alignments",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/beowulfey/linnaeo",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL3",
        "Operating System :: Windows",
    ],
    python_requires='>=3.7',
)
