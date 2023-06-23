import setuptools

with open("scrimp/README.md", "r") as fh:

    long_description = fh.read()

setuptools.setup(

    name="scrimp", # Replace with your username

    version="1.0.0",

    author="G.Ferraro",

    author_email="giuseppe.ferraro@isae-supaero.fr",

    description="<Template Setup.py package>",

    long_description=long_description,

    long_description_content_type="text/markdown",

    url="https://github.com/g-haine/scrimp",

    packages=setuptools.find_packages(),

    classifiers=[

        "Programming Language :: Python :: 3",

        "License :: OSI Approved :: GNU General Public License version 3",

        "Operating System :: OS Independent",

    ],

    python_requires='>=3.6',

)
