#!/bin/bash

# sphinx-apidoc -f -e --implicit-namespaces -o /docs/source

sphinx-build -b html docs/source/ docs/build/html

sphinx-build -b latex docs/source/ docs/build/latex

cd docs/build/latex

make
