#!/bin/bash

cd 1D

for i in $(ls -d */)
    do
        FILE="$i${i%%/}.ipynb"
        if [[ -f "$FILE" ]]
            then
                echo "$FILE already exists."
            else
                echo "$FILE does not exist. Creation!"
                touch "$FILE"
        fi
    done

cd ../2D

for i in $(ls -d */)
    do
        FILE="$i${i%%/}.ipynb"
        if [[ -f "$FILE" ]]
            then
                echo "$FILE already exists."
            else
                echo "$FILE does not exist. Creation!"
                touch "$FILE"
        fi
    done

cd ../3D

for i in $(ls -d */)
    do
        FILE="$i${i%%/}.ipynb"
        if [[ -f "$FILE" ]]
            then
                echo "$FILE already exists."
            else
                echo "$FILE does not exist. Creation!"
                touch "$FILE"
        fi
    done

cd ../finite_dim

for i in $(ls -d */)
    do
        FILE="$i${i%%/}.ipynb"
        if [[ -f "$FILE" ]]
            then
                echo "$FILE already exists."
            else
                echo "$FILE does not exist. Creation!"
                
                touch "$FILE"
        fi
    done
