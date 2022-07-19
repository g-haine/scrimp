# How to install SCRIMP

- Clone the GitHub repository in your folder `/path/to/scrimp/`:
```
cd /path/to/scrimp/
git clone https://github.com/g-haine/scrimp.git
```

- Install **GMSH** client (not from the package manager! Outdated version!): https://gmsh.info/

- Install **conda**: https://docs.conda.io

- Install **mutliphenics**: https://github.com/multiphenics/multiphenics in your /path/to/scrimp/src/modules
```
mkdir /path/to/scrimp/src/modules/multiphenics
cd /path/to/scrimp/src/modules/multiphenics
git clone https://github.com/multiphenics/multiphenics.git
python setup.py install
```

- Create the *scrimp* environment:
```
conda install --file "/path/to/your/scrimp/directory/data/scrimp.yml"
```

# How to get started using SCRIMP

- Activate the environment:
```
conda activate scrimp
```

- Import the module by starting your script.py with:
```
rootdir = "/path/to/scrimp/"
import sys
sys.path.append(rootdir)
import scrimp
scrimp.setup_path(rootdir)
```

# How to test SCRIMP

By running the provided Notebooks.

- In your scrimp folder `/path/to/scrimp/` (do not forget to activate the scrimp environment before):
```
jupyter notebook &
```

- In [Jupyter](https://jupyter.org/), choose your Notebook.

**The very first notebook to test: `/path/to/scrimp/notebooks/2D/wave.ipynb`**

___
