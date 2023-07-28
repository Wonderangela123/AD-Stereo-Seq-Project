conda create -n ADstereo python=3.8
conda activate ADstereo

git clone -b dev https://github.com/STOmics/stereopy.git
cd stereopy
python setup.py install

conda install -c conda-forge notebook
