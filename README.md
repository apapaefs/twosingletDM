# twosingletDM
TRSM + Dark Matter + ElectroWeak Baryogenesis vs. Higgs Boson Pair production

# Instructions:

## Download MG5, e.g. 2.9.22 (https://launchpad.net/mg5amcnlo)
- Copy ```loop_sm_twoscalar_generic.tar.gz``` into MG5_aMC_2_9_22/models and untar it: ```tar xvzf loop_sm_twoscalar_generic.tar.gz```
- Launch MG5 and generate the process: 
```
./bin/mg5_aMC
import model loop_sm_twoscalar_generic
generate g g > h h [noborn=QCD]
output gg_hh_twoscalar
launch
```
- Enter and proceed to next screen, edit run card and change to the desired beam energies. -
- Exit and edit ```MG5_aMC_2_9_22/input/mg5_configuration.txt``, changing:
```
automatic_html_opening = False
```
- In the ```generate_mg5_trsm_xsecs.py``` script, change the following to the absolute directory of MG5:
```
MGLocation = '/home/apapaefs/Projects/TwoSingletDM/twosingletDM/MG5_aMC_v2_9_22/'
```
- Make sure ```ProcLocation['hh'] = 'gg_hh_twoscalar/'``` is set to the directory in which you have outputted the process.
## Get HiggsTools: https://gitlab.com/higgsbounds/higgstools.git and compile it:
in the HiggsTools directory:
```
mkdir build; cd build; cmake ..; make -j12
```
- You may need to ```make install```, otherwise make sure HiggsTools is in your PYTHONPATH
- You also need to place the HiggssBounds and HiggsSignals datasets in the twosingletDM directories: https://gitlab.com/higgsbounds/hbdataset and https://gitlab.com/higgsbounds/hsdataset.
## Execute:
- To begin the scan, execute ```generate_trsm_points.py SEED``` where the SEED is an integer to be used as a seed for the random numbers. 


## Install micrOMEGAs 6.1.15 and create TRSM:

- ```cd DM``` and untar `micromegas_6.1.15.tar`: ```tar xvzf micromegas_6.1.15.tar```

- ```cd micromegas_6.1.15/``` then follow the READ.ME instructions I and V:
```[g]make``` and ```./newProject TRSM```

- copy the desired model files from DM/models into micromegas. I use h4GOn:
```cp ../models/h4GOn/* TRSM/work/models/```

-  in TRSM/ there is a main.c file and a Makefile. Compile with:
```make main=main.c``` or just ```make```
The setup is complete and and micromegas is used by codes in DM/example.

- To test manually, copy the example data point:
```cp DM/data.par micromegas_6.1.15/TRSM/```
```./main data.par```
compare the output with DM/example_test.out

## run micrOMEGAs with a steering code:
- check out DM/example_steer/README.md and create another directory for the study you'd like to perform.