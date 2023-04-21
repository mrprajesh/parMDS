# parMDS
Effective Parallelization of the Vehicle Routing Problem, a GECCO 2023 Paper.

# Requirements  
- parMDS: [NVIDIA's HPC SDK](https://developer.nvidia.com/hpc-sdk) . Has `nvc++`  compiler
- seqMDS: GCC 7+. Tested on GCC 9.3.1
- Should work on every Linux Distribution. Tested on Ubuntu 22.04.

# How to run

```
bash ./runRounds.sh
```

- Runs both seqMDS and parMDS on all 130 inputs using round conventions.
- Create .sol files for each instances.
- Create a folder `output*` and places all .log and .sol files.
- Create a time.txt file which contains Cost and Time of all instances.
- In addition to it, Cost and Time is printed on the terminal as well.


# LICENSE
[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)


Authors: 
 * Rajesh Pandian M | https://mrprajesh.co.in
 * Somesh Singh     | https://ssomesh.github.io
 * Rupesh Nasre     | www.cse.iitm.ac.in/~rupesh
 * N.S.Narayanaswamy| www.cse.iitm.ac.in/~swamy
