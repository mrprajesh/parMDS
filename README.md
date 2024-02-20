# parMDS
[![seqMDS](https://github.com/mrprajesh/parMDS/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/mrprajesh/parMDS/actions/workflows/ubuntu.yml)

Effective Parallelization of the Vehicle Routing Problem at [https://dl.acm.org/doi/10.1145/3583131.3590458](https://dl.acm.org/doi/10.1145/3583131.3590458).


## Publication

**Effective Parallelization of the Vehicle Routing Problem**, 
<ins>Rajesh Pandian M </ins>, Somesh Singh, Rupesh Nasre. and N.S. Narayanaswamy,
*Genetic and Evolutionary Computation Conference* **(GECCO)**, 2023. Pages 1036–1044, Lisbon, Portugal. 15-19 July, 2023.
 
 [(DOI)](https://dl.acm.org/doi/10.1145/3583131.3590458) [(Slides)](https://mrprajesh.co.in/pdfs/gecco-cvrp-v3.pdf) [(Video)](https://youtu.be/IWgqRR-UO6U) [(Code)](https://github.com/mrprajesh/parMDS)


### Quick 15-min Video

[![YouTube Link](http://img.youtube.com/vi/IWgqRR-UO6U/0.jpg)](http://www.youtube.com/watch?v=IWgqRR-UO6U "YouTube Link")



## Requirements  

- Should work on every Linux Distribution. Tested on Ubuntu 22.04.
- For seqMDS (sequential version): GCC 7+. Tested on GCC 9.3.1
- For parMDS (parallel OpenMP): We use `nvc++` compiler from [NVIDIA's HPC SDK](https://developer.nvidia.com/hpc-sdk).

## How to run

### Build and run the executables
```
## To compile. Runs toy example on parMDS & seqMDS 

make

## To build and run only seqMDS.
make seqMDS

## To run the executable

./seqMDS.out toy.vrp [-round 0 or 1 DEFAULT:1 means round it!]
./parMDS.out toy.vrp [-nthreads <n> DEFAULT is 20] [-round 0 or 1 DEFAULT:1]


## An example
./parMDS.out inputs/Antwerp1.vrp -nthreads 16 -round 1
./seqMDS.out toy.vrp -round 0

## Output Description

//One line stats printed on std err whereas the solution on std out
inputs/Antwerp1.vrp Cost 556105 548740 518042 Time(seconds) 0.599911 2.18137 2.19818 parLimit 16 VALID    
Route #1: ..
Route #2: ..
...
Route #k:
Cost <val>

```

### Build and run the paper's artifact

```
bash ./runRounds.sh
```


- Runs both seqMDS and parMDS on all 130 inputs using round conventions.
- Create .sol files for each instance.
- Create a folder `output*` and place all .log and .sol files.
- Create a time.txt file which contains the Cost and Time of all instances.
- In addition to it, Cost and Time are printed on the terminal as well.



## Authors 
 * Rajesh Pandian M | https://mrprajesh.co.in
 * Somesh Singh     | https://ssomesh.github.io
 * Rupesh Nasre     | www.cse.iitm.ac.in/~rupesh
 * N.S.Narayanaswamy| www.cse.iitm.ac.in/~swamy


# LICENSE
[![License](http://img.shields.io/:license-mit-blue.svg?style=flat-square)](http://badges.mit-license.org)
