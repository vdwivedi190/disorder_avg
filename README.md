# Disorder-averaging using transfer matrices

Introduction
------------

This repository contains the code to compute disorder-averaged Lyapunov exponents for two-dimensional tight-binding models using generalized transfer matrices (introduced in [Phys. Rev. B __93__, 134304 (2016)](https://doi.org/10.1103/PhysRevB.93.134304)). These can be further used to compute the transmission probability and hence conductance for the quasi-1d systems using the [Landauer formula](https://en.wikipedia.org/wiki/Landauer_formula). While the code for the computation of Lyapunov exponents should work for arbitarary systems, three examples are considered in more detail: 

__1. Chern insulator__: 
We consider the two-band model for a Chern insulator defined on the square lattice by the Bloch Hamiltonian

$$H(\mathbf{k}) = \sin k_x \sigma_1 + \sin k_y \sigma_2 + (2 - m - \cos k_x - \cos k_y) \sigma_3$$

which is gapless for $m=0,2,4$ and exhibits a Chern insulator phase for $0<m<2$ and $2<m<4$. We consider two kinds of on-site disorder: Anderson disorder and pseudomagnetic disorder. The former is realized by adding a random chemical potential $w\sigma_0$ at each site, while the latter is realized by adding $w\sigma_3$ on each site. In both cases, $w$ is drawn from a uniform distribution with mean 0 and variance $W^2/12$. 

__2. Second-order Kitaev spin liquid:__ 
We consider the generalized Kitaev model on the pentacoordinated Shastry-Sutherland lattice introduced in [Phys. Rev. B __98__, 054432 (2018)](https://doi.org/10.1103/PhysRevB.98.054432). At zero temperature, this reduces to a tight-binding model of Majorana fermions hopping in a static $\mathbb{Z}_2$ flux background, which realizes a second-order topological insulating phase. At finite temperatures, the background gauge field can fluctuate, leading to creation of vison pairs. We compute the effect of this _flux disorder_ on the transmission probability, which can be used as a probe to determine if the system is gapped/gapless. We find that a "thermal metal" phase emerges for strong enough fluctuations. 

__3. Topological insulator nanowires__ 
We consider nanowires of three-dimensional time-reversal invariant topological insulators (see [Phys. Rev. B __84__, 201105(R) (2011)](https://doi.org/10.1103/PhysRevB.84.201105) for the explicit 4-band model), which were proposed as a platform to realize Majorana fermions via proximity coupling to a s-wave superconductor. 3d topological insulators exhibt a Dirac cone on the surface protected by the bulk topology; however, for nanowires with a small cross section, the surface states acquire a finite-size gap. As one increases the Fermi level, new conduction channels along the nanowire are thus introduced, leading to steps in the conductance as a function of the Fermi level. We consider the effect of disorder on these nanowires, which has the counterintuitive effect of [_reducing_ the conductance](https://www.nature.com/articles/s41467-021-21230-3) when the Fermi level crosses a new surface band. The conductance eventually increases as the Fermi level is increased further.  



Usage
-----

The central function for this computation is `getLyapunovList` in the file `src/computeLyapunov.jl`, which computes the disorder-averaged Lyapunov exponents for a tight-binding Hamiltonian (described by the hopping matrix along the direction the transfer matrix acts in and the on-site matrix) for a given disorder type. The files `src/compute*.jl` construct these matrices given system parameters for various physical systems (as described below). To scan over a set of parameters, one can instead use `src/run*.jl`, which take as input text files with a list of parameters, generate the relevant matrices by invoking `src/compute*.jl`, and finally calls `getLyapunovList`. The computation can alternatively be performed on a HPC cluster by calling `src/run*_cluster.jl` instead (although this was written with a specific cluster in mind). 

The input text files can be generated from `genscript/genscript*.jl` by first setting the project name and requisite ranges for various system parameters therein. Running it then generates a directory with the project name (where the final computations and log files are stored) as well as a text file `job_list.txt`, which is intended as input to the `src/run*.jl` described above. The generating scripts also optionally create bash scripts that can be used to run the rest of the computation on a HPC cluster. The bash scripts are intended for the CHEOPS HPC cluster at University of Cologne, where this computation was performed. 
 
The basic workflow thus involves first editing the `genscript_*.jl` files with the ranges for various parameters one is interested in, as described at the beginning of these files. For instance, for CI (Chern Insulators), setting 
```julia 
projname = "proj_name"
``` 
inside `genscriptCI.jl` and then running it from the terminal as 
```bash 
username@~> julia disorder_avg/genscript/genscriptCI.jl  
```
then generates a directory named `proj_name` within the current directory containing a `job_list.txt` file. The computation of Lyapunov exponents then follows by a call to the corresponding `src/run*.jl` file with the directory created by `genscriptCI.jl` as input: 
```bash 
username@~> julia disorder_avg/src/runCI.jl ./proj_name
```
To run it on the cluster, one instead needs to first move the project folder to the location required by the cluster and run `runCI_cluster.jl`. The final data can be imported and plotted from the files `projname/lyaps/l_*.txt`. 
