# Read Me

This folder contains the codes for computing and testing the implementation of the eigenstates of the Schroedinger-Poisson problem.
Each of the following folders is authonomous, covering a specific part of the analysis.
- `ComputeStationary` contains the code to compute the $n$-th eigenstate, with fixed discretization step. Discretization steps are quite rough, to allow for faster computation times (about 60s). Results are saved as .mat files.
- `TestSensibility` contains the code to compute the $n$-th stationary state, varying the discretization step, so to test the sensibility of the algorithm to this parameter. Computation times are pretty wide, to allow for good resolutions. Results obtained from highest discretisations are used as a reference for elaboration and printing. Results are saved as .mat files.
- `PrintResults` contains the code to print to .dat file the results for the $n$-th eigenstate.  
- `ElaborateResults` contains the code to elaborate (plot and fit) the results obtained from `TestSensibility` or `ComputeStationary`, for the $n$-th stationary state. All plots are saved in the `Figures` subfolder.   
Notice that all the elaborations (except for slices plots) are repeated and refined in the `ExcitedScalings` code.
- `ExcitedScalings` contains the systematic analysis on the eigenstates to characterize how they scale with the excitation index $n$.

## 1. The `ComputeStationary` folder
This folder contains the code to compute the $n$-th eigenstate, with fixed discretization step. Discretization steps are quite rough, to allow for faster computation times (about 60s). Results are saved as .mat files.

The solver implements the procedure described in [*Bernstein, Giladi, Jones (1998)*](https://doi.org/10.1142/S0217732398002473), computing the spherically symmetric eigenstates of the Schroedinger-Poisson problem. The notation adopted in the implementation is slightly different with respect to that adopted in *Bernstein, Giladi, Jones (1998)*. We use:
```math
  \begin{cases}
      \triangle f = (\varepsilon + 2\phi ) f 
      \qquad \qquad \text{with } \int |f(x)|^2 d^3 x=4 \pi,
      \\
      \triangle \phi =  |f|^2 \\
  \end{cases}
```
and the `kthSolver.m` function computes the solution $(\varepsilon_n,f_n,\phi_n)$ relative to this notation.<br>
After computing the eigenstate, `main_computeStationary.m` allows to scale the results to match *Bernstein, Giladi, Jones (1998)*:
```math
  \begin{cases}
      \triangle f_B = 2(-\varepsilon_B +\phi_B ) f_B 
      \qquad \qquad \text{with } \int |f_B(x)|^2 d^3 x =1,
      \\
      \triangle \phi_B =  4\pi|f_B|^2 \\
  \end{cases}
```
for direct comparison with their tabulated results.<br>
Finally, `main_computeStationary.m` gives to option (activated by default) to scale the results to switch to norm=1 notation, i.e. to solve: 
```math
  \begin{cases}
      \triangle f = (\varepsilon + 2\phi ) f 
      \qquad \qquad \text{with } \int |f(x)|^2 d^3 x=1,
      \\
      \triangle \phi =  |f|^2 \\
  \end{cases}
```
The results are saved as .mat files in the desired notation. 

## 2. The `TestSensibility` folder
This folder contains the code to compute the $n$-th eigenstatestate varying the discretization step, so to test the sensibility of the algorithm to this parameter. Computation times are pretty wide, to allow for good resolutions. Results obtained from highest discretisations are used as a reference for elaboration and printing. Results are saved as .mat files in the desired subfolder.

The folder includes:
- `main_testSensibility.m` file. It computes the $n$-th eigenstate at five different resolutions, and saves all results in a .mat file. The solver is the same as in the `ComputeStationary` folder. As explained for that folder, the output of the solver refers to the norm=$4\pi$ notation. The results are then scaled to match the norm=$1$ notation (this step, performed by default, can optionally be avoided):  
```math
  \begin{cases}
      \triangle f = (\varepsilon + 2\phi ) f 
      \qquad \qquad \text{with } \int |f(x)|^2 d^3 x=1,
      \\
      \triangle \phi =  |f|^2 \\
  \end{cases}
```
- `main_reportSensibility.m` file. It analyzes the five different versions of the $n$-th eigenstate, each obtained with a different resolution, to check the sensibility of the solution to this parameter. The analysis includes:
  - plot, as the discretization step varies, of the eigenvalues;
  - plot, as the discretization step varies, of the correction ("error") on the eigenvalue, computed with respect to the value obtained with the previous (larger) discretization step;
  - plot, as the discretization step varies, of the percentage correction ("error") on the position of the local extrema, computed with respect to the value obtained with the previous (larger) discretization step; the corrections are averaged among all the local extrema and expressed relatively to the outermost local maximum;
  - plot, as the discretization step varies, of the radius enclosing the $95\%$ of the total mass;
  - plot, as the discretization step varies, of the time elapsed while computing the solution.

All the results are also printed to log files
The analysis is aimed at judging the quality of the single $n$-th solution, and does not cross-compare the results of the sensibility tests at different values of excitation index $n$.


## 3. The `PrintResults` folder
This folder contains the code to print to file the results obtained for the $n$-th eigenstate. <br>
The code loads the previously computed solutions from a given input folder. Results may include one solution per file (as from `ComputeStationary` analysis) or more solutions per file (as from `TestSensibility`, where several solutions with different refinements are saved in the same file). In the latter case, only the last (more refined) solution is used for printing. This is considered the preferred data in our analysis. <br>
Printed results are saved as .dat files in different subfolders:
  - `neNumData` contains the eigenvalues $\varepsilon(n)$ for different number of nodes $n$, in a unique file
  - `rfNumData` contains the eigenfunctions $f(r)$ for each $n$-th stationary state, in separate files.
  - `rpNumData` contains the potential function $\phi(r)$ for each $n$-th stationary state, in separate files.
  - `rvNumData` contains the velocity function $v(r)$ for each $n$-th stationary state, in separate files.
  - `rmNumData` contains the mass function $m(r)$ for each $n$-th stationary state, in separate files.

The code allows to select an output folder and it authonomously sets the required subfolders, adding an appropriate text comment in the header of the files.

## 4. The `ElaborateResults` folder
This folder contains the code to elaborate the results obtained from `TestSensibility` or `ComputeStationary`, for the $n$-th stationary state. All plots are saved in the desired subfolder. Computed (and plotted) quantities include:
  - velocity curve $v(r)$.
  - eigenfunction $f(r)$, in linear scale and in log-log scale.
  - potential $\phi(f)$.
  - sliced 3d or 2d colored plot of $|f(x)|$ or of $f^2(x)$.
  - fit of local extrema (peaks) of the eigenfunction, in linear scale and in log-log scale.
  - mid-range fitted velocity curve.

  Notice that all results elaborations (except for slices plots) are repeated and refined in the  `ExcitedScalings` folder. 
  
## 5. The `ExcitedScalings` folder
This folder contains the systematic analysis on the eigenstates to characterize how they scale with the excitation index $n$. 
- the `NumParam.py` file summarizes all the heuristic laws derived in the analysis and the associated numerical values.
- the `excitedScalings.ipynb` notebook collects all the code performing the analysis.

### The `excitedScalings.ipynb` notebook
The notebook uses the results printed by the `PrintResults` code and analyzes them to characterize the scaling of the eigenstates with the excitation index $n$. 
The first cell is necessary for the setting. Each other cell of this notebook is designed to be independent from the others. 
The structure is the following:
- eigenvalues 
- plot of the eigenstate
  - single eigenstate
  - sample eigenstates
  - magnified views
- eigenfunctions
  - nodes patterns scaling laws
  - local extrema scaling laws
  - approximate eigenfunction
- eigenpotential
- eigenvelocity
  - scaling laws
  - universality
  - approximate eigenvelocity     
- WKB approximation 
