# Read Me

This folder contains the codes for computing and testing the implementation of the eigenstates of the Schroedinger-Poisson problem.
1. `ComputeStationary` contains the code to compute the $n$-th eigenstate, with fixed discretization step. Discretization steps are quite rough, to allow for faster computation times (about 60s). Results are saved as .mat files.
2. `TestSensibility` contains the code to compute the $n$-th stationary state, varying the discretization step, so to test the sensibility of the algorithm to this parameter. Computation times are pretty wide, to allow for good resolutions. Results obtained from highest discretisations are used as a reference for elaboration and printing. Results are saved as .mat files.
3. `PrintResults` contains the code to print to .dat file the results for the $n$-th eigenstate.  
4. `ElaborateResults` contains the code to elaborate (plot and fit) the results obtained from `TestSensibility` or `ComputeStationary`, for the $n$-th stationary state. All plots are saved in the `Figures` subfolder.   
Notice that all the elaborations (except for slices plots) are repeated and refined in a separated PYTHON code. 
5. `WithSource` contains the code to analyze the excited stationary states in presence of one external source (`OneSource`).

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


  Notice that all results elaborations (except for slices plots) are repeated and refined in a separated `PYTHON` code. 


## 5. The `WithSource` folder
This folder contains the code to analyze the excited stationary states in presence of one external source (`OneSource`).

### 5.1. The `OneSource` folder
The `oneSource` folder computes the excited stationary states of the Schr\"odinger-Poisson problem in presence of  one external source. Three different  external density profiles are tested:
- an exponential profile: $\rho_{exp}(r) = a e^{-\frac{r}{r_0}}$
- a hard sphere profile: 
    $ \rho_{hs}(r) = 
    \begin{cases}
        a \qquad r\le r_0\\
        0 \qquad r > r_0\\        
    \end{cases} $
- a truncated Plummer potential (as in [Ji, Sin - 1994](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.50.3655)):
    $\rho_{tp}(r) = 
    \begin{cases}
        a \Big(\big(\frac{r}{r_0}\big)^2+1\Big)^{-5/2} \qquad\,\,\,\, r\le r_0\\
        \rho_{tp}(r_0) e^{r_0-r} \qquad\qquad\quad r > r_0\\ 
    \end{cases}$

For each of these profile, the two parameters $r_0$ and $a$ are varied, and the corresponding density profile, the eigenfunction, the eigenpotential, the eigenvelocity, and the total velocity are computed and plotted as follows:
- the eigenstate $(f_n(r),\phi_n(r))$ solves the stationary Schr\"odinger-Poisson problem with the unitary normalization on the eigenfunction:
```math
 \begin{cases}
        \triangle f = (\varepsilon + 2\phi ) f \\
        \triangle \phi =  f^2 + \rho \\
         \int f^2(x) \,d^3 x=1
    \end{cases}
```
- the eigenvelocity $v_n(r)$ and the velocity component due to the external source are computed as:
```math
v_n(r) = \sqrt{ \frac{ \int_0^{r} f_n^2(s)\,s^2 ds}{r}}\,,
\qquad 
v_\rho(r) = \sqrt{ \frac{ \int_0^{r} \rho(s)\,s^2 ds}{r}}\,,
```
- the total velocity reads:
```math
 v_{Tot}(r) =\sqrt{ v_n^2(r) + v_\rho^2(r)}
```

Computed quantities and plots are saved in dedicated subfolders.


The code is structured as follows:
- the `main.mat` file includes:
    - a parameter section
    - an optional sections for computing the desired quantities (density profiles, eigenstates, velocities) without source
    - plots in absence of external sources
    - an optional sections for computing the desired quantities (density profiles, eigenstates, velocities)  as $r_0$ varies   
    - plots as $r_0$ varies 
    - an optional sections for computing the desired quantities (density profiles, eigenstates, velocities)  as $a$ varies   
    - plots as $a$ varies.

  Computations are performed through the `computeStationaryWithSource` function, while the plots are build on loaded data, so they can be run without repeating the computations. 

- the `computeStationaryWithSource` function takes as input the type and parameters of the density and computes the density profile, solves the stationary SP with source to find the eigenfunction, eigenpotential and eigenvelocity, computes the source velocity contribution and the total velocity. 
The input parameters ($r_0$, $a$) are expressed according to the $4\pi$ normalization, as required in the solver (`nthSolver` function). The output of the solver and the original input data are then converted to unitary normalization, and the velocities are computed accordingly. 
The data (unitary normalization) are saved to a `.mat` file. 



