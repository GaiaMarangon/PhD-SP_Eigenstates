# Read Me

This folder contains the codes for computing and testing the implementation of the eigenstates of the Schroedinger-Poisson problem.
1. `ComputeStationary` contains the code to compute the $n$-th eigenstate, with fixed discretization step. Discretization steps are quite rough, to allow for faster computation times (about 60s). Results are saved as .mat files.
2. `TestSensitivity` contains the code to compute the $n$-th stationary state, varying the discretization step, so to test the sensitivity of the algorithm to this parameter. Computation times are pretty wide, to allow for good resolutions. Results obtained from highest discretisations are used as a reference for elaboration and printing. Results are saved as .mat files.
3. `PrintResults` contains the code to print to .dat file the results for the $n$-th eigenstate.  
4. `ElaborateResults` contains the code to elaborate (plot and fit) the results obtained from `TestSensitivity` or `ComputeStationary`, for the $n$-th stationary state. All plots are saved in the `Figures` subfolder.   
Notice that all the elaborations (except for slices plots) are repeated and refined in the `ExcitedScalings` code.
5. `ExcitedScalings` contains the systematic analysis on the eigenstates to characterize how they scale with the excitation index $n$.
6. `WithSource` contains the code to analyze the excited stationary states in presence of one external source (`OneSource`).

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

## 2. The `TestSensitivity` folder
This folder contains the code to compute the $n$-th eigenstatestate varying the discretization step, so to test the sensitivity of the algorithm to this parameter. Computation times are pretty wide, to allow for good resolutions. Results obtained from highest discretisations are used as a reference for elaboration and printing. Results are saved as .mat files in the desired subfolder.

The folder includes:
- `main_testSensitivity.m` file. It computes the $n$-th eigenstate at five different resolutions, and saves all results in a .mat file. The solver is the same as in the `ComputeStationary` folder. As explained for that folder, the output of the solver refers to the norm=$4\pi$ notation. The results are then scaled to match the norm=$1$ notation (this step, performed by default, can optionally be avoided):  
```math
  \begin{cases}
      \triangle f = (\varepsilon + 2\phi ) f 
      \qquad \qquad \text{with } \int |f(x)|^2 d^3 x=1,
      \\
      \triangle \phi =  |f|^2 \\
  \end{cases}
```
- `main_reportSensitivity.m` file. It analyzes the five different versions of the $n$-th eigenstate, each obtained with a different mesh size $\\{ h_i \\}_{i=1}^5$ (decreasing with $i$), to check the sensitivity of the solution to this parameter. Some reference quantities and compared, computing their difference for successive values of $h$, relative to the value of that quantity obtained with higher refinement:
  - Relative Difference of  the eigenvalues:
```math
\Delta \varepsilon_n (h_i) := \frac{\varepsilon_n(h_i)-\varepsilon_n(h_{i-1})}{\varepsilon_n(h_5)}\,,\qquad \text{for }i=2,\dots,5
```
  - Relative Difference of the outermost local extremum of the eigenfunction, $(r_{Out}, f_{Out}=|f(r_{Out})|)$:
```math
\Delta r_{Out} (h_i) := \frac{r_{Out}(h_i)-r_{Out}(h_{i-1})}{r_{Out}(h_5)}\,,
\qquad \text{for }i=2,\dots,5
```
```math
\Delta f_{Out} (h_i) := \frac{f_{Out}(h_i)-f_{Out}(h_{i-1})}{f_{Out}(h_5)}\,,
\qquad \text{for }i=2,\dots,5
```
  - Relative Difference of the outermost local extremum of the eigenvelocity, $(r_{Out}, v_{Out}=|v(r_{Out})|)$:
```math
\Delta r_{Out} (h_i) := \frac{r_{Out}(h_i)-r_{Out}(h_{i-1})}{r_{Out}(h_5)}\,,
\qquad \text{for }i=2,\dots,5
```
```math
\Delta v_{Out} (h_i) := \frac{v_{Out}(h_i)-v_{Out}(h_{i-1})}{v_{Out}(h_5)}\,,
\qquad \text{for }i=2,\dots,5
```
  - Relative Difference of  the innermost local extremum of the eigenvelocity, $(r_{inn}, v_{inn})$:
```math
\Delta r_{inn} (h_i) := \frac{r_{inn}(h_i)-r_{inn}(h_{i-1})}{r_{inn}(h_5)}\,,
\qquad \text{for }i=2,\dots,5
```
```math
\Delta v_{inn} (h_i) := \frac{v_{inn}(h_i)-v_{inn}(h_{i-1})}{v_{inn}(h_5)}\,,
\qquad \text{for }i=2,\dots,5
```
  - Relative Difference of Local Maxima of Eigenfunction (average):
```math
\frac{\sum_j |r_{max,j}(h_i)-r_{max,j}(h_{i-1})|}{ n_{max} }  \frac{1}{ r_{Out}}
```
```math
\frac{\sum_j |f_{max,j}(h_i)-f_{max,j}(h_{i-1})|}{ n_{max}  } \frac{1}{f_{0}}
```
  - Time elapsed while computing the solution.

All the results are plotted and printed to log files.
The analysis is aimed at judging the quality of the single $n$-th solution, and does not cross-compare the results of the sensitivity tests at different values of excitation index $n$.


## 3. The `PrintResults` folder
This folder contains the code to print to file the results obtained for the $n$-th eigenstate. <br>
The code loads the previously computed solutions from a given input folder. Results may include one solution per file (as from `ComputeStationary` analysis) or more solutions per file (as from `TestSensitivity`, where several solutions with different refinements are saved in the same file). In the latter case, only the last (more refined) solution is used for printing. This is considered the preferred data in our analysis. <br>
Printed results are saved as .dat files in different subfolders:
  - `neNumData` contains the eigenvalues $\varepsilon(n)$ for different number of nodes $n$, in a unique file
  - `rfNumData` contains the eigenfunctions $f(r)$ for each $n$-th stationary state, in separate files.
  - `rpNumData` contains the potential function $\phi(r)$ for each $n$-th stationary state, in separate files.
  - `rvNumData` contains the velocity function $v(r)$ for each $n$-th stationary state, in separate files.
  - `rmNumData` contains the mass function $m(r)$ for each $n$-th stationary state, in separate files.

The code allows to select an output folder and it authonomously sets the required subfolders, adding an appropriate text comment in the header of the files.

## 4. The `ElaborateResults` folder
This folder contains the code to elaborate the results obtained from `TestSensitivity` or `ComputeStationary`, for the $n$-th stationary state. All plots are saved in the desired subfolder. Computed (and plotted) quantities include:
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

### 5.1. The `excitedScalings.ipynb` notebook
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

## 6. The `WithSource` folder
This folder contains the code to analyze the excited stationary states in presence of one external source (`OneSource`).

### 6.1. The `OneSource` folder
The `oneSource` folder computes the excited stationary states of the Schroedinger-Poisson problem in presence of  one external source. Three different  external density profiles are tested:
- an exponential profile:
  ```math
  \rho_{exp}(r) = a e^{-\frac{r}{r_0}}
  ```
- a hard sphere profile: 
  ```math
   \rho_{hs}(r) = 
    \begin{cases}
        a \qquad r\le r_0\\
        0 \qquad r > r_0\\        
    \end{cases}
  ```
- a truncated Plummer potential (as in [Ji, Sin - 1994](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.50.3655)):
  ```math
  \rho_{tp}(r) = 
    \begin{cases}
        a \Big(\Big(\frac{r}{r_0}\Big)^2+1\Big)^{-5/2} \qquad\,\,\,\, r\le r_0\\
        \rho_{tp}(r_0) e^{r_0-r} \qquad\qquad\quad r > r_0\\ 
    \end{cases}
  ```

For each of these profile, the two parameters $r_0$ and $a$ are varied, and the corresponding density profile, the eigenfunction, the eigenpotential, the eigenvelocity, and the total velocity are computed and plotted as follows:
- the eigenstate $(f_n(r),\phi_n(r))$ solves the stationary Schroedinger-Poisson problem with the unitary normalization on the eigenfunction:
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



