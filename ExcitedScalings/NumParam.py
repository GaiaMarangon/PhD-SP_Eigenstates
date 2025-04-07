"""
NUMERICAL PARAMETERS:
Describe the dependence on n (number of nodes) of certain features of:
- the eigenvalue e_n(r), 
- the eigenstate f_n(r), 
- the potential phi_n(r), 
- the velocity curves v(r).

The analyzed features and the corresponding heuristic laws are as follows: 
1. EIGENVALUES:
- eigenvalues: eig(n) = b *(n+n0)^a

2. EIGENFUNCTION:
- radial position of innermost node: rOn(n) = an^2 +bn +c
- radial position of outermost node: rIn(n) = mn +q

- central value: f0(n)  = beta n^alpha
- starting point of interpolation: rStart(n) = mn +q
- value at starting point of interpolation: fStart(n) = beta ( rStart )^alpha 
- minimum point of interpolation: rMin(n) = an^2 +bn +c
- value at minimum point of interpolation: fMin(n) = beta ( rMin )^alpha 
- radial position of outermost local maxima: rOut(n) = an^2 +bn +c
- value at outermost local maxima: fMin(n) = beta ( rOut )^alpha 
- exponent parameter in fit of local maxima: a(n) = beta*n^alpha -1 
- amplitude parameter in fit of local maxima: b(n) = beta*n^alpha

3. POTENTIAL:
- central value: phi0(n) = beta n^alpha
- slope of short range log fit: m(n) = beta n^alpha
- intercept of short range log fit: q(n)= beta n^alpha
- exponent of long range exp fit: a(n) = beta n^alpha
- amplitude of long range exp fit: b(n)= beta n^alpha

4. VELOCITY:
- radius   of the outermost local maximum: rOut(n) = an^2 +bn +c
- velocity of the outermost local maximum: vOut(n) = beta ( an^2+bn+c )^alpha
- radius   of the innermost local maximum: rInn(n) = mn +q
- velocity of the innermost local maximum: vInn(n) = beta ( mn +q )^alpha
- slope     of the average linear growth in the mid-range: slope(n)  = beta n^alpha
- intercept of the average linear growth in the mid-range: interc(n) = beta n^alpha

All fits are performed on set of numerical points spanned by: nVec = np.concatenate( (  np.arange(5,31)  ,np.array([40,50,60,70,80])  )) 
For function f(r) interpolations, the power law fit is done with: percent=0.95
"""

#-- EIGENVALUE PARAMETERS -----------------------------------------------------------
a_eigval  = -2.011561835646143 
b_eigval  = 0.0012385560816246498
n0_eigavl =  0.77625101239998
#-- EIGENVALUE PARAMETERS -----------------------------------------------------------





#-- FUNCTION PARAMETERS -------------------------------------------------------------

#-- NODES ------------------
a_on = 129.5579226100235 
b_on = -124.797500257052 
c_on = 794.8017661949074 

m_in = 56.32983017266728 
q_in = 187.50353319128308 

#-- LOCAL MAXIMA ------------------
alpha_f0 = -1.8816527710950053
beta_f0 = 0.0006656232814537857

m_rStart = 82.45289820017001
q_rStart = 288.8862971738147

alpha_fStart = -2.1807604832275316
beta_fStart = 7.901127094028111

a_rMin = 105.15141802166829
b_rMin = 207.39656202468987
c_rMin = -372.1593926191778

alpha_fMin = -1.5070079778077055
beta_fMin = 0.4210641640905319

a_rOut_f = 130.9108780617742
b_rOut_f = 53.53048936165272
c_rOut_f = 340.3178320999541

alpha_fOut = -1.4081570481819001
beta_fOut = 0.2045280226107504

alpha_alpha = -0.30304933762326847
beta_alpha = 0.2822924971086222

alpha_beta = -0.8489942388176521
beta_beta = 0.005298487552979472

m_beta = 0.010627460785413349
q_beta = 0.009877677361839338
#-- FUNCTION PARAMETERS -------------------------------------------------------------



#-- POTENTIAL PARAMETERS ------------------------------------------------------------
alpha_phi0 = -1.81650481e+00
beta_phi0  = -1.02181702e-03

alpha_m = -1.97378494e+00
beta_m  = 2.20333088e-04

alpha_q = -1.81468612e+00
beta_q  = -1.86434655e-03

alpha_a = -1.83548740e+00
beta_a  = -3.49529818e-03

alpha_b = -1.84333716e+00
beta_b  = 7.80929594e-04
#-- POTENTIAL PARAMETERS ------------------------------------------------------------



#-- VELOCITY PARAMETERS -------------------------------------------------------------
# - parameters for:   rOut(n) = an^2 +bn +c
a_rOut = 132.98418795903186
b_rOut = 244.6101536051851
c_rOut = -185.32041639205738

# - parameters for:   vOut(n) = beta ( an^2+bn+c )^alpha
alpha_vOut = -0.49569090723881976
beta_vOut  = 0.2657948172660721

# - parameters for:   rInn(n) = mn +q
m_rInn = 40.092911709640575
q_rInn = 124.33889924566756

# - parameters for:   vInn(n) = beta ( mn +q )^alpha
alpha_vInn = -1.1435286503688644
beta_vInn  = 1.9410777328447457

# - parameters for:   slope(n) = beta n^alpha
alpha_slope = -2.86468831e+00
beta_slope  = 2.81644734e-05

# - parameters for:   interc(n) = beta n^alpha
alpha_interc = -9.58831961e-01
beta_interc = 1.26195529e-02
#-- VELOCITY PARAMETERS -------------------------------------------------------------





