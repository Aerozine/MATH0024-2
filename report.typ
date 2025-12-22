//#set text(font: "Iosevka")
#set page(numbering: "1")
#set heading(numbering: "1.")
#set math.equation(numbering: "(1)", number-align: bottom)
#import "@preview/dashy-todo:0.1.3": todo
#import "@preview/zero:0.5.0": ztable
#import "@preview/theorion:0.4.1": *
#import cosmos.clouds: *
#show: show-theorion
#let theorem = theorem.with(fill: blue.lighten(85%))
#let theorem-box = theorem-box.with(fill: blue.lighten(85%))
/// Custom block style
#let theorem = theorem.with(radius: 5pt)
#let theorem-box = theorem-box.with(radius: 5pt)
#import "@preview/physica:0.9.7": *
#set document(
  title: [MATH0024 PDEs Homework 2],
)
#title()
Loïc Delbarre (S215072)
= Fourier Transform Verification of the Fundamental Solution for the 3D Wave Equation
The goal is to prove that 
$expval(Phi_t\,accent(phi,hat))=expval(accent(Phi,hat)_t\,phi)$ for all Schwartz functions $phi$ on $RR^3$ \
#let bxi = [$norm(bold(xi))$]
By definition, \
$ #expval[$Phi_t, accent(phi,hat)$]&=#expval[$Phi_t, integral_(RR^3) exp(-i bold(x) dot bold(xi)) phi(bold(xi)) d bold(xi) $] \
&= t/(4 pi) integral_(norm(y)=1) integral_(RR^3) exp(-i c t bold(y) dot bold(xi)) phi(bold(xi)) d bold(xi)d S_y \
&= integral_(RR^3) t/(4pi) integral_(norm(y)=1) exp(- i c t bold(y) dot bold(xi))d S_y phi(bold(xi))  d bold(xi)
\
&=integral_(RR^3) t/(4pi) integral_(norm(y)=1)
exp(-i c t y_3 norm(bold(xi))) d S_y phi(bold(xi))d bold(xi)
$
By replacing it with spherical coordinates the inner integral becomes
$$
$ integral_(norm(y)=1)
exp(-i c t y_3 norm(bold(xi))) d S_y &=
integral_0^(2pi) integral_0^(pi) exp(-i c t cos(theta) norm(bold(xi))) sin(theta) d theta d chi \

&=

2pi integral_0^(pi) exp(-i c t cos(theta) #bxi) sin(theta) d theta

$
By performing a substitution $p=cos(theta)$,where  $d p=-sin(theta)d theta$, we obtain
$
2pi integral_(-1)^(1) exp(-i c t p #bxi)  d p  &=
2pi [ frac(exp(-i c t p #bxi),- i c t p #bxi )]_(-1)^(1)
\
&=frac(2pi,-i c t #bxi) (e^(- i c t #bxi)- e^( i c t #bxi)) \
&= frac(4pi sin(c t #bxi),c t #bxi)
$

#quote-box[
#set text(fill: black)
Therefore
$
#expval[$Phi_t,accent(phi,hat)$]&= integral_(RR^3) t/(4pi) dot frac(4pi sin(c t #bxi),c t #bxi) phi(bold(xi)) d bold(xi) \
&= integral_(RR^3)   frac(sin(c t #bxi),c #bxi) phi(bold(xi)) d bold(xi)
&= #expval[$accent(Phi,hat)_t,phi$]
$
]
= Absolute Stability Analysis of the Forward Euler Method
This question will consider the forward Euler method for the initial-value problem
$
cases(
  frac(d u,d t)(t) = lambda u(t) #h(1cm) & "for" t>0 \,&  #h(1cm) lambda in CC,
  u(0)=1 & "at" t = 0
)
$
The exact solution is given by 
$
u(t)= e^(lambda t)
$
The forward euler method is given by
$
u_(n+1)= (1 + triangle t dot lambda)u_n
$

#theorem-box(title:"A time marching method is absolutely stable")[
 for a specific timestep $triangle t$ if its application to this particular IVP leads for this timestep to a numerical solution with the same asymptotic behaviour. \
$ u_n -> 0 "as" n->+infinity "when" Re(lambda)<0 $
]
\
For the forward Euler method with initial condition, 
$
abs(u_n)= abs(1+triangle t dot lambda)^n ->0 "as" n-> +infinity
$
This requires
$
abs(1+triangle t dot lambda)<1
$
The region of absolute stability ( ie where the condition is respected) is a disk in the complex plane centered in (-1,0) with a radius of 1.\
For $lambda=-5$, the absolute stability inequation gives us:
$
-1 <= 1-5 triangle t <= 1
$
Where 
- The upperbound is trivial for a positive $triangle t$
- the lower bound gives $triangle t<=0.4$

#figure(
  image("graph/Q2S.svg",width:80%),
  caption:[Representation of the forward Euler method and the representation of the absolute error associated to the method for $triangle t$ in the stability interval.],
)

#figure(
  image("graph/Q2U.svg",width:80%),
  caption:[Representation of the forward Euler method and the representation of the absolute error associated to the method for unstable and stable marginally $ triangle t$],
)
#quote-box[
#set text(fill: black)
The term $abs(1+ triangle t dot lambda)$ acts as the amplification factor.
In the case of small value, the error will decay.
In the case of $triangle t=0.4$, the scheme is marginally stable, the oscillation does not even grow or decay. 
After this limit value the amplification increases the amplitude of the oscillations.
The absolute error explodes.
]

= Spectral Analysis and Wave Propagation in Rectangular Waveguides
== Spectral solution
The spectral problem equation is 
$
- frac(partial ^2 phi.alt,partial  x^2) = lambda phi.alt
$
Three cases must be discussed based on the sign of $lambda$
- Case $lambda <0$
Let $lambda=-alpha^2$, with $alpha >0$ \
In this particular case the general solution is given by:
$
phi(x)= A cosh(alpha x ) + B sinh(alpha x)
$
By applying the boundary condition
$
cases(
phi(- L/2)= A cosh(- alpha L/2)+B sinh(- alpha L/2)=0,
phi(L/2)= A cosh(alpha L/2)+ B sinh(alpha L/2)=0
)
$
Knowing that $sinh$ is a odd function and that $cosh$ is even:
$ cases(
phi(- L/2)= A cosh(alpha L/2)-B sinh(alpha L/2)=0,
phi(L/2)= A cosh(alpha L/2)+ B sinh(alpha L/2)=0
)
$
By rearranging, this gives
$
cases(
  2 A cosh(alpha L/2)=0
,
  2 B sinh(alpha L/2)=0
)
$
That concludes that $A=B=0$ and in extenso $phi(x)=0$
- Case $lambda=0$
In this case the general solution is 
$
phi(x)=A x+b
$
By applying bondary condition;
$
cases( 
  phi(-L/2)= -A L/2 +B =0
,
  phi(L/2)= A L/2 +B =0
)
$
This situation also concludes that $A=B=0$ and then that $phi (x)=0$
- Case $lambda >0$

In this last case, the general solution is
$
phi.alt(x)= A sin(sqrt(lambda)x) + B cos(sqrt(lambda)x)
$
Using the bondary condition,
$
cases(phi.alt(-L/2)=0,phi.alt(L/2)=0)
$
There is two sub-cases emerging 
- $B=0$  and $sin(sqrt(lambda)L/2)=0$
This conditions requires that 
$
A sin(sqrt(lambda)L/2)=0 \
lambda_m=(frac(2 m pi, L))^2
$
With $m in NN_0^+$
The value m=0 is excluded as it gives the trivial solution, and negative values give the same eigenvalues as positive ones.

- $A=0$ and $cos(frac(sqrt(lambda)L,2))=0$
In this case 
$
lambda_k=(frac((2k+1)pi,L))^2
$
with $k in NN_0^+$
To satisfy both boundary conditions with the domain centered at the origin, we use a phase-shifted representation.The complete set of eigenfunctions are :
$
phi.alt_m (x)=sin(frac( m pi ,L)(x+L/2))
$

By normalizing //( i.e. $integral_(-L/2)^(L/2) phi.alt_m^2(x) d x = 1$)
$
integral^(L/2)_(-L/2) sin^2(frac(m pi, L)(x+L/2)) d x = L/2
$

#quote-box[
#set text(fill: black)
The solution to satisfy the bondary condition is then given by
$
phi.alt_m (x)=sqrt(2/L) sin(frac( m pi ,L)(x+L/2))
$
For $m in NN_0^+$
]

With the same analysis the y-direction problem gives with 
$
kappa_n= (frac(n pi, H))^2
$

#quote-box[
#set text(fill: black)
The normalized form is given by
$
psi_n (y)= sqrt(2/H) sin(frac(n pi, H)(y+H/2))
$
]
== Derivation in the sense of the theory of distribution
For all test function $phi in RR$
$
#expval[$frac(d,d z) G_(m n), phi$]&= -#expval[$G_(m n),frac(d,d z) phi$] \
&=- integral_(-infinity)^(+ infinity) G_(m n) (z) frac(d,d z) phi d z
&= - integral_(-infinity)^0 B_(m n) e^(-i k_(f,m n)z) frac(d,d z) phi d z
\
&&- integral^(+infinity)_0 B_(m n) e^(i k_(f,m n)z) frac(d,d z) phi d z
$

The integration by part returns
$
= - [B_(m n) e^(-i k_(f,m n) z) phi]_(-infinity)^0
+integral_(-infinity)^0 B_(m n)(-i k_(f,m n)) e^(-i k_(f,m n)z) phi d z
\
- [B_(m n) e^(i k_(f,m n) z)phi ]^(+infinity)_0
+integral^(+infinity)_0 B_(m n)(i k_(f,m n)) e^(i k_(f,m n)z) phi d z
$
Since $phi$ is a Schwartz function with compact support,it vanishes for $abs(x)$ sufficiently large,
$
= -B_(m n ) phi(0)
+B_(m n ) phi(0)

+integral_(-infinity)^0 B_(m n)(-i k_(f,m n)) e^(-i k_(f,m n)z) phi d z \

+integral^(+infinity)_0 B_(m n)(i k_(f,m n)) e^(i k_(f,m n)z) phi d z
$

#quote-box[
#set text(fill: black)
  The first order derivative is then 
$
#expval[$frac(d,d z) G_(m n),phi$]=
integral_(-infinity)^0 B_(m n)(-i k_(f,m n)) e^(-i k_(f,m n) z) phi d z 
+integral^(+infinity)_0 B_(m n)(i k_(f,m n)) e^(i k_(f,m n)z) phi d z
$
]
The second order derivative can be obtain in the same way 
$
#expval[$frac(d^2,d z^2) G_(m n), phi$]&= #expval[$G_(m n),frac(d^2,d z^2) phi$]  \
&= integral_(-infinity)^(+ infinity) G_(m n) (z) frac(d^2,d z^2) phi d z \
&=  integral_(-infinity)^0 B_(m n) e^(-i k_(f,m n)z) frac(d^2,d z^2) phi d z 
+ integral^(+infinity)_0 B_(m n) e^(i k_(f,m n)z) frac(d^2,d z^2) phi d z
$
 By decomposing the left integral (negative domain)
 $
 integral_(-infinity)^0 B_(m n) e^(- i k_(f,m n)z) frac(d^2,d z^2) phi d z&=
 [B_(m n)e^(- i k_(f,m n)z)frac(d,d z) phi]_(- infinity)^0 \
 &- integral^0_(-infinity) (-i k_(f,n m)B_(m n)e^(-i k_(f,m n)z))frac(d,d z)phi d z \
 =B_(m n) frac(d,d z) phi(0) + i k_(f,m n) B_(m n) phi(0) 
 & + integral_(-infinity)^0 (-k^2_(f,m n))B_(m n) e^(- i k_(f,n m)z)phi d z
 $
 By doing the same thing developpement for the right integral (positive domain)
$
 integral^(+infinity)_0 B_(m n) e^(i k_(f,m n)z) frac(d^2,d z^2) phi d z 
 = - B_(m n) frac(d, d z) phi(0) + i k_(f,m n) B_(m n) phi(0) \
 + integral_0^(+infinity) (-k^2_(f,m n))B_(m n) e^(i k_(f,m n)z) phi d z
$
By adding both contribution 

#quote-box[
#set text(fill: black)
  The second order derivative is then 
  $
#expval[$frac(d^2,d z^2) G_(m n),phi$]
= - k_(f,m n)^2 
#expval[$G_(m n),phi$]+2 i k_(f, m n) B_(m n) phi(0)
  $
]
== Evaluation of $B_(m n)$
By substituing in the equation
$
-c^2  [-k_(f,n m)^2 G_(m n)+ 2 i k_(f,n m)B_(m n) delta(z)]
+ ( - omega_f^2+ c^2(lambda_m + kappa_n))G_(m n)= phi.alt_m (0)psi_n (0)delta(z)
$

Knowing that $k^2_(f,m n)= (omega_f^2)/c^2-(lambda_m+kappa_n)$
$
-2i c^2 k_{f,m n}B_{m n} delta(z)= phi.alt_m (0)psi_n (0)delta(z)
$


#quote-box[
#set text(fill: black)
The value found is: 
$
B_(m n)= i frac(phi.alt_m (0)psi_n (0),2 c^2 k_(f, m n))
$
]
= Von Neuman Stability Analysis of the Upwind Finite Difference Scheme
== Numerical Implementation and Error Observation <4A>

The upwind method is defined as 
$
cases(u_j^(n+1)= u_j^n - c (triangle t)/h (u_j^n-u^n_(j-1)),
u_j^0=u_0(x_j)
)
$
Where $nu = (c triangle t)/h$ \
By assuming a solution on the form:
$
u_j=gamma^n e^(i xi j h)
$
The schema can be rewritten as
$ gamma ^(n+1) e^(i xi j h)
= gamma^n e^(i xi j h) - nu gamma^n(
  e^(i xi j h) - e^(i xi (j-1)h)
)
\

gamma = 1 - nu (1 - e^(-i xi h))
$
Since the transport equation propagates the initial profile unchanged at the velocity of $c=1$, the exact solution at time $t$ is given by
$
u(x,t)= sin(2 pi (x-t)).
$
This represents a rightward-propagating wave that retains its shape and amplitude, exhibiting no dispersion.
The comparison between the exact and numerical solutions in @Q4A reveals characteristic numerical artifacts.\
In Case 1, 
 the numerical solution exhibits moderate amplitude reduction and slight phase lag relative to the exact solution, indicating the presence of both dissipation and dispersion errors. The dissipation error manifests as artificial damping of the wave amplitude, while the dispersion error causes different Fourier components to propagate at slightly different velocities, resulting in the observed phase shift.
\
Case 2, demonstrates significantly more pronounced numerical errors: the wave amplitude is severely attenuated and the phase lag is substantially larger, reflecting the coarser spatial resolution and larger time step. These observations align with the theoretical understanding that numerical errors accumulate more rapidly when the grid spacing increases, as the discrete approximation of spatial derivatives becomes less accurate.

#figure(
  image("graph/Q4A.svg",width:100%),
  caption:[Comparison of exact and numerical solution for the transport equation at t≃1using the upwind finite difference method.]
)<Q4A>

== Spectral Analysis of Dissipation and Dispersion Errors <4B>
The Von Neumann analysis results from provide quantitative explanation for these observations. The amplification error plot ($e_a$) shows that $abs(gamma(xi)) < 1$ for all non-zero wavenumbers, confirming artificial dissipation across all Fourier modes, with higher wavenumbers experiencing greater damping.This directly accounts for the amplitude reduction observed in @4A. The dispersion error plot ($e_d$) reveals that $e_d (xi) < 1$ for $xi != 0$, indicating that numerical wave components propagate slower than the exact solution, which explains the phase lag. Crucially, both error measures deteriorate more severely in Case 2 due to the coarser discretization (larger $h$ means $xi h$ spans a wider range of the error spectrum), confirming that the increased spatial step size amplifies both dissipation and dispersion effects. 

#figure(
image("graph/Q4B.svg",width:100%),
  caption:[Von Neumann analysis of the upwind method, with dissipation $e_a$ and dispersion $e_d$ errors.],
)<Q4B>

= Use of AI
AI was used to summarise the instructions for question 3.
