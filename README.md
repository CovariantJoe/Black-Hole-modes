# Black-Hole-modes

In very, very short, Post-merger (collision) Black-hole radiation is described by complex frequencies known as Quasi Normal Modes (QNM). The real part describes oscillation while the imaginary corresponds to damping.

In the context of Schwarzschild post-mergers (non-rotating, non-charged black holes) QNMs are usually computed by either the WKB or Leaver's method.

Main.F90: Fortran90 code to find the QNMs by either method (WKB is 2nd order accurate) (Own implementations). Includes a program to call the solver of the Zerilly eq. (used to evolve the perturbation once the QNM is known).

ZerilliSolver.F90: RK45 solver (own implementation) for said equation.

Configs.nml: Parameter configuration for QNMs and the solver.

QNM.pdf: Brief article-like comparison of the results with known literature.
