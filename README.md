# bouncesolver

A Fortran library for computing bounce solutions in scalar field theories using the
**Path Deformation Method (PDM)**. The bounce is the O(3)- or O(4)-symmetric instanton
that connects a false vacuum to a true vacuum in field space. Its Euclidean action S
determines quantum and thermal tunneling rates.

## Physics background

The bounce satisfies the equation of motion

$$\frac{d^2\phi}{d\rho^2} + \frac{\alpha}{\rho}\frac{d\phi}{d\rho} = \frac{\partial V}{\partial \phi}$$

with boundary conditions $d\phi/d\rho|_{\rho=0} = 0$ and $\phi(\rho\to\infty) = \phi_\text{false}$.
Here $\alpha = 2$ gives the O(3)-symmetric solution (thermal tunneling in 3D) and
$\alpha = 3$ gives the O(4)-symmetric solution (quantum tunneling in 4D).

For multi-field potentials the library maintains a 1D path through field space and
iteratively deforms it until the path aligns with the true bounce trajectory. Along
the current path the bounce equation reduces to a 1D problem that is solved by the
**overshoot/undershoot (bisection) method**: the initial field value at $\rho=0$ is
bisected until the solution reaches $\phi_\text{false}$ at large $\rho$.

## Requirements

- A Fortran compiler supported by [FPM](https://fpm.fortran-lang.org/) (gfortran ≥ 10
  or Intel ifx recommended).
- FPM ≥ 0.10.

All library dependencies are fetched automatically by FPM from their git repositories.

## Build

```bash
fpm build                    # default build
fpm build --profile release  # optimised: -O3, -march=native
fpm build --profile debug    # bounds-checks, backtracing, -Wall, -fcheck=all
```

Run the test suite (always use the debug profile for bounds-checks):

```bash
fpm test --profile debug                      # all tests
fpm test pdm_thickwalled --profile debug      # one specific test
```

Quad precision is opt-in at compile time:

```bash
fpm build --flag "-DQUAD"
```

## Using bouncesolver in your FPM project

Add the library as a git dependency in your `fpm.toml`:

```toml
[dependencies]
bouncesolver = { git = "https://gitlab.com/thomas.biekoetter/bouncesolver.git", profile = "release" }
```

Then use the public module in your source:

```fortran
use bouncesolver__pdm,    only : solver
use bouncesolver__config, only : wp        ! working precision kind
```

## User interface

The entire public API is the derived type `solver` in `bouncesolver__pdm`. Constructing
it runs the full algorithm and returns the result.

### Constructor

```fortran
pd = solver(d, V, phi_false, phi_true [, options])
```

| Argument | Type | Description |
|---|---|---|
| `d` | `integer` | Number of scalar fields |
| `V` | `procedure(real(wp) function(real(wp)(:)))` | Potential energy function |
| `phi_false` | `real(wp)(d)` | Location of the false vacuum |
| `phi_true` | `real(wp)(d)` | Location of the true vacuum |

All remaining arguments are **optional keyword arguments**:

| Keyword | Default | Description |
|---|---|---|
| `alpha` | `2` | Symmetry: `2` = O(3) (thermal), `3` = O(4) (quantum) |
| `n_odeint` | `2000` | Number of integration points along ρ |
| `maxiter` | `100` | Maximum number of path-deformation iterations |
| `deform_eps` | `0.02` | Step size for path deformation |
| `Rforce_threshold` | `0.1` | Convergence criterion: ratio \|Nforce\|/\|Pforce\| |
| `rho_max_fac` | `20.0` | Sets ρ\_max = rho\_max\_fac / m\_false |
| `smoothing` | `.false.` | Apply viscosity smoothing to path after deformation |
| `init_path_mode` | `1` | `1` = straight-line start; `2` = follow potential well |
| `barrier_buffer` | `1e-3` | Tolerance for locating the potential barrier top |
| `verbose_level` | `1` | `0` = silent, `1` = iterations, `2` = shooting details |

### Checking the result

After construction inspect `pd%exit_status`:

- `exit_status == 0`: converged successfully.
- `exit_status > 0`: reached maximum iterations (result is the best-so-far solution).
- `exit_status < 0`: a fatal error occurred.

Call `pd%print_exit_status()` to print a human-readable message.

### Key output fields

| Field | Type | Description |
|---|---|---|
| `pd%S` | `real(wp)` | Euclidean action of the best solution |
| `pd%Sp` | `real(wp)` | Potential part of the action |
| `pd%Sk` | `real(wp)` | Kinetic part of the action |
| `pd%Rforce_best` | `real(wp)` | Best force ratio \|Nforce\|/\|Pforce\| achieved |
| `pd%phi(n,d)` | `real(wp)` | Field values along the bounce path |
| `pd%rho(n)` | `real(wp)` | Radial coordinate ρ grid |
| `pd%xb(n)` | `real(wp)` | Arc-length coordinate x(ρ) along the bounce |
| `pd%pot(n)` | `real(wp)` | Potential V evaluated on the bounce |

### CSV export (optional)

The module `bouncesolver__export` provides convenience routines that write the bounce
data to CSV files:

```fortran
use bouncesolver__export, only : csv_x_of_rho, csv_pot_of_rho, csv_pot_of_x
use bouncesolver__export, only : csv_phi_of_rho, csv_phi_of_x, csv_forces_of_x

call csv_x_of_rho(pd,    "path/to/x_of_rho.csv")    ! rho, x, xdot
call csv_pot_of_rho(pd,  "path/to/pot_of_rho.csv")   ! rho, V
call csv_pot_of_x(pd,    "path/to/pot_of_x.csv")     ! x, V
call csv_phi_of_rho(pd,  "path/to/phi_of_rho.csv")   ! rho, phi_1, phi_2, ...
call csv_phi_of_x(pd,    "path/to/phi_of_x.csv")     ! x, phi_1, phi_2, ...
call csv_forces_of_x(pd, "path/to/forces_of_x.csv")  ! x, Nforce, Pforce components
```

## Example: thick-walled bounce in a 2D potential

This example is the test program `test/pdm_thickwalled.f90`, and is based on an example
from the CosmoTransitions manual ([arXiv:1109.4189](https://arxiv.org/abs/1109.4189)).
It computes the bounce for the two-field potential

$$V(\phi_1,\phi_2) = (\phi_1^2+\phi_2^2)\bigl[a(\phi_1-1)^2 + b(\phi_2-1)^2 - c\bigr]$$

with parameters $a=1.8$, $b=0.2$, $c=0.1$, which has a false vacuum near the origin
and a true vacuum near $(1,1)$.

```fortran
program example

  use bouncesolver__config, only : wp
  use bouncesolver__pdm,    only : solver
  use gradmin__descent,     only : minimize   ! optional: refine vacuum locations

  implicit none

  real(wp), parameter :: a = 1.8e0_wp
  real(wp), parameter :: b = 0.2e0_wp
  real(wp), parameter :: c = 0.1e0_wp

  integer,  parameter :: d = 2
  real(wp)            :: phi_false(d) = [0.0e0_wp, 0.0e0_wp]
  real(wp)            :: phi_true(d)  = [1.0e0_wp, 1.0e0_wp]
  type(solver)        :: pd

  real(wp) :: phi_min(d), V_min

  ! Refine vacuum locations with gradient descent
  call minimize(V, phi_true,  phi_min, V_min, maxiter=10000, mode=1)
  phi_true  = phi_min
  call minimize(V, phi_false, phi_min, V_min, maxiter=10000, mode=1)
  phi_false = phi_min

  ! Solve for the bounce (O(3) symmetry, thermal tunneling)
  pd = solver(          &
    d, V,               &
    phi_false, phi_true,&
    alpha          = 2, &
    n_odeint       = 2000, &
    verbose_level  = 1, &
    smoothing      = .true., &
    rho_max_fac    = 40.0e0_wp, &
    init_path_mode = 2, &
    maxiter        = 40)

  if (pd%exit_status >= 0) then
    write(*,*) "Action S =", pd%S
  else
    call pd%print_exit_status()
  end if

contains

  function V(phi) result(y)
    real(wp), intent(in) :: phi(:)
    real(wp) :: y
    y = (phi(1)**2 + phi(2)**2) * (a*(phi(1)-1.0e0_wp)**2 + b*(phi(2)-1.0e0_wp)**2 - c)
  end function V

end program example
```

### Interpreting the output

The solver prints a summary at the end of each run:

```
========================================
Result:
       Action =    <S>    <Sp>    <Sk>
  Force ratio =    <Rforce>
   Iteration  =  <best> / <total>
========================================
```

A force ratio well below 0.1 indicates a well-converged result. If `exit_status == 1`
the maximum iteration count was reached — the best solution found is still returned and
may be acceptable depending on the force ratio achieved.

## Potential gotchas

- **Locate the minima precisely.** The solver checks that `phi_true` and `phi_false` are
  distinct and uses the curvature at both minima to estimate the integration range. Pass
  the gradient-descent-refined locations rather than rough guesses.
- **`rho_max_fac`** controls how far out in ρ the integration runs. The default (20) is
  usually fine for thick walls; thin-wall transitions may need a larger value (40–100).
- **`init_path_mode = 2`** initialises the path by following the bottom of the potential
  in all fields except the first, which can help convergence for curved tunnelling paths.
- **Quad precision** is available for high-precision requirements by compiling with
  `-DQUAD` (see `src/config.F90`).

## License

GNU Affero General Public License v3.0 (AGPL-3.0). See `LICENSE`.
