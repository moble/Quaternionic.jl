# Functions of time

When using quaternions to represent rotations and orientations, we frequently
model dynamical systems, which means that the quaternions must be regarded as
functions of time.  The preceding functions can all be applied at each instant
of time, but we also need to be deal with the *change* of quaternions over
time, for which there are several important techniques:

  * [Interpolation](@ref) — Taking discretely sampled quaternionic time series and
    interpolating to different samples, and possibly differentiating.
  * [Angular velocity](@ref) — Both differentiation of a quaternionic function of time
    to get angular velocity, and integration of angular velocity to get an
    orientation as a function of time.
  * [Minimal rotation](@ref) — Finding the *least* dynamical motion that can achieve
    pointing in a certain direction.
  * [Derivatives and gradients](@ref) — Thanks to the [chain
    rule](https://en.wikipedia.org/wiki/Chain_rule#Multivariable_case),
    differentiating many nontrivial quaternionic functions of time will also
    involve differentiating with respect to components of quaternionic
    arguments.


## Interpolation

Component-wise interpolation of quaternions does not generally yield good
results when the quaternions are interpreted as rotations.  The basic reason is
that rotations correspond to *unit* quaternions, but component-wise
interpolation does not respect this constraint.  There are two specialized
functions for dealing with this problem.  The first is [`slerp`](@ref), which
is an abbreviation of "Spherical Linear intERPolation", and is the direct
analog of standard linear interpolation of functions ℝ → ℝ.  The second is
[`squad`](@ref), which is an abbreviation of "Spherical QUADrangle
interpolation", and is more analogous to cubic interpolation by Bézier splines.
The first is relatively fast but has discontinuous first derivatives at the
input points, while the second is somewhat slower but has continuous first and
second derivatives.

In both cases, it is important for extraneous sign flips to be eliminated
before passing quaternions to the interpolating functions.  For this purpose,
there is the [`unflip`](@ref) utility function, which can also be called
automatically by passing the corresponding keywords to `slerp` and `squad`.

As noted [below](#Gradients), `slerp` and `squad` can also be differentiated
analytically (or with `ForwardDiff`).

```@autodocs
Modules = [Quaternionic]
Pages   = ["interpolation.jl"]
```


## Angular velocity

Given a quaternionic function of time ``q(t)``, we can — in principle —
differentiate it to find ``\dot{q}(t)``.  This is related to the more familiar 
[angular velocity](https://en.wikipedia.org/wiki/Angular_velocity)
``\vec{\omega}`` by the equation
```math
\vec{\omega} = 2 \dot{q}\, q^{-1}.
```
While ``\dot{q}(t)`` can certainly be useful, it "lives" in a four-dimensional
space, which means that it takes a little more work to relate to our
three-dimensional space.  On the other hand, the direction of ``\vec{\omega}``
represents the instantaneous axis of rotation in our three-dimensional space,
and its magnitude describes the rate of rotation about that axis.  More
importantly, many results from physics describe rotations in terms of
``\vec{\omega}``, rather than ``\dot{q}(t)``.  So generally, we need a way to
go from ``q(t)`` to ``\vec{\omega}(t)`` — and back.

The details of this process and some examples are discussed in full in [this
paper](https://arxiv.org/abs/1604.08139).

Taking ``q(t)`` and obtaining ``\vec{\omega}(t)`` is fairly straightforward:
one simply needs to obtain ``\dot{q}(t)`` and apply the equation above.  If
``q(t)`` is known analytically, it may be possible to compute ``\dot{q}``
directly (as in the example from Sec. 6.2 of [the
paper](https://arxiv.org/abs/1604.08139)).  But ordinarily, this is a difficult
task.  For numerical functions automatic differentiation can be used to obtain
a numerical result for ``\dot{q}`` (see [below](#Gradient)).  Or, if ``q`` is
discretely sampled, [`∂squad∂t`](@ref) can be used to find the derivative of
the interpolant at any instant.

Going the other way, obtaining ``q(t)`` from ``\vec{\omega}(t)``, is more
delicate — though still possible.  It requires integrating the ordinary
differential equation
```math
\dot{q} = \frac{1}{2} \vec{\omega}\, q,
```
which is just a rearrangement of the equation above to solve for ``\dot{q}``.
Standard theorems on differential equations tell us that this equation has a
solution for ``q(t)`` as long as we know ``\vec{\omega}`` as a function of time
(and possibly of ``q``, but not ``\dot{q}``), and supply an initial value for
``q`` at some instant of time.  Note that [the
paper](https://arxiv.org/abs/1604.08139) noted above found that integrating
this equation without any restriction on the norm of ``q`` is generally the
most accurate and efficient choice.  However, with loose integration
tolerances, numerical error may result in the output ``q(t)`` having norm
significantly different from 1.  In this case, to use ``q`` as a rotation,
you can normalize the result, or simply apply it by "true" conjugation with the
inverse as in
```math
q\, \vec{v}\, q^{-1},
```
rather than using ``\bar{q}`` in the last factor.

If `DifferentialEquations.jl` is available, we can solve the ODE with code like
this:
```julia
using Quaternionic
using DifferentialEquations

# Make up some simple problem to solve
Ω⃗ = randn(QuatVecF64)
ω⃗(t) = Ω⃗
q₀ = randn(RotorF64)
tspan = (0.0, 100.0)

# Construct the ODE and `ODEProblem`
angular_velocity_ode(q, ω⃗, t) = ω⃗(t) * q / 2
angular_velocity_problem = ODEProblem(angular_velocity_ode, Quaternion(q₀), tspan, ω⃗)

# Now, solve it
q = solve(angular_velocity_problem, Vern9(), abstol=1e-12, reltol=0)
```
Here, we have constructed a very simple `ω⃗` function — which can actually be
integrated analytically — but it could be any callable that returns a `QuatVec`
(or `Quaternion`).  Note that the initial condition `q₀` has to be a
`Quaternion` on input to `ODEProblem` because `DifferentialEquations.jl` tries
to construct a zero element of that type, which is impossible for a `Rotor`.
This initial condition is also used for the output type and — as mentioned
above — numerical errors will tend to move the norm away from 1, so a `Rotor`
would not make sense anyway.  The choices of `Vern9` and `abstol=1e-12` here
are quite stringent, and may be overkill for many problems.  It is recommended
to set `reltol=0` in all cases.

The `q` that results from `solve` contains the time steps `q.t` at which the
solution was found along with the solution values `q.u` at each of those time
steps, but can also be used as a function of time to compute the value `q(t)`
at an arbitrary time in `tspan`.  In this simple case, we can compute the exact
value of, and compare it to the result of integration:
```julia
q_exact = @. exp([Ω⃗] * q.t / 2) * q₀
maximum(distance.(q.u, q_exact))
```
Typical numbers for this maximum error are roughly `1e-14`, though they
increase as `tspan` is increased.

Also note that this simple version involves allocations at each time step.  If
this is problematic, it should be easy to [define in-place
updating](https://diffeq.sciml.ai/dev/tutorials/ode_example/#Example-2:-Solving-Systems-of-Equations)
versions of `ω⃗` and `angular_velocity_ode` to eliminate allocations, using
vectors to hold the numbers instead of quaternions.


## Minimal rotation

One common problem arises when a system must be rotated to align some axis with
some direction, though the rotation of the system *about* that axis is
irrelevant.  To be specific, suppose we want to rotate our basis vectors
``\hat{x}, \hat{y}, \hat{z}`` so that ``\hat{z}`` points in a particular
direction.  A naive approach may be to determine the direction in terms of
spherical coordinates ``(\theta, \phi)``, and then use the rotation determined
by Euler angles ``(\alpha, \beta, \gamma) = (\phi, \theta, 0)``.  However, as
with most applications of Euler angles, this is a terrible idea.  The resulting
orientation will be extremely sensitive to the direction whenever it happens to
be near the poles.  In such cases, the angular velocity of the system will be
very high — potentially infinite, in principle.

Fortunately, it is possible to take *any* rotor ``R_\mathrm{axis}(t)`` that
aligns the axis correctly, and compute another rotation that also aligns the
axis, but has the smallest possible angular velocity.  This is called the
"minimal rotation".  The general problem is discussed (by way of a very
specific physical situation) in detail in Sec. III and Appendix B of [this
paper](https://arxiv.org/abs/1110.2965).  Essentially, we solve the equation
```math
\dot{\gamma} = 2 \left[\dot{R}_\mathrm{axis}\, 𝐤\, \bar{R}_\mathrm{axis}
\right]_w
```
(where the subscript ``w`` just takes the scalar component) for ``\gamma(t)``.
We then construct the rotor
```math
R(t) = R_\mathrm{axis}(t)\, \exp\left[\gamma(t) 𝐤 / 2 \right]
```
This rotor also aligns the axis correctly, but otherwise has the smallest
possible angular velocity.  Here, ``R_\mathrm{axis}`` may be constructed in any
convenient way, including using spherical coordinates; the resulting ``R(t)``
will be independent of such poor life choices.


## Derivatives and gradients

It can be very useful to compute the derivative of various functions with
respect to their arguments — for example, when computing angular velocity of a
`squad` interpolant, one needs to use the chain rule, and therefore needs each
of the derivatives of `exp`, `log`, etc.  Here, we provide the essential
(nontrivial) derivatives, treating each quaternionic argument as a series of
four real arguments.  For each input argument, the output is generally a
quaternion; interpreting those outputs as also being a series of four real
quantities, these derivatives could also be thought of as [Jacobian
matrices](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant) of the
relevant functions, though the actual return types are collections of
`Quaternion` objects.

Because of massive simplifications that result when using the right types,
these derivatives are more strict about input types than the main functions
themselves.  For example, the derivatives of `exp` are defined only for
`QuatVec` arguments; the derivatives of `log` are defined only for `Rotor`
arguments; etc.

Note that gradients can also be calculated automatically using
[`ForwardDiff.jl`](https://juliadiff.org/ForwardDiff.jl/).[^2] For example, we
could compute the derivative of `slerp` with respect to the `y` component of
the first input quaternion as
```julia
∂slerp∂q₁y(q₁, q₂, τ) = ForwardDiff.derivative(ϵ->slerp(q₁+ϵ*imy, q₂, τ), 0)
```
This is equal to
```julia
slerp∂slerp(q₁, q₂, τ)[2][3]
```
though `slerp∂slerp` computes the value and all derivatives of `slerp`
simultaneously, and does at least as fast as — and likely faster than — most AD
systems would.  Nonetheless, it is useful to know that `ForwardDiff` can
process functions *involving* `Quaternionic` methods.

[^2]:
    As of this writing, `Zygote.jl` does *not* work.  I'm not sure why, but I'm
    guessing it's related to [`Zygote`'s troubles with complex
    numbers](https://github.com/FluxML/Zygote.jl/issues/342).  I don't really
    understand AD terminology, but [this
    comment](https://discourse.julialang.org/t/automatic-differentiation-of-complex-valued-functions/30263/8)
    suggests forward-mode AD is better for this kind of thing anyway.  In any
    case, pull requests for improving this package's interaction with AD are
    certainly welcome.

```@autodocs
Modules = [Quaternionic]
Pages   = ["gradients.jl"]
```
