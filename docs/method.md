# Method and variable contract

The public input is a set of first-order free-wave components

```text
g, h, a, b, kx, ky, Ux, Uy
```

using the same first-nine-input layout as `mf12_spectral_coefficients`.
At evaluation time the components are deposited on the requested FFT grid as
the positive analytic first-order elevation spectrum. Every nonlinear field
is generated internally from that spectrum.

The Green--Laplace construction separates the parent-frequency decay from
the output-wavenumber response. A prescribed Gauss--Laguerre rule then turns
the continuous integral into a fixed sequence of Fourier multipliers,
inverse transforms, products and forward transforms. There are no explicit
parent-pair or parent-triple loops inside the released GL executors.

`eta` is the physical free-surface elevation. `psi` is the physical
free-surface potential

```text
psi = phi(x,z=eta,t).
```

The flat trace `Phi = phi(x,z=0,t)` is not interchangeable with `psi`.
Although flat-potential states occur inside the order-three forcing graph,
they are not exposed as public surface-potential components.

The returned total fields are the sum of the explicitly listed components in
`audit.included_components`. Version 0.1 contains the positive pure-sum
linear, second-order and third-order components only; it is not a silently
truncated representation of every signed sector.
