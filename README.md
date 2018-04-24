# IR
Simulator for IR spectrum of small molecules

Developed by Chelsea Sweet, Dina Moustafa, Kelly Moran, and Jonathan J. Foley.

Uses the finite-difference method in the time-domain to solve the time-dependent Schrodinger equation for a diatomic molecule within the Harmonic Oscillator model with dipolar coupling to a laser field. The shape and intensity of the laser field can be adjusted, and the response of the molecule via its time-dependent dipole moment and the related absorption spectra can be computed. A simple and relatively inefficient implementation of the discrete Fourier transform is used to compute the absorption spectrum from the time-dependent dipole moment.
