# Exact one-loop check for the (p,q) string equation

`exact_pure.py p q` proves, in exact arithmetic over Q[w]/(f(w)) (one root of f per leading symbol), that the
hbar^2 coefficient of the pure branch solution of [P,(L^q)_+] = kappa equals p*gamma1 = p(p-1)(q-1)/(24(p+q-1))
on every leading symbol. It exits with status 0 only if every check passes. Dependencies: sympy, mpmath, python-flint;
Singular is used for the Groebner step when installed.
