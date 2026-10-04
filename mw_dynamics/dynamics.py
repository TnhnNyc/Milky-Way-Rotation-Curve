from __future__ import annotations
import numpy

def fit_polynomial_derivative(x, y, weights, degree, min_points=None,
                              eval_at=None):

    x = numpy.asarray(x, dtype=float)
    y = numpy.asarray(y, dtype=float)
    weights = numpy.asarray(weights, dtype=float)

    if eval_at is None:
        eval_at = x

    eval_at = numpy.asarray(eval_at, dtype=float)

    if min_points is None:
        min_points = degree + 2

    # Control whether there is Nan or inf or -inf expression; output: True/False 
    finite = numpy.isfinite(x) & numpy.isfinite(y) & numpy.isfinite(weights) & (weights > 0)

    if finite.sum()< min_points:
        NaN = numpy.full(eval_at.shape, numpy.nan, dtype=float)
        return NaN, NaN

    # Fit weighted polynomial and calculate analytics derivative
    coeffs = numpy.polyfit(x[finite], y[finite], deg=degree, w=weights[finite])
    poly_equation = numpy.poly1d(coeffs)
    d_poly_equation = poly_equation.deriv()

    print("---Fitting Polynomial---\n", poly_equation)
    print("\n---Analytics Derivative---", d_poly_equation)

    return poly_equation(eval_at), d_poly_equation(eval_at)
        
