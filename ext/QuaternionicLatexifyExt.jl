module QuaternionicLatexifyExt

import Quaternionic: AbstractQuaternion, QuatVec
using Latexify

function latexraw_component(x::Number)
    # Utility function to print one component of a quaternion in raw LaTeX.  Latexify
    # writes NaN as plain text, which would be set in math italics.
    if x isa AbstractFloat && isnan(x)
        "\\mathrm{NaN}"
    else
        Latexify.latexify(x, env=:raw, bracket=true)
    end
end

function _pm_latex(x::Number)
    # Utility function to print a component of a quaternion in LaTeX
    s = latexraw_component(x)
    if s[1] ∉ "+-"
        s = "+" * s
    end
    if occursin(r"[+^/-]", s[2:end])
        if s[1] == '+'
            s = " + " * "\\left(" * s[2:end] * "\\right)"
        else
            s = " + " * "\\left(" * s * "\\right)"
        end
    else
        s = " " * s[1] * " " * s[2:end]
    end
    s
end

function latexraw_quaternion(q::AbstractQuaternion)
    # Build the LaTeX form of `q`, without math delimiters.  As in the plain-text form, the
    # scalar part of a `QuatVec` is omitted.
    string(
        q isa QuatVec ? "" : latexraw_component(q[1]),
        _pm_latex(q[2]), "\\,\\mathbf{i}",
        _pm_latex(q[3]), "\\,\\mathbf{j}",
        _pm_latex(q[4]), "\\,\\mathbf{k}"
    )
end

# This recipe makes `latexify(q)` agree with the `text/latex` form below, and lets
# Latexify's own methods for arrays and tables display arrays of quaternions.  Without the
# `env` setting, Latexify would return the `LaTeXString` without math delimiters.
Latexify.@latexrecipe function f(q::AbstractQuaternion)
    env --> :inline
    return Latexify.LaTeXStrings.LaTeXString(latexraw_quaternion(q))
end

function Base.show(io::IO, ::MIME"text/latex", q::AbstractQuaternion)
    print(io, Latexify.LaTeXStrings.latexstring(latexraw_quaternion(q)))
end

end  # module
