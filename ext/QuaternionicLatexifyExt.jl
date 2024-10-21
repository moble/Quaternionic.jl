module QuaternionicLatexifyExt

import Quaternionic: AbstractQuaternion, QuatVec
isdefined(Base, :get_extension) ? (using Latexify) : (using ..Latexify)

function _pm_latex(x::Number)
    # Utility function to print a component of a quaternion in LaTeX
    s = Latexify.latexify(x, env=:raw, bracket=true)
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

function Base.show(io::IO, ::MIME"text/latex", q::AbstractQuaternion)
    s = Latexify.LaTeXStrings.latexstring(
        q isa QuatVec ? "" : Latexify.latexify(q[1], env=:raw, bracket=true),
        _pm_latex_latexify(q[2]), "\\,\\mathbf{i}",
        _pm_latex_latexify(q[3]), "\\,\\mathbf{j}",
        _pm_latex_latexify(q[4]), "\\,\\mathbf{k}"
    )
    print(io, s)
end

end  # module
