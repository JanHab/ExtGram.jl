#!/usr/bin/env julia
#
# Write the quadrature angle sets of src/angles.jl as LaTeX tables, several (θ_k, φ_k)
# pairs per row, in two styles:
#
#   angles_table.tex        θ_k and φ_k in separate, decimal-aligned siunitx S columns
#   angles_table_pairs.tex  every pair printed in brackets as (θ_k, φ_k)
#
#     julia miscellaneous/angles_to_tex.jl            # writes both files into miscellaneous/
#     julia miscellaneous/angles_to_tex.jl outdir     # ... into outdir
#
# Each file is a standalone document (booktabs + siunitx + longtable) with one table per
# dictionary: slab sets and full 3V sets. A table with more than LONGTABLE_ABOVE rows is a
# longtable (breaks across pages, each M block kept together); shorter ones are floats.
# For a paper, copy the part between the BEGIN TABLES / END TABLES markers, or \input it.

using Printf
include(joinpath(@__DIR__, "..", "src", "angles.jl"))   # SlabAngles, FullAngles

const OUTDIR       = isempty(ARGS) ? (@__DIR__) : ARGS[1]
const STYLES       = [:columns => "angles_table.tex", :pairs => "angles_table_pairs.tex"]
const PER_ROW      = 3          # (θ, φ) pairs per table row
const SHOW_INDEX   = true       # k column in front of every pair
const INDEX_BASE   = 0          # first index shown: (θ_0, φ_0), (θ_1, φ_1), ...
const GAP          = raw"@{\hspace{2em}}"               # separator between pair groups
const SUFFIX_ORDER = ["Max", "Arc", "Con", "Bad", "MinCond"]  # block order inside a table
# Text printed for a key suffix in the block header; unknown suffixes appear verbatim.
const LABELS = Dict("Max" => "Max", "Arc" => "Arc", "Con" => "Con", "Bad" => "Bad",
                    "MinCond" => "MinCond")
const LONGTABLE_ABOVE = 45      # tables with more rows than this become a longtable
const EOL = " \\\\"             # LaTeX row terminator

const SIGFIGS    = 6            # significant figures printed (the source values have at most 6)
const ZERO_BELOW = 1e-300       # |x| below this prints as 0: subnormal leftovers such as
                                # 3.23939e-315 are zero to machine precision and would otherwise
                                # print as 3.23939 x 10^-315. Set to 0 to print them verbatim.
const FMT = Printf.Format("%.$(SIGFIGS)g")
numstr(x) = Printf.format(FMT, abs(x) < ZERO_BELOW ? 0.0 : float(x))   # 3.14159, 1.63196e-09

# The tables to emit, in order: (dictionary, caption, label). Empty dictionaries are skipped.
const TABLES = [
    (SlabAngles,
     raw"Slab quadrature directions $(\theta_k, \phi_k)$ of \texttt{SlabAngles}, in radians: " *
     raw"$\theta_k \in [0, 2\pi)$ is the azimuth, $\phi_k \in [0, \pi]$ the polar angle " *
     raw"and $M$ the moment order.",
     "tab:slab-angles"),
    (FullAngles,
     raw"Full 3V quadrature directions $(\theta_k, \phi_k)$ of \texttt{FullAngles}; " *
     raw"conventions as in Table~\ref{tab:slab-angles}.",
     "tab:full-angles"),
]

# ---- the two styles ---------------------------------------------------------------------
struct Style
    colspec::String     # column spec of one group (index column + pair)
    header::String      # header cells of one group
    ncols::Int          # number of columns in one group
    cell::Function      # (θ, φ) -> the pair's cell text (may contain '&')
end

# A siunitx S column drops the exponent of an entry such as 1.63196e-09 when table-format
# has no exponent part, so those go in as braced \num text.
numcell(x) = (s = numstr(x); occursin('e', s) ? "{\\num{$s}}" : s)

function style(kind)
    idx_spec, idx_head, idx_n = SHOW_INDEX ? ("r ", raw"$k$ & ", 1) : ("", "", 0)
    if kind == :columns
        return Style(idx_spec * "S S", idx_head * raw"{$\theta_k$} & {$\phi_k$}", idx_n + 2,
                     (θ, φ) -> numcell(θ) * " & " * numcell(φ))
    elseif kind == :pairs
        return Style(idx_spec * "l", idx_head * raw"$(\theta_k, \phi_k)$", idx_n + 1,
                     (θ, φ) -> "\$(\\num{$(numstr(θ))},\\,\\num{$(numstr(φ))})\$")
    end
    error("unknown style $kind")
end

# ---- table emission ---------------------------------------------------------------------
"Split a key like \"M6_Arc\" into (6, \"Arc\")."
function parse_key(key)
    m = match(r"^M(\d+)_(.+)$", key)
    m === nothing && error("key $key is not of the form M<order>_<suffix>")
    return parse(Int, m[1]), String(m[2])
end

"Keys of `d` sorted by suffix (SUFFIX_ORDER first, unknown suffixes after) and then by M."
function ordered_keys(d)
    rank(s) = something(findfirst(==(s), SUFFIX_ORDER), length(SUFFIX_ORDER) + 1)
    sortkey(k) = (M = parse_key(k); (rank(M[2]), M[2], M[1]))
    return sort(collect(keys(d)); by = sortkey)
end

nrows(d)    = sum(1 + cld(length(v), PER_ROW) for v in values(d))   # block header + data rows
colspec(st) = "@{} " * join(fill(st.colspec, PER_ROW), " $GAP ") * " @{}"
header(st)  = join(fill(st.header, PER_ROW), " & ") * EOL

# `keep` = true keeps the block on one page (longtable honours `\\*` as "no page break here").
function emit_block(io, indent, st, M, suffix, pairs; keep)
    label  = get(LABELS, suffix, suffix)
    nlines = cld(length(pairs), PER_ROW)
    term(last) = keep && !last ? EOL * "*" : EOL
    println(io, indent, raw"\multicolumn{", PER_ROW * st.ncols, raw"}{@{}l}{$M = ", M, raw"$, ", label, "}", term(false))
    for (line, first) in enumerate(1:PER_ROW:length(pairs))
        groups = String[]
        for i in first:(first + PER_ROW - 1)
            if i <= length(pairs)
                θ, φ = pairs[i]
                push!(groups, (SHOW_INDEX ? string(i - 1 + INDEX_BASE, " & ") : "") * st.cell(θ, φ))
            else                                    # pad an incomplete last row
                push!(groups, join(fill("", st.ncols), " & "))
            end
        end
        println(io, indent, join(groups, " & "), term(line == nlines))
    end
end

function emit_blocks(io, indent, st, d; keep)
    previous_suffix = nothing
    for key in ordered_keys(d)
        M, suffix = parse_key(key)
        if previous_suffix !== nothing   # thin space within a strategy, rule between strategies
            println(io, indent, suffix == previous_suffix ? raw"\addlinespace" : raw"\midrule")
        end
        emit_block(io, indent, st, M, suffix, d[key]; keep)
        previous_suffix = suffix
    end
end

function emit_table(io, st, d; caption, label)
    n = PER_ROW * st.ncols
    if nrows(d) > LONGTABLE_ABOVE
        println(io, raw"\begin{longtable}{", colspec(st), "}")
        println(io, raw"  \caption{", caption, raw"}\label{", label, "}", EOL)
        println(io, raw"  \toprule")
        println(io, "  ", header(st))
        println(io, raw"  \midrule")
        println(io, raw"  \endfirsthead")
        println(io, raw"  \multicolumn{", n, raw"}{@{}l}{\tablename~\thetable{} (continued)}", EOL)
        println(io, raw"  \toprule")
        println(io, "  ", header(st))
        println(io, raw"  \midrule")
        println(io, raw"  \endhead")
        println(io, raw"  \midrule")
        println(io, raw"  \multicolumn{", n, raw"}{@{}r}{\footnotesize continued on next page}", EOL)
        println(io, raw"  \endfoot")
        println(io, raw"  \bottomrule")
        println(io, raw"  \endlastfoot")
        emit_blocks(io, "  ", st, d; keep = true)
        println(io, raw"\end{longtable}")
    else
        println(io, raw"\begin{table}[htbp]")
        println(io, raw"  \centering")
        println(io, raw"  \caption{", caption, "}")
        println(io, raw"  \label{", label, "}")
        println(io, raw"  \begin{tabular}{", colspec(st), "}")
        println(io, raw"    \toprule")
        println(io, "    ", header(st))
        println(io, raw"    \midrule")
        emit_blocks(io, "    ", st, d; keep = false)
        println(io, raw"    \bottomrule")
        println(io, raw"  \end{tabular}")
        println(io, raw"\end{table}")
    end
end

const PREAMBLE = raw"""
% Generated by miscellaneous/angles_to_tex.jl from src/angles.jl -- regenerate, do not edit.
\documentclass[10pt]{article}
\usepackage[a4paper, margin=2cm]{geometry}
\usepackage{booktabs}
\usepackage{longtable}
\setlength{\LTcapwidth}{\textwidth}   % longtable captions as wide as float captions
\usepackage{siunitx}
\sisetup{
  table-format = 1.6,   % S columns: one integer digit, up to six decimals (wider entries widen the column)
  group-digits = none,  % print 3.14159, not 3.141 59
  % Uncomment for a fixed number of decimals (1.63196e-9 then prints as 0.0000):
  % round-mode = places, round-precision = 4, exponent-mode = fixed, fixed-exponent = 0,
}
\begin{document}
"""

for (kind, filename) in STYLES
    out = joinpath(OUTDIR, filename)
    st = style(kind)
    open(out, "w") do io
        print(io, PREAMBLE)
        println(io, "% --- BEGIN TABLES (style: $kind) ---")
        first = true
        for (d, caption, label) in TABLES
            isempty(d) && continue
            first || println(io)
            emit_table(io, st, d; caption, label)
            first = false
        end
        println(io, "% --- END TABLES ---")
        println(io, raw"\end{document}")
    end
    println("wrote $out")
end
