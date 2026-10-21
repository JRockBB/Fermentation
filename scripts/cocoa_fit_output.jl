# CSV checkpoint columns are explicitly ordered by name. The parameter loader
# and model use different NamedTuple orders, so positional serialization is unsafe.
using TOML

function cocoa_fit_initial(published, path)
    saved = TOML.parsefile(path)
    parameters = get(saved,"parameters",nothing)
    parameters isa AbstractDict && Set(keys(parameters)) == Set(String.(keys(published))) ||
        throw(ArgumentError("initial_file must contain all 24 kinetics and sigma in [parameters]"))
    return NamedTuple{keys(published)}(Tuple(parameters[String(name)] for name in keys(published)))
end

function write_cocoa_draw_header(io, names)
    println(io,"chain,draw,",join(String.(names),','),",divergent,tree_depth")
end

function write_cocoa_draw(io,chain,draw,parameters,diagnostics,names)
    println(io,join((chain,draw,(parameters[name] for name in names)...,
        diagnostics.divergent,diagnostics.tree_depth),','))
end
