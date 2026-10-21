const _PARAMETER_FIELDS = ("value", "units", "source", "location", "notes", "confidence")
const _CONFIDENCE_LEVELS = ("high", "medium", "low", "placeholder")
const _ProvenanceRow = NamedTuple{
    (:name, :units, :source, :location, :notes, :confidence),
    Tuple{Symbol, String, Union{Missing, String}, Union{Missing, String}, String, String},
}

function _parameter_error(path, name, message)
    throw(ArgumentError("Parameter $(repr(name)) in $(repr(path)): $message"))
end

function _parameter_text(entry, field, path, name; required=true)
    if !haskey(entry, field)
        required && _parameter_error(path, name, "missing required field `$field`")
        return missing
    end
    text = entry[field]
    text isa AbstractString ||
        _parameter_error(path, name, "`$field` must be a string")
    isempty(strip(text)) &&
        _parameter_error(path, name, "`$field` must not be blank")
    return text
end

"""
    load_params(path::AbstractString) -> (values=NamedTuple, provenance=Vector)

Read one top-level TOML table per scalar parameter. Return a `NamedTuple` of
values and a vector of named provenance rows, both sorted by parameter name.
Rows have columns `name`, `units`, `source`, `location`, `notes`, and `confidence`.

Every entry requires nonblank `units` and a `confidence` of `"high"`, `"medium"`,
`"low"`, or `"placeholder"`. Non-placeholder entries also require a finite
numeric `value` and nonblank `source` and `location` strings. Boolean values are
not accepted as numbers. `notes` is an optional string, defaulting to `""`.

Placeholder entries may omit `source` and `location`, represented as `missing`
in provenance. For an unavailable value, omit `value` or write
`value = "missing"`, set `confidence = "placeholder"`, and include a `TODO` in
`notes`. Explicit numeric placeholders are preserved with their confidence tag.

Values retain their parsed numeric types and source units. This function checks
metadata, not the validity of citations, and does not perform unit conversions.
Unknown fields and invalid entries raise `ArgumentError` naming the file and
parameter. TOML syntax errors and file-read errors propagate from `TOML.parsefile`.
"""
function load_params(path::AbstractString)
    tables = TOML.parsefile(path)
    names = sort!(collect(keys(tables)))
    parameter_values = Any[]
    provenance = _ProvenanceRow[]

    for name in names
        isempty(strip(name)) && _parameter_error(path, name, "name must not be blank")
        entry = tables[name]
        entry isa AbstractDict ||
            _parameter_error(path, name, "expected a TOML table")
        for field in sort!(collect(keys(entry)))
            field in _PARAMETER_FIELDS ||
                _parameter_error(path, name, "unknown field $(repr(field))")
        end

        units = _parameter_text(entry, "units", path, name)
        confidence = _parameter_text(entry, "confidence", path, name)
        confidence in _CONFIDENCE_LEVELS ||
            _parameter_error(path, name, "`confidence` must be high, medium, low, or placeholder")
        placeholder = confidence == "placeholder"
        source = _parameter_text(entry, "source", path, name; required=!placeholder)
        location = _parameter_text(entry, "location", path, name; required=!placeholder)
        notes = get(entry, "notes", "")
        notes isa AbstractString || _parameter_error(path, name, "`notes` must be a string")

        value = get(entry, "value", missing)
        unavailable = ismissing(value) || (value isa AbstractString && value == "missing")
        if unavailable
            placeholder ||
                _parameter_error(path, name, "missing `value` requires `confidence = \"placeholder\"`")
            occursin(r"\bTODO\b", notes) ||
                _parameter_error(path, name, "missing `value` requires a TODO in `notes`")
            value = missing
        elseif !(value isa Real) || value isa Bool || !isfinite(value)
            _parameter_error(path, name, "`value` must be a finite number or a missing placeholder")
        end

        push!(parameter_values, value)
        push!(provenance, (name=Symbol(name), units=units, source=source,
                          location=location, notes=notes, confidence=confidence))
    end

    values = NamedTuple{Tuple(Symbol.(names))}(Tuple(parameter_values))
    return (; values, provenance)
end
