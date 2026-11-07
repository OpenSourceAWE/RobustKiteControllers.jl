# Copyright (c) 2026 Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

"""
Compare the masses of the V3Kite model with those of AWETrim's quasi-steady model and
print them side by side: canopy, bridle, KCU, the airborne total and the tether mass per
metre.

V3Kite side: the `SystemStructure` of `project` (default: the project selected in
`data/gui.yaml`), built as `init` builds it before settling, without compiling a model.
The canopy is the `extra_mass` spread over the wing points, the KCU that of point 1. The
bridle is its line segments (density × cross-section × `l0`, as SymbolicAWEModels
computes it) plus the `extra_mass` of the remaining points (pulleys).

AWETrim side: the system file the runs send, `awetrim_system_config()` (AWETrim's own
file with the bridle mass set to `AWETRIM_BRIDLE_MASS`; `system_config` overrides it),
read with AWETrim's own Python in `awetrim_dir`. `canopy`, `bridle` and `kcu` are the
file's components; `wing in trim` is the point mass AWETrim's trim actually uses, which
includes the bridle since AWETrim `aa25dfc` (2026-10-29) and is the canopy alone before.

    include("examples/validate_masses.jl")
    run_example("validate_masses.jl"; project = "system_reelout_cabauw.yaml")
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(joinpath(@__DIR__))
end

using RobustKiteControllers: run_example, script_inputs, project_file, read_gui_field,
                             awetrim_system_config, AWETRIM_DIR
using V3Kite
using V3Kite: V3SettleConfig, build_settling_struct, struc_geometry_path,
              aero_geometry_path, load_kite, project_data_path, resolve_aero_mode
import SymbolicAWEModels
import JSON3
using Printf

(; project, awetrim_dir, system_config) = script_inputs(@__FILE__,
    (; project = something(read_gui_field("project"), "system_reelout_maasvlakte.yaml"),
       awetrim_dir = AWETRIM_DIR, system_config = awetrim_system_config()))

"Python run in AWETrim's venv: the masses of a system file as JSON on the last line"
const AWETRIM_MASSES_PY = """
import json, sys, warnings
warnings.simplefilter("ignore")
from awetrim.utils.config_paths import LEI_V3_SYSTEM_CONFIG
from awetrim.utils.system_config import get_kite, get_kcu_mass
from awetrim.system import factory
path = sys.argv[1] if len(sys.argv) > 1 else str(LEI_V3_SYSTEM_CONFIG)
cfg = factory._load_yaml(factory.Path(path))
kite = get_kite(cfg)
bridle = ((kite.get("bridle") or {}).get("structure") or {}).get("mass", 0.0) or 0.0
mass_wing, _, mass_kcu, _, _, d, rho, _, _ = factory._extract_params_awesio(cfg)
print(json.dumps(dict(system=path, canopy=kite["wing"]["structure"]["mass"],
    bridle=float(bridle), kcu=float(get_kcu_mass(cfg) or 0.0), wing_in_trim=mass_wing,
    kcu_in_trim=mass_kcu, tether_diameter=d, tether_density=rho)))
"""

"""
    v3kite_masses(project) -> NamedTuple

Masses [kg] of V3Kite's structure for `project`: `canopy`, `bridle_lines`,
`bridle_points`, `kcu`, and the tether's `tether_per_m` [kg/m].
"""
function v3kite_masses(project)
    path = project_file(project)
    data_path = project_data_path(path, nothing)
    kite_set = load_kite(path; data_path)
    config = V3SettleConfig(project = path, kite_set = kite_set, v_wind = 7.0,
                            tether_length = 150.0)
    sys, _ = build_settling_struct(config; data_path,
        source_struc = struc_geometry_path(path; data_path),
        source_aero = aero_geometry_path(path; data_path,
                                         aero_mode = resolve_aero_mode(kite_set)))
    seg_mass(s) = 2 * SymbolicAWEModels.segment_half_mass(s.l0, s.diameter, s.density)
    tether_segs = Set(Iterators.flatten(t.segment_idxs for t in sys.tethers))
    tether = [s for s in sys.segments if s.idx in tether_segs]
    wing_pts = Set(Iterators.flatten(st.point_idxs for st in sys.stations))
    return (canopy = sum(sys.points[i].extra_mass for i in wing_pts),
            bridle_lines = sum(seg_mass(s) for s in sys.segments if !(s.idx in tether_segs)),
            bridle_points = sum(p.extra_mass for (i, p) in enumerate(sys.points)
                                if i != 1 && !(i in wing_pts); init = 0.0),
            kcu = sys.points[1].extra_mass,
            tether_per_m = sum(seg_mass, tether) / sum(s.l0 for s in tether))
end

"""
    awetrim_masses(awetrim_dir, system_config) -> JSON3.Object

Masses of AWETrim's system file, from `AWETRIM_MASSES_PY` run in `awetrim_dir/venv`.
"""
function awetrim_masses(awetrim_dir, system_config)
    python = joinpath(awetrim_dir, "venv", "bin", "python")
    isfile(python) || error("No AWETrim venv at $python; run bin/install_awetrim.")
    args = isnothing(system_config) ? String[] : [system_config]
    out = read(Cmd(`$python -c $AWETRIM_MASSES_PY $args`; dir = awetrim_dir), String)
    return JSON3.read(last(split(strip(out), '\n')))
end

v3 = v3kite_masses(project)
at = awetrim_masses(awetrim_dir, system_config)
v3_bridle = v3.bridle_lines + v3.bridle_points
v3_above = v3.canopy + v3_bridle
at_tether_per_m = at.tether_density * π * at.tether_diameter^2 / 4

println("\nV3Kite:  $project")
println("AWETrim: $(at.system)\n")
@printf("%-26s %10s %10s %10s\n", "mass [kg]", "V3Kite", "AWETrim", "diff")
row(name, a, b) = @printf("%-26s %10.3f %10.3f %+10.3f\n", name, a, b, b - a)
row("canopy", v3.canopy, at.canopy)
@printf("%-26s %10.3f\n", "  bridle lines", v3.bridle_lines)
@printf("%-26s %10.3f\n", "  bridle point masses", v3.bridle_points)
row("bridle", v3_bridle, at.bridle)
row("KCU", v3.kcu, at.kcu)
row("airborne (file)", v3_above + v3.kcu, at.canopy + at.bridle + at.kcu)
row("wing in trim", v3_above, at.wing_in_trim)
row("airborne in trim", v3_above + v3.kcu, at.wing_in_trim + at.kcu_in_trim)
row("tether [kg/m]", v3.tether_per_m, at_tether_per_m)
nothing
