# Copyright (c) 2026 Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

"""
Compare the lift and drag of the V3Kite model with AWETrim's reduced-order model (ROM)
at the same power-tape length and apparent wind speed, and save the result to
`output/lift_drag_validation.json` for `show_validation.jl`.

V3Kite side: the kite is parked at 150 m as in `simple_auto_parking.jl` (heading PID,
winch length loop) and flies one block per point of `points`, `(tape, v_a)` with the
power-tape length in m and the target apparent wind speed in m/s. In each block the
ground wind is first trimmed for `t_trim` seconds until `v_a` reaches the target (rate
limited, so the kite stays parked), then held for `t_avg` seconds, over which lift, wing
drag, `span_mean_aoa` and `v_a` are averaged. Lift is V3Kite's `lift_drag`: the wing's
force perpendicular to its apparent wind. Drag is put on AWETrim's basis, the VSM wing drag
plus the drag of ALL bridle lines (`drag_parts`); `lift_drag`'s own drag, which
holds only the bridle drag lumped onto the wing nodes, and the tether and KCU drag are
saved next to it.

AWETrim side: `POST /aero_forces` of the AWETrim server at `url` (needs AWETrim with that
endpoint; restart the server after switching branch), at the MEASURED mean `v_a` of each
block, for the ROM `AWETRIM_ROM` (`src/awetrim_client.jl`), at the ROM's own
gravity-free trim.

`fly = false` re-queries AWETrim for the V3Kite results already saved, without flying:

    include("examples/validate_lift_drag.jl")
    run_example("validate_lift_drag.jl"; fly = false)
    run_example("validate_lift_drag.jl"; points = [(1.6, 12.0)], url = "http://127.0.0.1:8001")
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(joinpath(@__DIR__))
end

using RobustKiteControllers
using RobustKiteControllers: run_example, script_inputs, package_git_state, AWETRIM_URL,
                             AWETRIM_ROM
using V3Kite
using V3Kite: init, step!
using DiscretePIDs: set_K!
using HTTP, JSON3
using OrderedCollections: OrderedDict
using LinearAlgebra: dot, norm, normalize
import SymbolicAWEModels
using Statistics
using Printf
import Dates
include(joinpath(@__DIR__, "winch_adapter.jl"))  # winch_torque!: V3Kite is torque-only

@info "validate_lift_drag.jl: lift and drag of the parked V3Kite vs. AWETrim's ROM at the \
       same power tape and apparent wind speed, saved for show_validation.jl."

"Where the results are saved"
const RESULT_FILE = normpath(joinpath(@__DIR__, "..", "output", "lift_drag_validation.json"))
"The ROMs queried: label => `rom_config_path`; only the one every optimizer request uses"
const ROMS = [splitext(AWETRIM_ROM)[1] => AWETRIM_ROM]
"The parked setup of `simple_auto_parking.jl`, so its settled state is reused"
const PARK = (project = "system_reelout.yaml", v_wind = 9.51, tether_length = 150.0,
              depower = 0.25, dt = 0.05 / 3, body_damping = [0.0, 0.0, 40.0],
              damping_per_stiffness = 0.002)
"Heading PID of `simple_auto_parking.jl`, gain scheduled with 1/v_a"
const HEADING = (P = 1.1, I = false, D = 0.15, v_app_ref = 13.1, v_app_min = 5.0,
                 max_steering = 0.175)
"Wind trim: gain [1/s], max rate [m/s²] and the time constant of the v_a filter [s]"
const WIND_TRIM = (gain = 0.4, max_rate = 0.5, tau = 2.0)

(; points, t_trim, t_avg, fly, url) = script_inputs(@__FILE__,
    (; points = [(1.58, 15.3), (1.58, 8.8), (1.76, 15.3), (1.76, 8.7)],
       t_trim = 30.0, t_avg = 10.0, fly = true, url = AWETRIM_URL))

"""
    check_aero_forces_endpoint(url)

Fail before flying when the AWETrim server at `url` does not answer or has no
`/aero_forces` (an older AWETrim, or a server started before the branch switch).
"""
function check_aero_forces_endpoint(url)
    spec = try
        String(HTTP.get(url * "/openapi.json"; request_timeout = 5).body)
    catch err
        error("No AWETrim server at $url ($(typeof(err))); start it with bin/run_server start.")
    end
    occursin("/aero_forces", spec) ||
        error("The AWETrim server at $url has no /aero_forces: check out AWETrim's branch \
               feat/server-aero-forces (until it is merged) and restart the server \
               (bin/run_server restart).")
    nothing
end

"""
    drag_parts(s) -> NamedTuple

Drag [N] along the wing's apparent wind, split as AWETrim needs it: `vsm` (the VSM wing),
`bridle` (all bridle lines), `bridle_in_wing` (the part of it `lift_drag` counts as wing
drag: SymbolicAWEModels lumps half of each segment's drag onto each end point, and
`lift_drag` adds the wing nodes' share), `tether` and `kcu` (the points with an area).
Segment drag is recomputed with SymbolicAWEModels' law at the mean wind and velocity of
the two end points, so it matches the model to about 2 %.
"""
function drag_parts(s)
    sys = s.sys
    wing = sys.wings[1]
    va_b = wing.va_b
    va_hat = V3Kite.calc_R_b_w(sys) * (va_b / norm(va_b))
    tether_segs = Set(Iterators.flatten(t.segment_idxs for t in sys.tethers))
    bridle = bridle_in_wing = tether = 0.0
    for seg in sys.segments
        a, b = sys.points[seg.point_idxs[1]], sys.points[seg.point_idxs[2]]
        a.is_wing_node && b.is_wing_node && continue  # wing frame: no line drag
        u = normalize(b.pos_w - a.pos_w)
        va = 0.5 * (a.wind_vec + b.wind_vec - a.vel_w - b.vel_w)
        rho = SymbolicAWEModels.air_density(s.sam.am, 0.5 * (a.pos_w[3] + b.pos_w[3]))
        f = 0.5 * rho * s.set.cd_tether * norm(va) * seg.len * seg.diameter *
            (va - dot(va, u) * u)
        d = dot(f, va_hat)
        if seg.idx in tether_segs
            tether += d
        else
            bridle += d
            bridle_in_wing += 0.5 * d * (a.is_wing_node + b.is_wing_node)
        end
    end
    kcu = sum(sys.points) do p
        p.area > 0 || return 0.0
        rho = SymbolicAWEModels.air_density(s.sam.am, p.pos_w[3])
        dot(SymbolicAWEModels.point_drag_force(p.wind_vec - p.vel_w, rho, p.drag_coeff, p.area),
            va_hat)
    end
    return (; vsm = dot(wing.aero_force_b, va_b / norm(va_b)), bridle, bridle_in_wing, tether, kcu)
end

"""
    fly_points(points; t_trim, t_avg) -> Vector{OrderedDict}

Fly one parked block per `(tape, v_a)` point in a single V3Kite run and return the
block means (forces in N, angles in degrees). `drag` is on AWETrim's basis (VSM wing +
all bridle lines), `drag_lift_drag` is V3Kite's `lift_drag` drag.
"""
function fly_points(points; t_trim, t_avg)
    s = init(PARK.v_wind, PARK.tether_length; body_start_damping = PARK.body_damping,
             body_sim_damping = 0.8 .* PARK.body_damping,
             damping_per_stiffness = PARK.damping_per_stiffness,
             depower_setpoint = PARK.depower, dt = PARK.dt,
             sim_time = length(points) * (t_trim + t_avg) + 1.0,
             system_yaml = PARK.project, use_turbulence = 0.0, remake_model = false)
    l0 = s.sys_state.l_tether[1]
    wpc = WinchPosController(WCSettings(; dt = s.dt); dt = s.dt)
    pid = create_heading_pid(; K = HEADING.P, Ti = HEADING.I, Td = HEADING.D, dt = s.dt,
                             umin = -HEADING.max_steering, umax = HEADING.max_steering)
    v_wind = PARK.v_wind
    v_app_f = s.sys_state.v_app
    n_trim, n_avg = round(Int, t_trim / s.dt), round(Int, t_avg / s.dt)
    results = OrderedDict{String, Any}[]
    for (tape, v_target) in points
        rel_depower = depower_length_to_percentage(tape) / 100
        acc = Dict(k => Float64[] for k in (:v_app, :lift, :drag, :drag_lift_drag, :drag_vsm,
                                            :drag_bridle, :drag_bridle_in_wing, :drag_tether,
                                            :drag_kcu, :aoa, :elevation, :force, :v_wind))
        for i in 1:(n_trim + n_avg)
            v_app = s.sys_state.v_app
            v_app_f += s.dt / WIND_TRIM.tau * (v_app - v_app_f)
            if i <= n_trim
                dv = WIND_TRIM.gain * v_wind * (v_target / v_app_f - 1)
                v_wind += s.dt * clamp(dv, -WIND_TRIM.max_rate, WIND_TRIM.max_rate)
            end
            set_K!(pid, HEADING.P * HEADING.v_app_ref / max(v_app, HEADING.v_app_min), 0.0,
                   s.sys_state.heading)
            rel_steering = pid(0.0, s.sys_state.heading, 0.0)
            step!(s; rel_depower, rel_steering, v_wind_gnd = v_wind,
                  set_torque = winch_torque!(wpc, s, l0))
            i <= n_trim && continue
            lift, drag = lift_drag(s)
            parts = drag_parts(s)
            push!(acc[:v_app], s.sys_state.v_app)
            push!(acc[:lift], lift)
            push!(acc[:drag], parts.vsm + parts.bridle)
            push!(acc[:drag_lift_drag], drag)
            for (k, v) in pairs(parts)
                push!(acc[Symbol(:drag_, k)], v)
            end
            push!(acc[:aoa], rad2deg(span_mean_aoa(s.sys)))
            push!(acc[:elevation], rad2deg(s.sys_state.elevation))
            push!(acc[:force], winch_force(s))
            push!(acc[:v_wind], v_wind)
        end
        res = OrderedDict{String, Any}("tape" => tape, "v_app_target" => v_target,
            "rel_depower" => rel_depower)
        for (k, v) in acc
            res[string(k)] = mean(v)
        end
        res["lift_std"] = std(acc[:lift])
        res["v_app_std"] = std(acc[:v_app])
        @printf("V3Kite  tape %.2f m  v_a %5.2f m/s (target %5.2f)  L %6.0f N  D %5.0f N (lift_drag %5.0f N)  AoA %5.2f°  wind %5.2f m/s\n",
                tape, res["v_app"], v_target, res["lift"], res["drag"], res["drag_lift_drag"],
                res["aoa"], res["v_wind"])
        push!(results, res)
    end
    return results
end

"""
    aero_forces(url, tape, v_app; rom = nothing) -> Dict

One `POST /aero_forces` at the ROM's trim. A request the server rejects returns
`Dict("error" => detail)` instead of throwing.
"""
function aero_forces(url, tape, v_app; rom = nothing)
    req = Dict{String, Any}("input_depower" => tape, "speed_apparent_wind" => v_app)
    isnothing(rom) || (req["rom_config_path"] = rom)
    try
        resp = HTTP.post(url * "/aero_forces", ["Content-Type" => "application/json"];
                         body = JSON3.write(req))
        return JSON3.read(resp.body, Dict{String, Any})
    catch err
        err isa HTTP.StatusError || rethrow()
        return Dict{String, Any}("error" => String(copy(err.response.body)))
    end
end

"""
    query_awetrim(v3kite, url) -> Vector{OrderedDict}

For each V3Kite block: the ROM forces at its tape and mean `v_a`, per ROM of `ROMS`,
at the ROM's trim (`trim`).
"""
function query_awetrim(v3kite, url)
    map(v3kite) do res
        OrderedDict{String, Any}(label => OrderedDict{String, Any}(
            "trim" => aero_forces(url, res["tape"], res["v_app"]; rom))
            for (label, rom) in ROMS)
    end
end

check_aero_forces_endpoint(url)
if fly
    v3kite = fly_points(points; t_trim, t_avg)
    flown = Dict("date" => string(Dates.now()), "git" => package_git_state(),
                 "park" => PARK, "t_trim" => t_trim, "t_avg" => t_avg)
else
    isfile(RESULT_FILE) || error("No $RESULT_FILE to re-query; run with fly = true first.")
    saved = JSON3.read(read(RESULT_FILE, String), Dict{String, Any})
    v3kite, flown = saved["v3kite"], saved["flown"]
end
awetrim = query_awetrim(v3kite, url)

for (res, aw) in zip(v3kite, awetrim), (label, _) in ROMS
    trim = aw[label]["trim"]
    fmt(r) = haskey(r, "error") ? "      error" : @sprintf("%5.0f / %4.0f", r["lift"], r["drag"])
    @printf("AWETrim %-16s tape %.2f m  v_a %5.2f m/s  trim %s N (%5.2f°)\n",
            label, res["tape"], res["v_app"], fmt(trim), get(trim, "angle_of_attack", NaN))
end

mkpath(dirname(RESULT_FILE))
open(RESULT_FILE, "w") do io
    JSON3.pretty(io, OrderedDict("flown" => flown, "queried" => string(Dates.now()),
                                 "url" => url, "roms" => first.(ROMS),
                                 "v3kite" => v3kite, "awetrim" => awetrim))
end
@info "validate_lift_drag: wrote $RESULT_FILE; show it with include(\"examples/show_validation.jl\")."
nothing
