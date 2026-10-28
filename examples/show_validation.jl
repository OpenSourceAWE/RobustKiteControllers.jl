# Copyright (c) 2026 Uwe Fechner
# SPDX-License-Identifier: MPL-2.0

"""
Show the lift/drag comparison saved by `validate_lift_drag.jl`
(`output/lift_drag_validation.json`) as an HTML page in the default browser: one table
per AWETrim ROM, each row a V3Kite block (power tape, measured apparent wind speed),
with V3Kite's lift and drag next to the ROM's at its own trim, and the deviation of the
ROM from V3Kite in percent. Rows outside the ROM's validity band are marked.

    include("examples/show_validation.jl")
"""

using Pkg
if Base.active_project() != joinpath(@__DIR__, "Project.toml")
    Pkg.activate(joinpath(@__DIR__))
end

using JSON3
using Printf

"The results `validate_lift_drag.jl` writes"
const RESULT_FILE = normpath(joinpath(@__DIR__, "..", "output", "lift_drag_validation.json"))
# A snap-confined browser cannot read dotdirs under $HOME, so the page goes to $HOME itself.
"The page written and opened"
const HTML_FILE = joinpath(homedir(), "lift_drag_validation.html")

const CSS = """
:root { --bg:#fff; --fg:#1d1d1f; --muted:#6e6e73; --line:#d9d9de; --head:#f3f3f6;
        --good:#1a7f37; --warn:#b26a00; --bad:#c62828; }
@media (prefers-color-scheme: dark) {
  :root { --bg:#16161a; --fg:#ececf1; --muted:#9a9aa3; --line:#34343c; --head:#202027;
          --good:#4cc26a; --warn:#e0a03a; --bad:#ef6b6b; } }
body { background:var(--bg); color:var(--fg); margin:0 auto; max-width:1100px; padding:24px 16px;
       font:15px/1.45 system-ui, -apple-system, "Segoe UI", sans-serif; }
h1 { font-size:22px; margin:0 0 4px; } h2 { font-size:17px; margin:28px 0 8px; }
.meta { color:var(--muted); font-size:13px; }
.wrap { overflow-x:auto; }
table { border-collapse:collapse; width:100%; font-variant-numeric:tabular-nums; }
th, td { padding:6px 10px; border-bottom:1px solid var(--line); text-align:right; white-space:nowrap; }
thead th { background:var(--head); font-weight:600; }
thead tr:first-child th { text-align:center; }
th.grp { border-left:2px solid var(--line); } td.grp { border-left:2px solid var(--line); }
.good { color:var(--good); } .warn { color:var(--warn); } .bad { color:var(--bad); font-weight:600; }
.out { color:var(--bad); }
ul { color:var(--muted); font-size:13px; padding-left:18px; }
"""

esc(x) = replace(string(x), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;")
num(x, d = 0) = x isa Number ? @sprintf("%.*f", d, x) : "–"
ld(r) = haskey(r, "error") ? "<span class=\"bad\" title=\"$(esc(r["error"]))\">error</span>" :
        "$(num(r["lift"])) / $(num(r["drag"]))"

"Deviation of `x` from the V3Kite value `ref` in percent, coloured by size"
function deviation(x, ref)
    (x isa Number && ref isa Number && ref != 0) || return "–"
    d = 100 * (x - ref) / ref
    cls = abs(d) < 5 ? "good" : abs(d) < 15 ? "warn" : "bad"
    return "<span class=\"$cls\">$(@sprintf("%+.0f", d)) %</span>"
end

"""
    rom_table(io, label, v3kite, awetrim)

One table: the V3Kite blocks against the ROM `label`.
"""
function rom_table(io, label, v3kite, awetrim)
    first_ok = something(findfirst(aw -> !haskey(aw[label]["trim"], "error"), awetrim), 0)
    rom_file = first_ok == 0 ? "" : basename(string(awetrim[first_ok][label]["trim"]["rom_config"]))
    println(io, "<h2>ROM: $(esc(label)) <span class=\"meta\">$(esc(rom_file))</span></h2>")
    println(io, """<div class="wrap"><table><thead>
      <tr><th colspan="2"></th><th colspan="2" class="grp">V3Kite</th>
          <th colspan="4" class="grp">AWETrim, ROM trim</th><th class="grp"></th></tr>
      <tr><th>tape [m]</th><th>v<sub>a</sub> [m/s]</th>
          <th class="grp">AoA [°]</th><th>L / D [N]</th>
          <th class="grp">AoA [°]</th><th>L / D [N]</th><th>ΔL</th><th>ΔD</th>
          <th class="grp">in ROM band</th></tr>
      </thead><tbody>""")
    for (res, aw) in zip(v3kite, awetrim)
        trim = aw[label]["trim"]
        valid = get(trim, "within_validity", nothing)
        band = get(trim, "input_depower_validity", nothing)
        valid_cell = valid === nothing ? "–" : valid ? "yes" :
            "<span class=\"out\" title=\"tape band $(esc(band))\">no</span>"
        println(io, "<tr><td>$(num(res["tape"], 2))</td><td>$(num(res["v_app"], 1))</td>",
                "<td class=\"grp\">$(num(res["aoa"], 2))</td><td>$(num(res["lift"])) / $(num(res["drag"]))</td>",
                "<td class=\"grp\">$(num(get(trim, "angle_of_attack", nothing), 2))</td><td>$(ld(trim))</td>",
                "<td>$(deviation(get(trim, "lift", nothing), res["lift"]))</td>",
                "<td>$(deviation(get(trim, "drag", nothing), res["drag"]))</td>",
                "<td class=\"grp\">$valid_cell</td></tr>")
    end
    println(io, "</tbody></table></div>")
end

"""
    drag_table(io, v3kite)

The V3Kite drag split per block, along the wing's apparent wind; skipped for results
saved before the split was recorded.
"""
function drag_table(io, v3kite)
    all(res -> haskey(res, "drag_vsm"), v3kite) || return
    println(io, """<h2>V3Kite drag split</h2><div class="wrap"><table><thead><tr>
      <th>tape [m]</th><th>v<sub>a</sub> [m/s]</th><th class="grp">VSM wing [N]</th>
      <th>bridle lines [N]</th><th>compared D [N]</th><th class="grp">bridle in <code>lift_drag</code> [N]</th>
      <th><code>lift_drag</code> D [N]</th><th class="grp">tether [N]</th><th>KCU [N]</th>
      </tr></thead><tbody>""")
    for res in v3kite
        println(io, "<tr><td>$(num(res["tape"], 2))</td><td>$(num(res["v_app"], 1))</td>",
                "<td class=\"grp\">$(num(res["drag_vsm"]))</td><td>$(num(res["drag_bridle"]))</td>",
                "<td>$(num(res["drag"]))</td>",
                "<td class=\"grp\">$(num(res["drag_bridle_in_wing"]))</td><td>$(num(res["drag_lift_drag"]))</td>",
                "<td class=\"grp\">$(num(res["drag_tether"]))</td><td>$(num(res["drag_kcu"]))</td></tr>")
    end
    println(io, "</tbody></table></div>")
end

"""
    write_validation_html(data, file)

Write the page for the parsed `lift_drag_validation.json` `data` to `file`.
"""
function write_validation_html(data, file)
    flown = data["flown"]
    git = flown["git"]
    open(file, "w") do io
        println(io, """<!DOCTYPE html><html lang="en"><head><meta charset="UTF-8">
          <meta name="viewport" content="width=device-width, initial-scale=1">
          <title>Lift/Drag Validation</title><style>$CSS</style></head><body>
          <h1>Lift and drag: V3Kite vs. AWETrim ROM</h1>
          <p class="meta">Flown $(esc(flown["date"])) (RobustKiteControllers $(esc(git["hash"]))
          $(esc(git["status"]))), AWETrim queried $(esc(data["queried"])) at $(esc(data["url"])).
          Parked at $(num(flown["park"]["tether_length"])) m, $(esc(flown["park"]["project"])),
          averaged over $(num(flown["t_avg"])) s per point.</p>""")
        for label in data["roms"]
            rom_table(io, label, data["v3kite"], data["awetrim"])
        end
        drag_table(io, data["v3kite"])
        println(io, """<ul>
          <li>V3Kite L: <code>lift_drag</code>, the wing's force perpendicular to its apparent
              wind. V3Kite D: the VSM wing drag plus the drag of all bridle lines, AWETrim's
              basis; <code>lift_drag</code>'s own D holds only the bridle share lumped onto the
              wing nodes.</li>
          <li>AWETrim L / D: the ROM wing (bridle lines included in D, KCU not), at the
              measured mean v<sub>a</sub> of the V3Kite block.</li>
          <li>ROM trim: the ROM's own gravity-free angle of attack. ΔL, ΔD: AWETrim relative
              to V3Kite.</li>
        </ul></body></html>""")
    end
    return file
end

"""
    open_in_browser(file)

Open `file` in the default browser on Linux, Windows or macOS; print its path if that fails.
"""
function open_in_browser(file)
    cmd = Sys.iswindows() ? `cmd /c start "" $file` : Sys.isapple() ? `open $file` : `xdg-open $file`
    try
        run(pipeline(cmd; stdout = devnull, stderr = devnull); wait = false)
    catch
        println("Could not open a browser; open $file manually.")
    end
    nothing
end

isfile(RESULT_FILE) || error("No $RESULT_FILE; run examples/validate_lift_drag.jl first.")
data = JSON3.read(read(RESULT_FILE, String), Dict{String, Any})
@info "show_validation: wrote $(write_validation_html(data, HTML_FILE))"
open_in_browser(HTML_FILE)
nothing
