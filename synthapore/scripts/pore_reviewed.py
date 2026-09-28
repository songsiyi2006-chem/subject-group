"""Analytic 2D probe-center geometry and synthetic BET inverse benchmarks.

No atomistic POP, material mass, GCMC, experimental adsorption or pore-size
inversion is present. All geometric lengths are angstrom; areas are angstrom^2.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / "synthapore/source/synthapore_engine.py"
OUTPUT = ROOT / "synthapore/results/pore"
NA = 6.02214076e23
N2_CROSS_SECTION_M2 = 0.162e-18
STP_MOLAR_VOLUME_CM3_MOL = 22414.0  # stated 273.15 K, 1 atm convention
AREA_PER_MONOLAYER_CM3 = NA * N2_CROSS_SECTION_M2 / STP_MOLAR_VOLUME_CM3_MOL
WINDOWS = [(0.01, 0.30), (0.05, 0.30), (0.10, 0.35), (0.01, 0.45),
           (0.05, 0.50), (0.15, 0.70), (0.35, 0.95), (0.01, 0.95)]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("No rows supplied")
    keys = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=keys, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def positive(value: float, name: str) -> float:
    if not np.isfinite(value) or value <= 0:
        raise ValueError(name + " must be positive and finite")
    return float(value)


def nonnegative(value: float, name: str) -> float:
    if not np.isfinite(value) or value < 0:
        raise ValueError(name + " must be nonnegative and finite")
    return float(value)


def grid_points(cell_A: float, resolution: int, rule: str = "midpoint") -> np.ndarray:
    positive(cell_A, "cell_A")
    if not isinstance(resolution, (int, np.integer)) or isinstance(resolution, bool) or resolution < 2:
        raise ValueError("resolution must be an integer >= 2")
    if rule == "midpoint":
        axis = (np.arange(resolution) + 0.5) * cell_A / resolution - cell_A/2
    elif rule == "source_endpoints":
        axis = np.linspace(-cell_A/2, cell_A/2, resolution)
    else:
        raise ValueError("Unknown grid rule")
    x, y = np.meshgrid(axis, axis, indexing="xy")
    return np.column_stack((x.ravel(), y.ravel()))


def _points(points) -> np.ndarray:
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError("points must be a finite (n,2) array")
    return points


def channel_clearance(points, shape: str, inradius_A: float = 9.8) -> np.ndarray:
    """Nearest wall distance inside the channel; hex outside values are plane slacks."""
    pts = _points(points)
    radius = positive(inradius_A, "inradius_A")
    if shape == "circle":
        return radius - np.linalg.norm(pts, axis=1)
    if shape == "hexagon":
        angles = np.arange(6) * np.pi/3
        normals = np.column_stack((np.cos(angles), np.sin(angles)))
        return radius - np.max(pts @ normals.T, axis=1)
    raise ValueError("shape must be circle or hexagon")


def channel_area(shape: str, inradius_A: float, probe_A: float = 0.0) -> float:
    radius = positive(inradius_A, "inradius_A")
    probe = nonnegative(probe_A, "probe_A")
    remaining = max(radius-probe, 0.0)
    if shape == "circle":
        return float(np.pi*remaining**2)
    if shape == "hexagon":
        return float(2*np.sqrt(3)*remaining**2)
    raise ValueError("Unknown channel shape")


def minimum_image(delta, cell_A: float) -> np.ndarray:
    length = positive(cell_A, "cell_A")
    values = np.asarray(delta, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Displacements must be finite")
    return (values + length/2) % length - length/2


def periodic_disk_clearance(points, cell_A: float = 26.0,
                            center_A=(12.4, 2.5), radius_A: float = 3.0) -> np.ndarray:
    points = _points(points)
    center = np.asarray(center_A, dtype=float)
    if center.shape != (2,) or not np.isfinite(center).all():
        raise ValueError("center_A must be a finite pair")
    positive(radius_A, "radius_A")
    if radius_A >= positive(cell_A, "cell_A")/2:
        raise ValueError("Benchmark solid disk radius must be less than half the cell")
    return np.linalg.norm(minimum_image(points-center, cell_A), axis=1) - radius_A


def disk_square_intersection_area(radius_A: float, cell_A: float) -> float:
    """Exact disk area inside its centered minimum-image square, including overlap regime."""
    r = nonnegative(radius_A, "radius_A")
    half = positive(cell_A, "cell_A")/2
    if r <= half:
        return float(np.pi*r*r)
    if r >= np.sqrt(2)*half:
        return float(cell_A**2)
    segment = r*r*np.arccos(half/r) - half*np.sqrt(r*r-half*half)
    return float(np.pi*r*r - 4*segment)


def periodic_accessible_area(cell_A: float, radius_A: float, probe_A: float) -> float:
    if not 0 < radius_A < positive(cell_A, "cell_A")/2:
        raise ValueError("Require 0 < disk radius < cell/2")
    nonnegative(probe_A, "probe_A")
    return float(cell_A**2 - disk_square_intersection_area(radius_A+probe_A, cell_A))


def clearance_cdf(distance_A, shape: str, inradius_A=9.8, cell_A=26.0, disk_radius_A=3.0) -> np.ndarray:
    values = np.asarray(distance_A, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError("Distances must be finite")
    if shape in ("circle", "hexagon"):
        radius = positive(inradius_A, "inradius_A")
        return 1-(1-np.clip(values/radius, 0, 1))**2
    if shape == "periodic_disk":
        total = periodic_accessible_area(cell_A, disk_radius_A, 0)
        return np.array([0.0 if t <= 0 else 1-periodic_accessible_area(cell_A, disk_radius_A, float(t))/total for t in values])
    raise ValueError("Unknown shape")


def bet_volume(relative_pressure, monolayer_cm3_g=145.0, C=115.0) -> np.ndarray:
    p = np.asarray(relative_pressure, dtype=float)
    if not np.isfinite(p).all() or np.any(p <= 0) or np.any(p >= 1):
        raise ValueError("Relative pressures must satisfy 0 < p < 1")
    positive(monolayer_cm3_g, "monolayer_cm3_g")
    positive(C, "C")
    return monolayer_cm3_g*C*p / ((1-p)*(1+(C-1)*p))


def source_synthetic_volume(relative_pressure) -> np.ndarray:
    p = np.asarray(relative_pressure, dtype=float)
    low = bet_volume(p)
    return np.where(p < 0.35, low, 145.0*3.8*(1-np.exp(-4.2*p)))


def fit_bet(relative_pressure, volume_cm3_g, lower=0.01, upper=0.30) -> dict:
    """OLS BET linearization with selected consistency diagnostics, not full certification."""
    p, v = np.asarray(relative_pressure, dtype=float), np.asarray(volume_cm3_g, dtype=float)
    if p.ndim != 1 or p.shape != v.shape or not np.isfinite(p).all() or not np.isfinite(v).all():
        raise ValueError("Pressure and volume must be finite matching vectors")
    if np.any(p <= 0) or np.any(p >= 1) or np.any(v <= 0):
        raise ValueError("Require 0 < p < 1 and positive adsorption volume")
    if not np.isfinite([lower, upper]).all() or not 0 < lower < upper < 1:
        raise ValueError("Require 0 < lower < upper < 1")
    selected = (p >= lower) & (p <= upper)
    x, observed = p[selected], v[selected]
    if len(x) < 3 or len(np.unique(x)) != len(x):
        raise ValueError("At least three distinct pressure points are required per fit")
    order = np.argsort(x)
    x, observed = x[order], observed[order]
    y = x/(observed*(1-x))
    design = np.column_stack((np.ones_like(x), x))
    intercept, slope = np.linalg.lstsq(design, y, rcond=None)[0]
    yhat = design @ np.array([intercept, slope])
    denominator = slope+intercept
    vm = 1/denominator if denominator != 0 else None
    c = 1+slope/intercept if intercept != 0 else None
    valid = vm is not None and c is not None and vm > 0 and c > 0 and intercept > 0
    total_ss = float(np.sum((y-y.mean())**2))
    mono_p = float(1/(1+np.sqrt(c))) if valid else None
    monotonic_rouquerol = bool(np.all(np.diff(observed*(1-x)) >= -1e-12))
    result = {"lower_requested": lower, "upper_requested": upper, "points": len(x),
        "lowest_used_p": float(x[0]), "highest_used_p": float(x[-1]),
        "intercept": float(intercept), "slope": float(slope),
        "fitted_Vm_cm3_g": float(vm) if vm is not None else None,
        "fitted_C": float(c) if c is not None else None,
        "positive_parameters": bool(valid),
        "transformed_R2": float(1-np.sum((y-yhat)**2)/total_ss) if total_ss > 0 else None,
        "transformed_RMSE": float(np.sqrt(np.mean((y-yhat)**2))),
        "v_times_one_minus_p_nondecreasing": monotonic_rouquerol,
        "BET_monolayer_relative_pressure": mono_p,
        "monolayer_pressure_inside_window": bool(x[0] <= mono_p <= x[-1]) if mono_p is not None else False,
        "apparent_area_m2_g": float(vm*AREA_PER_MONOLAYER_CM3) if valid else None,
        "source_conversion_area_m2_g": float(vm*4.353) if valid else None,
        "volume_RMSE_cm3_g": float(np.sqrt(np.mean((bet_volume(x, vm, c)-observed)**2))) if valid else None,
        "truth_Vm_relative_error": float((vm-145)/145) if vm is not None else None,
        "truth_C_relative_error": float((c-115)/115) if c is not None else None,
        "evidence": "Synthetic inverse fit; diagnostics are not a complete BET validity assessment"}
    return result


def source_class_audit() -> dict:
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    node = next(n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == "ReticularPOPEngine")
    namespace = {"np": np, "Dict": dict}
    # Execute only this immutable, dependency-free numerical class, never the source master workflow.
    exec(compile(ast.Module(body=[node], type_ignores=[]), "<source-pore-class>", "exec"), namespace)
    generator = namespace["ReticularPOPEngine"].generate_hexagonal_pop_lattice
    calls = []
    outputs = []
    for cell in (24.5, 26.0):
        for probe in (0.0, 1.82, 4.0, 10.0):
            output = generator(unit_cell_size_angstrom=cell, pore_probe_radius=probe)
            outputs.append(output)
            calls.append({"cell_A": cell, "probe_A": probe,
                "source_reported_porosity_pct": output["calculated_void_porosity_pct"],
                "source_reported_BET_area_m2_g": output["bet_specific_surface_area_m2_g"]})
    loads = [n.id for n in ast.walk(node) if isinstance(n, ast.Name) and isinstance(n.ctx, ast.Load)]
    boundary = 0.35
    left = float(bet_volume(boundary))
    right = float(145*3.8*(1-np.exp(-4.2*boundary)))
    return {"class_only_execution_count": len(calls), "calls": calls,
        "probe_argument_read_in_body": "pore_probe_radius" in loads,
        "source_mask_shape": "circle: X^2+Y^2 < 9.8^2 on an endpoint-including square grid",
        "topology_claim_supported_by_graph": False, "atomic_framework_coordinates_supplied": False,
        "atomic_masses_or_material_mass_used": False, "PSD_computed": False,
        "source_area_derivation": "145 cm3(STP)/g prescribed monolayer times prescribed 4.353 conversion; independent of cell geometry",
        "source_type_IV_label_validated": False,
        "source_unit_cell_26_output": outputs[4],
        "branch_boundary_p": boundary, "BET_left_limit_cm3_g": left,
        "empirical_right_value_cm3_g": right, "jump_cm3_g": right-left,
        "relative_jump_vs_left": (right-left)/left,
        "interpretation": "Discontinuity imposed by switching formula, not evidence of capillary condensation or hysteresis"}


def run_benchmark(output: Path = OUTPUT, pilot: bool = False) -> dict:
    started = time.perf_counter()
    output.mkdir(parents=True, exist_ok=True)
    cell, channel_radius, disk_radius = 26.0, 9.8, 3.0
    resolutions = [40, 80] if pilot else [40, 80, 160, 320, 640]
    probes = [0.0, 1.2, 1.82, 3.0, 6.0, 9.8, 10.0]
    periodic_probes = [0.0, 1.82, 3.0, 6.0]
    geometry, endpoint_rows, cdf_rows, histogram_rows = [], [], [], []
    cached = {}
    field_points = 0
    for n in resolutions:
        points = grid_points(cell, n)
        for shape in ("circle", "hexagon", "periodic_disk"):
            clearance = periodic_disk_clearance(points, cell, radius_A=disk_radius) if shape == "periodic_disk" else channel_clearance(points, shape, channel_radius)
            field_points += len(points)
            if n == resolutions[-1]:
                cached[shape] = clearance
            for probe in periodic_probes if shape == "periodic_disk" else probes:
                count = int(np.count_nonzero(clearance > probe))
                area = cell**2*count/len(points)
                exact = periodic_accessible_area(cell, disk_radius, probe) if shape == "periodic_disk" else channel_area(shape, channel_radius, probe)
                geometry.append({"geometry": shape, "resolution": n, "cell_A": cell,
                    "channel_inradius_or_disk_radius_A": disk_radius if shape == "periodic_disk" else channel_radius,
                    "probe_A": probe, "grid_points": len(points), "accessible_points": count,
                    "midpoint_area_A2": area, "analytic_area_A2": exact,
                    "signed_area_error_A2": area-exact, "abs_area_error_A2": abs(area-exact),
                    "grid_accessible_fraction": area/cell**2, "analytic_accessible_fraction": exact/cell**2,
                    "abs_fraction_error": abs(area-exact)/cell**2})
        for source_cell in (24.5, 26.0):
            source_points = grid_points(source_cell, n, "source_endpoints")
            midpoint = grid_points(source_cell, n)
            exact_fraction = np.pi*channel_radius**2/source_cell**2
            for rule, pts in (("source_endpoints", source_points), ("midpoint", midpoint)):
                fraction = float(np.mean(np.sum(pts**2, axis=1) < channel_radius**2))
                endpoint_rows.append({"cell_A": source_cell, "resolution": n, "rule": rule,
                    "estimated_circle_area_fraction": fraction, "analytic_circle_area_fraction": exact_fraction,
                    "signed_error_percentage_points": 100*(fraction-exact_fraction)})
    for shape, clearance in cached.items():
        free = clearance[clearance > 0]
        maximum = np.sqrt(2)*cell/2-disk_radius if shape == "periodic_disk" else channel_radius
        edges = np.linspace(0, maximum, 41)
        counts, _ = np.histogram(free, bins=edges)
        exactcdf = clearance_cdf(edges, shape, channel_radius, cell, disk_radius)
        for i, value in enumerate(edges):
            estimated = float(np.mean(free <= value))
            cdf_rows.append({"geometry": shape, "resolution": resolutions[-1], "clearance_A": float(value),
                "empirical_area_weighted_CDF": estimated, "analytic_CDF": float(exactcdf[i]),
                "absolute_CDF_error": abs(estimated-exactcdf[i])})
        for i, count in enumerate(counts):
            histogram_rows.append({"geometry": shape, "lower_clearance_A": float(edges[i]),
                "upper_clearance_A": float(edges[i+1]), "count": int(count), "grid_void_points": len(free),
                "empirical_probability": float(count/len(free)), "analytic_probability": float(exactcdf[i+1]-exactcdf[i]),
                "interpretation": "2D point-to-wall clearance distribution, not experimental or atomistic PSD"})
    source_audit = source_class_audit()
    source_p = np.array(source_audit["source_unit_cell_26_output"]["relative_pressures_p_p0"])
    source_v = np.array(source_audit["source_unit_cell_26_output"]["n2_adsorption_isotherm_cm3_g"])
    data_rows = [{"point": i, "relative_pressure": float(p), "source_rounded_volume_cm3_g": float(v),
                  "source_unrounded_volume_cm3_g": float(source_synthetic_volume(p)),
                  "BET_only_control_cm3_g": float(bet_volume(p)), "branch": "BET" if p < .35 else "empirical_saturation"}
                 for i, (p, v) in enumerate(zip(source_p, source_v))]
    fits = []
    for label, observed in [("source_rounded", source_v), ("source_unrounded", source_synthetic_volume(source_p)),
                            ("BET_exact_control", bet_volume(source_p))]:
        for lower, upper in WINDOWS:
            fits.append(dict(data=label, **fit_bet(source_p, observed, lower, upper)))
    noise_rows, noise_fits = [], []
    seeds = list(range(8 if pilot else 64))
    low_p = source_p[source_p <= .30]
    for seed in seeds:
        rng = np.random.default_rng(202609280+seed)
        noise_v = bet_volume(low_p)*np.exp(rng.normal(-0.5*.01**2, .01, len(low_p)))
        noise_rows.extend({"seed_index": seed, "actual_seed": 202609280+seed, "point": i,
            "relative_pressure": float(p), "noisy_volume_cm3_g": float(v),
            "noise_model": "Independent multiplicative lognormal, log SD 0.01, mean factor 1"}
            for i, (p, v) in enumerate(zip(low_p, noise_v)))
        for lower, upper in WINDOWS[:2]:
            noise_fits.append(dict(seed_index=seed, **fit_bet(low_p, noise_v, lower, upper)))
    noise_summary = []
    for lower, upper in WINDOWS[:2]:
        group = [f for f in noise_fits if f["lower_requested"] == lower]
        valid = [f for f in group if f["positive_parameters"]]
        for param in ("fitted_Vm_cm3_g", "fitted_C", "apparent_area_m2_g"):
            vals = np.array([f[param] for f in valid])
            noise_summary.append({"lower": lower, "upper": upper, "parameter": param,
                "total_fits": len(group), "positive_fits": len(valid),
                "mean": float(np.mean(vals)), "sample_sd": float(np.std(vals, ddof=1)),
                "p025": float(np.quantile(vals, .025)), "median": float(np.median(vals)), "p975": float(np.quantile(vals, .975)),
                "interpretation": "Synthetic-noise repetition quantiles conditional on model, not experimental confidence intervals"})
    tables = {"geometry_grid_convergence": geometry, "endpoint_bias": endpoint_rows,
        "clearance_CDF": cdf_rows, "clearance_histogram": histogram_rows, "source_isotherm_points": data_rows,
        "BET_window_fits": fits, "noise_observations": noise_rows, "noise_BET_fits": noise_fits, "noise_summary": noise_summary}
    for name, rows in tables.items():
        write_csv(output/(name+".csv"), rows)
    write_json(output/"source_pore_audit.json", source_audit)
    rng = np.random.default_rng(8211)
    points = rng.uniform(-cell/2, cell/2, (64, 2))
    center = np.array([12.4, 2.5])
    shifts = np.array([(i*cell, j*cell) for i in (-1,0,1) for j in (-1,0,1)])
    brute = np.min(np.linalg.norm(points[:, None, :]-center[None, None, :]-shifts[None, :, :], axis=2), axis=1)-disk_radius
    mi_error = float(np.max(np.abs(brute-periodic_disk_clearance(points, cell))))
    checks = {"minimum_image_vs_nine_images_max_A": mi_error,
              "source_probe_has_no_effect": not source_audit["probe_argument_read_in_body"] and len({r["source_reported_porosity_pct"] for r in source_audit["calls"] if r["cell_A"] == 26}) == 1,
              "exact_control_max_Vm_abs_error_cm3_g": max(abs(f["fitted_Vm_cm3_g"]-145) for f in fits if f["data"] == "BET_exact_control"),
              "exact_control_max_C_abs_error": max(abs(f["fitted_C"]-115) for f in fits if f["data"] == "BET_exact_control"),
              "finest_max_area_fraction_error": max(r["abs_fraction_error"] for r in geometry if r["resolution"] == resolutions[-1]),
              "clearance_CDF_max_error": max(r["absolute_CDF_error"] for r in cdf_rows)}
    summary = {"schema_version": 1, "pilot": pilot,
        "evidence": "Executed analytic 2D geometry and synthetic inverse-model benchmarks; no physical porous material validated",
        "geometry": {"cell_A": cell, "channel_inradius_A": channel_radius, "disk_radius_A": disk_radius,
            "periodic_disk_center_A": center.tolist(), "grid_rule": "Equal-area cell centers, strict clearance > probe",
            "resolutions": resolutions, "channel_probes_A": probes, "periodic_probes_A": periodic_probes,
            "exact_channel_area": "circle pi*(R-rp)^2; regular hexagon 2*sqrt(3)*(a-rp)^2; zero once rp >= inradius",
            "hexagon_convention": "apothem/inradius 9.8 A, circumradius 2*a/sqrt(3); six half-plane normals",
            "periodic_geometry": "One artificial solid disk in a periodic square; excluded distance expanded by probe radius",
            "clearance_definition": "Area-weighted distance of free 2D points to nearest wall; not a histogram of pore diameters"},
        "adsorption": {"source_prescribed_Vm_cm3_g": 145.0, "source_prescribed_C": 115.0,
            "branch_switch": .35, "STP_molar_volume_cm3_mol": STP_MOLAR_VOLUME_CM3_MOL,
            "STP_convention": "273.15 K and 1 atm, rounded ideal-gas molar volume 22414 cm3/mol",
            "nitrogen_cross_section_m2": N2_CROSS_SECTION_M2, "Avogadro_mol_1": NA,
            "derived_area_factor_m2_per_cm3_STP": AREA_PER_MONOLAYER_CM3,
            "source_area_factor": 4.353, "regression": "Unweighted OLS on p/[V*(1-p)] versus p, intercept and slope retained",
            "validity_scope": "Positive parameters, V*(1-p) monotonicity and monolayer pressure placement are diagnostics only; no complete standard-compliance claim",
            "noise": "64 seeds (pilot 8), 1% lognormal multiplicative error, positive mean-one factor; no inferred measurement noise"},
        "counts": {"clearance_fields": 3*len(resolutions), "clearance_point_evaluations": field_points,
            "probe_masks": len(geometry), "endpoint_or_midpoint_bias_cases": len(endpoint_rows),
            "clearance_CDF_rows": len(cdf_rows), "clearance_histogram_rows": len(histogram_rows),
            "source_class_calls": source_audit["class_only_execution_count"], "source_isotherm_points": len(source_p),
            "deterministic_BET_fits": len(fits), "noise_seeds": len(seeds), "noise_observations": len(noise_rows),
            "noise_BET_fits": len(noise_fits), "total_BET_fits": len(fits)+len(noise_fits),
            "GCMC_runs": 0, "atomistic_frameworks": 0, "experimental_measurements": 0},
        "numerical_checks": checks,
        "limitations": ["2D analytic geometry is not a chemical hcb/sql framework; no connectivity/chemistry assignment is made.",
            "Area of allowed probe centers differs from probe-occupiable area, surface area and 3D volume.",
            "Clearance histograms are not experimental or atomistic pore-size distributions.",
            "BET areas derive from prescribed synthetic mass-normalized uptake, not from atomistic mass or geometry.",
            "Good R2 or positive BET C does not validate a fitting window or a physical material.",
            "Noise quantiles are conditional synthetic-repeat variability; no adsorption/desorption hysteresis or experiment is supplied."],
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "elapsed_seconds": time.perf_counter()-started},
        "sources": [{"url": "https://www.zeoplusplus.org/examples.html", "role": "Definitions distinguish probe-accessible centers, occupied volume, surface, and PSD; software not executed"},
            {"url": "https://nvlpubs.nist.gov/nistpubs/Legacy/SP/nistspecialpublication960-17.pdf", "role": "BET linearization, limitations and conventional N2 cross-section"},
            {"url": "https://mail.goldbook.iupac.org/terms/view/13997", "role": "BET area is isotherm-derived with conventional molecular cross-section; search-index evidence, main host returned 403"}],
        "license_note": "No external structure, dataset or model downloaded; original source retained elsewhere. Analytic benchmark data generated locally.",
        "source_sha256": {"synthapore/scripts/pore_reviewed.py": digest(Path(__file__)), "synthapore/source/synthapore_engine.py": digest(SOURCE)},
        "output_sha256": {p.name: digest(p) for p in sorted(output.iterdir()) if p.is_file() and p.name != "summary.json"}}
    write_json(output/"summary.json", summary)
    print(json.dumps({"counts": summary["counts"], "checks": checks}, indent=2), flush=True)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    run_benchmark(args.output_dir or (OUTPUT/"pilot" if args.pilot else OUTPUT), args.pilot)
