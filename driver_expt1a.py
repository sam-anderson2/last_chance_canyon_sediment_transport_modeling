"""
driver_expt1a.py ― Last-Chance Creek single-storm driver
Author: Sam Anderson
========================================================
This script couples Landlab **OverlandFlow** (hydraulics) with 
**RiverBedDynamics** (mixed-grain sediment transport).

I use it to drive one storm over a DEM with spatially variable surface GSDs.  
The run streams stats to disk (low RAM), is NaN-robust, and produces a single,
coherent set of artifacts that tell me — in plain terms — **what moved, where,
when, and by how much**.

----------------------------------------------------------------------------
QUICK START / CLI EXAMPLES
----------------------------------------------------------------------------

# LC-1 demo, one storm, with quick-look PNGs
python driver_expt1a.py --rainfall working_climate/15min_1000yr.txt --plots

# Batch run: every rainfall file on all CPU cores
python run_many_storms.py --sim_max_t 25000 --plots

# Switch to LC-3 instead of LC-1
python driver_expt1a.py --z_dem filled_lc3_dem.asc \
                        --gsd_loc_raster lc3_gsd_locations.asc \
                        --gsd_file GSDs/LC3_grain_size_dist_with_max.txt

----------------------------------------------------------------------------
EXPECTED INPUTS (disk layout and file formats)
----------------------------------------------------------------------------
project_root/
│
├─ driver_expt1a.py          (← this script)       ┐
├─ run_many_storms.py                                │  Python scripts
│                                                   ┘
├─ filled_lc1_dem.asc   filled_lc3_dem.asc          # DEMs (ESRI ASCII grid)
├─ lc1_gsd_locations.asc   lc3_gsd_locations.asc    # raster: node→GSD section (0/1/2)
│
├─ GSDs/LC1_grain_size_dist_with_max.txt
│   GSDs/LC3_grain_size_dist_with_max.txt           # 4-col text table:
│                                                   #   size_mm, upstream%, mid%, downstream%
│
└─ working_climate/*.txt   # rainfall hyetographs (2-col text):
                           #   time_s   intensity_m/s

RUN OUTPUTS (core + diagnostics)  — created in run_<tag>/ or ./output/
----------------------------------------------------------------------------
File                                       | Units   | What it shows (in plain English)
-------------------------------------------|---------|----------------------------------------------
elev_change.asc                            | m       | Net bed elevation change per node (erosion negative, deposition positive).
D10_change.asc                             | mm      | Change in the 10th-percentile surface grain size (fines shift). Positive = coarser.
D50_change.asc                             | mm      | Change in median surface grain size. Positive = coarsening; negative = fining.
D90_change.asc                             | mm      | Change in the 90th-percentile surface grain size (coarse tail behavior).
Dmax_change.asc                            | mm      | Change in the coarsest present size on the surface (before vs after).
gsd_section_thresholds_from_table.txt      | –       | For each section: max input size present in table and max eligible given FRAC_MIN.
gsd_init_diagnostics.txt                   | –       | For each section: largest size with any>0, >FRAC_MIN, and >PRES_STRICT at t=0.
sampled_links_validation.txt               | –       | Which hard-coded link IDs were kept/dropped for this grid (bounds check).
channel_links_section_check.txt            | –       | Any requested channel links that didn’t belong to their intended section (and were dropped).
channel_probe_summary.txt                  | mixed   | For each probed link: mean t, depth, τ, qb, and D50 over the recorded series.
channel_probe_<sec>_link<ID>_timeseries.txt| mixed   | Full time series at each probed link: t_s, depth_m, τ_Pa, qb_m2s, and D50_mm.
channel_probe_<sec>_link<ID>_before.txt    | –       | Per-link GSD snapshot at t=0 (fraction + CDF); verifies initial conditions.
channel_probe_<sec>_link<ID>_after.txt     | –       | Per-link GSD snapshot at t=end (fraction + CDF); direct before/after comparison.
channel_probe_<sec>_link<ID>_diff.txt      | –       | Per-link Δfraction (after − before) by size bin.
channel_probe_<sec>_link<ID>_diff_b_minus_a.txt | –  | Per-link Δfraction (before − after) by size bin.
section_<section>_gsd_before.txt           | –       | Section-mean surface GSD at t=0 (fraction + CDF).
section_<section>_gsd_after.txt            | –       | Section-mean surface GSD at t=end (fraction + CDF).
section_<section>_diff.txt                 | –       | Section-mean Δfraction (after − before) by size bin.
section_<section>_diff_b_minus_a.txt       | –       | Section-mean Δfraction (before − after) by size bin.
section_<section>_largest_moved_mm_timeseries.txt | mm | Largest grain size moved per step for that section.
section_<section>_largest_moved_mm_final.txt | mm    | Final (run-max) largest grain size moved for that section.
run_health.json                            | –       | Run metadata: storm file, t_end, max depth, any RBD steps, and link counts per section.

----------------------------------------------------------------------------
WHAT THE SCRIPT DOES (workflow)
----------------------------------------------------------------------------
• I set switches for bedload law, thresholds, rainfall file, and DEM/GSD paths.  
• I build a Landlab grid from a DEM + a raster of section classes (0=up,1=mid,2=down).  
• I initialize OverlandFlow with a rainfall hyetograph (time, intensity).  
• I couple hydraulics to RiverBedDynamics to compute shear stress and bedload.  
• I track the largest grain moved per link and aggregate section-level diagnostics.  
• I stream summaries and write post-event rasters and text tables for GIS/plots.  
• Optional PNGs provide quick reality checks (but slow the run).

----------------------------------------------------------------------------
PERFORMANCE NOTES
----------------------------------------------------------------------------
• On my Razer 15 laptop, LC-3 takes ~2 days with a “seeded bed” and ~4 days if unseeded.  
• LC-1 is faster. PNG plotting adds noticeable overhead, so I usually disable --plots in batches.
"""


# ── ONE-STOP USER SWITCHBOARD ────────────────────────────────────────────────
# I'm putting all the knobs and switches I routinely touch in one place so I can
# quickly re-run experiments without hunting through the code. If I need to run
# Parker instead of Wilcock & Crowe or seed an initial surface, I just toggle here.
BEDLOAD_EQUATION = "WilcockAndCrowe"      # My bedload equation choice: "WilcockAndCrowe" or "Parker"
SEED_SURFACE     = False                   # If True, I initialize every surface GSD bin to a uniform 1/N

# These thresholds gate whether grains are considered present and whether hydraulics are active.
FRAC_MIN         = 1e-6                    # I ignore any GSD bin with fraction below this tiny cutoff
DEPTH_MIN        = 0.000005                # Lowered so low-intensity storms can couple

# LC-1 vs LC-3 selector (I keep these three filenames in sync when switching study site)
DEM_NAME       = "filled_lc1_dem.asc"      # DEM for LC-1. If I switch to LC-3: "filled_lc3_dem.asc"
GSD_TABLE      = "GSDs/LC1_grain_size_dist_with_max.txt"  # GSD CDF table (mm + 3 columns for sections)
GSD_LOC_RASTER = "lc1_gsd_locations.asc"   # Raster with codes 0/1/2 mapping nodes to upstream/mid/downstream

# I use these section labels throughout to group outputs and write file names consistently.
SECTS = ("upstream", "intermediate", "downstream")

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# For reproducible channel diagnostics, I hard-code a few link IDs per section
# for LC3 and LC1. I validate them against the active grid so any out-of-range
# IDs are logged and dropped cleanly (no crashes).
SAMPLED_LINKS = {
    "LC3": {
        "upstream":     [14774, 22260, 30084],
        "intermediate": [40275, 41481, 43010],
        "downstream":   [46743, 60313, 67629],
    },
    "LC1": {
        "upstream":     [11471, 32127, 59956],
        "intermediate": [60120, 62572, 65185],
        "downstream":   [65348, 68282, 71051],
    },
}
# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

WARMUP_S          = 5.0
SUSTAIN_STEPS     = 1
USE_RH_TAU        = True
TRACK_STRATIGRAPHY= True

# ── IMPORTS ───────────────────────────────────────────────────────────────────
import matplotlib, math, random, warnings, numpy as np, json, re
matplotlib.use("Agg")
from pathlib import Path
from shutil  import rmtree
from os      import mkdir, getcwd
import matplotlib.pyplot as plt
from numpy.ma import masked_equal
from landlab.io import read_esri_ascii, write_esri_ascii
from landlab.components import OverlandFlow, RiverBedDynamics, FlowAccumulator, FlowDirectorSteepest
from rainfall_manager import rainfall_data, update_rainfall_intensity
from update_time      import update_time
warnings.filterwarnings("ignore", category=UserWarning)

# -----------------------------------------------------------------------------
# Small filesystem helper: I create a clean output directory per run tag.
# ⚙️ EDITED: now robust to race conditions / missing dirs.
def _clean(root: Path, tag: str) -> Path:
    """Create fresh output directory 'run_<tag>' (or 'output' if tag is empty)."""
    out = Path(root) / ("output" if tag == "" else f"run_{tag}")
    try:
        if out.exists():
            rmtree(out)
    except Exception:
        pass
    out.mkdir(parents=True, exist_ok=True)
    if not out.exists():
        out.mkdir(parents=True, exist_ok=True)
    print(f"[OUT DIR] {out.resolve()}")   # Optional breadcrumb
    return out

# Recreate directory if deleted mid-run
def _ensure_dir(p):
    Path(p).mkdir(parents=True, exist_ok=True)

# -----------------------------------------------------------------------------
def _ws_key_from_dem(name: str) -> str:
    nm = name.lower()
    return "LC3" if "lc3" in nm else "LC1"

def _choose_sampled_links_for_grid(g, dem_name: str):
    key = _ws_key_from_dem(dem_name)
    picks = SAMPLED_LINKS.get(key, {})
    return {s: np.asarray(picks.get(s, []), dtype=int) for s in SECTS}

def _validate_links(g, picks_by_sec, out_dir):
    valid_max = g.number_of_links
    cleaned = {}
    out_dir = Path(out_dir)
    with open(out_dir / "sampled_links_validation.txt", "w", encoding="utf-8") as fh:
        fh.write(f"# number_of_links = {valid_max}\n")
        for sec, arr in picks_by_sec.items():
            arr = np.asarray(arr, int)
            ok = arr[(arr >= 0) & (arr < valid_max)]
            bad = arr[(arr < 0) | (arr >= valid_max)]
            cleaned[sec] = ok
            fh.write(f"{sec}: kept {ok.tolist()} ; dropped {bad.tolist()}\n")
    return cleaned

# -----------------------------------------------------------------------------
# Convert a fraction vector + size bins to a Dx percentile (ascending sizes).
def _Dx(vec, sizes, p):
    n = min(len(vec), len(sizes))
    v = np.asarray(vec[:n]); s = np.asarray(sizes[:n])
    if p >= 1:
        return float(s[-1])
    if v.sum() == 0:
        return float(s[0])
    return float(np.interp(p, np.cumsum(v)/v.sum(), s))

# -----------------------------------------------------------------------------
# Weighted quantiles (I use this to summarize moved-size distributions by mass).
def _weighted_quantile(x, w, q):
    x = np.asarray(x, float); w = np.asarray(w, float)
    m = w > 0
    if not np.any(m):
        return np.nan if np.isscalar(q) else np.full_like(q, np.nan, dtype=float)
    x, w = x[m], w[m]
    order = np.argsort(x)
    x, w = x[order], w[order]
    cum_w = np.cumsum(w); total = cum_w[-1]
    probs = (cum_w - 0.5*w) / total
    return np.interp(q, probs, x)

# -----------------------------------------------------------------------------
# Online (one-pass) stats with optional light subsampling for percentiles.
class OnlineStats:
    def __init__(self, sf=0.01):
        self.n=0; self.mean=0.; self.M2=0.; self.vmin=math.inf; self.vmax=-math.inf
        self.sf=sf; self.sample=[]
    def add(self, x):
        if not (np.isfinite(x)):  # ignore non-finite
            return
        self.n+=1; d=x-self.mean; self.mean+=d/self.n; self.M2+=d*(x-self.mean)
        self.vmin=min(self.vmin,x); self.vmax=max(self.vmax,x)
        if random.random()<self.sf:
            self.sample.append(x)
    def summary(self):
        var=self.M2/(self.n-1) if self.n>1 else 0.
        p=(0.,0.,0.)
        if self.sample:
            p=np.percentile(self.sample,[5,50,95])
        return dict(min=self.vmin,max=self.vmax,mean=self.mean,
                    median=p[1],p5=p[0],p95=p[2],std=math.sqrt(var))

# -----------------------------------------------------------------------------
def _grid(dem, nodata, gsd_loc):
    g, z = read_esri_ascii(dem, name="topographic__elevation")
    nodemask = (z == nodata)
    g.add_field("valid_mask", (~nodemask).astype(np.uint8), at="node")
    g.at_node["topographic__elevation"] = masked_equal(z, nodata).filled(z.mean())

    g.add_zeros("surface_water__depth", at="node")
    g.add_zeros("surface_water__velocity", at="link")
    g.add_zeros("surface_water__discharge", at="link")
    g.at_node["surface_water__depth"] += 1e-3

    _, loc = read_esri_ascii(gsd_loc, name="bed_surf__gsd_loc_node")
    loc = np.asarray(loc, int)
    loc[loc < 0] = 0
    loc[loc > 2] = 2
    g.add_field("bed_surf__gsd_loc_node", loc, at="node", clobber=True)
    g.add_zeros("max_size_moved_mm", at="link")
    return g

# -----------------------------------------------------------------------------
def _labels_from_gsd_loc(g):
    code_to_name = {0: "upstream", 1: "intermediate", 2: "downstream"}
    node_codes = np.asarray(g.at_node["bed_surf__gsd_loc_node"], int)
    node_lab = np.array([code_to_name.get(int(c), "downstream") for c in node_codes], dtype="<U12")
    z = g.at_node["topographic__elevation"]
    n0 = g.nodes_at_link[:, 0]; n1 = g.nodes_at_link[:, 1]
    dwn = np.where(z[n0] <= z[n1], n0, n1)
    link_codes = node_codes[dwn]
    link_lab = np.array([code_to_name.get(int(c), "downstream") for c in link_codes], dtype="<U12")
    return node_lab, link_lab, link_codes

# ─────────────────────────────────────────────────────────────────────────────
def run_simulation(
        z_dem=DEM_NAME, rainfall_file="working_climate/15min_1000yr.txt",
        gsd_file=GSD_TABLE, gsd_loc_raster=GSD_LOC_RASTER,
        nodata=-9999, max_dt=0.5,
        sim_max_t=25000,
        cwd=getcwd(), tag="", plots=False):

    cwd = Path(cwd)
    rain_stem = Path(rainfall_file).stem
    run_tag = tag if (tag and tag != rain_stem) else rain_stem
    out = _clean(cwd, run_tag)

    g   = _grid(cwd/z_dem, nodata, cwd/gsd_loc_raster)

    # Quick peek at the section-class raster
    if "bed_surf__gsd_loc_node" in g.at_node:
        _vals = np.asarray(g.at_node["bed_surf__gsd_loc_node"]).astype(int)
        _uniq, _cnts = np.unique(_vals, return_counts=True)
        print("\n[DEBUG] bed_surf__gsd_loc_node — initial summary")
        print(f"  shape={_vals.shape}, dtype={_vals.dtype}")
        print(f"  unique={_uniq.tolist()}, counts={_cnts.tolist()}")
        print(f"  min={_vals.min()}, max={_vals.max()}, n_zero={int((_vals==0).sum())}, n_nonzero={int((_vals!=0).sum())}")
        print(f"  head20={_vals[:20].tolist()}")
    else:
        print("\n[DEBUG] bed_surf__gsd_loc_node field NOT FOUND in g.at_node")

    of = OverlandFlow(g, mannings_n=0.025, alpha=0.25,
                      rainfall_intensity=0., steep_slopes=True, theta=1.0)
    FlowAccumulator(g, flow_director=FlowDirectorSteepest(g)).run_one_step()

    fixed = np.zeros(g.number_of_nodes, int)
    for idx in (300, 301, 302):
        if idx < g.number_of_nodes:
            fixed[idx] = 1

    raw = np.loadtxt(cwd/gsd_file)
    sizes_mm_col = raw[:, 0].copy()
    site0 = raw[:, 1].copy()
    site1 = raw[:, 2].copy()
    site2 = raw[:, 3].copy()
    raw_gsd = np.column_stack([sizes_mm_col, site0, site1, site2])

    node_lab, link_lab, link_class = _labels_from_gsd_loc(g)
    sec_counts = {s: int((node_lab == s).sum()) for s in SECTS}
    print("[Section node counts (from raster)]", sec_counts)
    g.add_field("bed_surf__gsd_loc_link", link_class.astype(int), at="link", overwrite=True)

    rbd = RiverBedDynamics(
        g, gsd=raw_gsd,
        outlet_boundary_condition="fixedValue",
        bed_surf__elev_fix_node=fixed,
        bedload_equation=BEDLOAD_EQUATION,
        use_hydraulics_radius_in_shear_stress=bool(USE_RH_TAU),
        track_stratigraphy=bool(TRACK_STRATIGRAPHY),
        bed_surf__gsd_loc_node=g.at_node["bed_surf__gsd_loc_node"],
    )

    assert "topographic__elevation" in g.at_node
    assert "surface_water__depth"   in g.at_node
    assert "surface_water__discharge" in g.at_link
    assert "drainage_area" in g.at_node

    SIZE_MM = raw_gsd[:, 0]
    sizes_descending = np.all(np.diff(SIZE_MM) < 0)
    if sizes_descending:
        SIZE_MM = SIZE_MM[::-1]
    if SEED_SURFACE:
        rbd._bed_surf__gsd_link[:] = 1/len(SIZE_MM)

    # thresholds from table
    try:
        gsd_view   = raw_gsd[::-1, :] if sizes_descending else raw_gsd.copy()
        sizes_view = gsd_view[:, 0]
        site_cols  = {"upstream": 1, "intermediate": 2, "downstream": 3}
        SECTION_THRESHOLDS_FROM_TABLE = {}
        eps = 1e-9
        for sec in SECTS:
            col_idx = site_cols[sec]
            cdf = gsd_view[:, col_idx].astype(float)
            idxs = np.where(cdf < (100.0 - eps))[0]
            max_any = float(sizes_view[idxs[-1]]) if idxs.size else float(sizes_view[-1])
            SECTION_THRESHOLDS_FROM_TABLE[sec] = {
                "max_input_any_mm": max_any,
                "max_eligible_FRAC_MIN_mm": max_any,
            }
        with open(out / "gsd_section_thresholds_from_table.txt", "w", encoding="utf-8") as fth:
            fth.write("section max_input_any_mm max_eligible_FRAC_MIN_mm\n")
            for sec in SECTS:
                th = SECTION_THRESHOLDS_FROM_TABLE.get(sec, {"max_input_any_mm": float("nan"),
                                                             "max_eligible_FRAC_MIN_mm": float("nan")})
                fth.write(f"{sec} {th['max_input_any_mm']:.6f} {th['max_eligible_FRAC_MIN_mm']:.6f}\n")
    except Exception:
        SECTION_THRESHOLDS_FROM_TABLE = {}
        with open(out / "gsd_section_thresholds_from_table.txt", "w", encoding="utf-8") as fth:
            fth.write("# NOTE: could not compute thresholds directly from the table; proceeding without.\n")

    QB_EPS = 1e-10

    def _D50_from_frac(frac_row: np.ndarray, sizes_mm: np.ndarray) -> float:
        f = np.clip(frac_row, 0.0, None)
        s = f.sum()
        if s <= 0:
            return float(sizes_mm[0])
        cdf = np.cumsum(f) / s
        return float(np.interp(0.5, cdf, sizes_mm[:len(f)]))

    def _is_cumulative(arr_row: np.ndarray) -> bool:
        if arr_row.size == 0: return False
        if np.any(arr_row < -1e-12): return False
        nondec = np.all(np.diff(arr_row) >= -1e-9)
        ends_one = 0.95 <= arr_row[-1] <= 1.05 or 95 <= arr_row[-1] <= 105
        starts_ge0 = arr_row[0] >= -1e-9
        return bool(nondec and ends_one and starts_ge0)

    def link_frac() -> np.ndarray:
        arr = rbd._bed_surf__gsd_link
        if sizes_descending:
            arr = arr[:, ::-1]
        sample = arr[0] if arr.shape[0] else np.array([])
        arr = np.asarray(arr, float)
        arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
        if _is_cumulative(sample):
            if sample[-1] > 1.5:
                arr = arr / 100.0
            cum = arr
            fr = np.diff(np.c_[np.zeros((cum.shape[0], 1)), cum], axis=1)
            fr = np.clip(fr, 0.0, None)
        else:
            fr = np.clip(arr, 0.0, None)
        s = fr.sum(axis=1, keepdims=True)
        s[s <= 0.0] = 1.0
        return fr / s

    gsd0_link = link_frac().copy()
    gsd_prev  = gsd0_link.copy()

    rain_raw = np.loadtxt(cwd / rainfall_file)
    rain = rain_raw[np.all(np.isfinite(rain_raw), axis=1)]
    if rain.size == 0:
        raise ValueError(f"Rainfall file {rainfall_file} has no finite rows.")
    rtime, rint, ridx, curr_r, rval = rainfall_data(rain, g)

    t_end = float(np.nanmax(rtime)) if np.size(rtime) else float(sim_max_t)
    t_end += 180.0
    t_end = float(min(t_end, sim_max_t))

    z0 = g.at_node["topographic__elevation"].copy()

    n0 = g.nodes_at_link[:, 0]; n1 = g.nodes_at_link[:, 1]
    elev = g.at_node["topographic__elevation"]
    dwn = np.where(elev[n0] <= elev[n1], n0, n1)
    dwn_sec = np.array([node_lab[i] for i in dwn], dtype="<U12")
    export_link = np.array([(link_lab[i]!="downstream") and (dwn_sec[i]!=link_lab[i])
                            for i in range(g.number_of_links)])

    PRES_STRICT = 0.002
    channel_links = _choose_sampled_links_for_grid(g, z_dem)
    channel_links = _validate_links(g, channel_links, out)

    # Enforce that requested channel links actually belong to their declared section
    def _enforce_channel_links_section(link_lab, picks_by_sec, out_dir):
        kept = {s: [] for s in SECTS}
        out_dir = Path(out_dir)
        with open(out_dir / "channel_links_section_check.txt", "w", encoding="utf-8") as fh:
            for sec in SECTS:
                for lid in picks_by_sec.get(sec, []):
                    true_sec = link_lab[int(lid)]
                    if true_sec == sec:
                        kept[sec].append(int(lid))
                    else:
                        fh.write(f"dropped link {lid} (intended {sec}, actual {true_sec})\n")
            for sec in SECTS:
                fh.write(f"kept {sec}: {kept[sec]}\n")
        return {s: np.asarray(v, dtype=int) for s, v in kept.items()}

    channel_links = _enforce_channel_links_section(link_lab, channel_links, out)

    # Helper to save a single link’s GSD snapshot (both fractions and CDF) for before/after
    def _write_link_gsd_snapshot(path, sizes_mm, frac_row):
        f = np.clip(frac_row, 0.0, None)
        s = f.sum()
        if s > 0: f = f / s
        c = np.cumsum(f)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("# size_mm  frac  cdf\n")
            for sz, ff, cc in zip(sizes_mm[:len(f)], f, c):
                fh.write(f"{sz:8.3f} {ff:9.6f} {cc:9.6f}\n")

    # Save BEFORE snapshots for every selected channel link
    gsd0_link_for_probe = gsd0_link.copy()
    for sec in SECTS:
        for lid in channel_links[sec]:
            _write_link_gsd_snapshot(out / f"channel_probe_{sec}_link{lid}_before.txt",
                                     SIZE_MM, gsd0_link_for_probe[lid, :])

    # --- Time-series collectors for selected channel links --------------------
    chan_ts = {s: {int(lid): {"t": [], "tau": [], "qb": [], "depth": [], "D50": []}
                   for lid in channel_links[s]} for s in SECTS}

    # Section-level accumulators
    ts      ={s:[] for s in SECTS}            # per-step largest moved size (mm) after warm-up
    ts_moved={s:[] for s in SECTS}            # per-step union of moved sizes
    lag     ={s:None for s in SECTS}; tlag=lag.copy()
    mass    ={s:None for s in SECTS}          # integrated mass moved by size (per section)
    dur     ={s:{} for s in SECTS}            # per-size cumulative duration of motion
    tau_stats={s:OnlineStats() for s in SECTS}
    qb_stats ={s:OnlineStats() for s in SECTS}
    tau_stats_wet={s:OnlineStats() for s in SECTS}
    qb_stats_wet ={s:OnlineStats() for s in SECTS}

    # Per-link streaming stats and transport integral (∫qb dt per link)
    link_tau_stats = [OnlineStats(sf=0.0) for _ in range(g.number_of_links)]
    link_qb_stats  = [OnlineStats(sf=0.0) for _ in range(g.number_of_links)]
    link_total_qb  = np.zeros(g.number_of_links)  # ∫ qb dt (m^2) per link

    # Run-long criteria: (a) moved on ≥2 links in section; (b) exported across a boundary
    two_links_mass = {s: None for s in SECTS}
    export_mass    = {s: None for s in SECTS}

    # RBD expects previous-step velocities; sustain counters allocated later
    rbd._surface_water__velocity_prev_time_link = np.zeros(g.number_of_links)
    moving_steps = None  # shape (n_links, n_sizes), integers

    # ==== T=0 DIAGNOSTICS with STRICT presence threshold ====
    def _max_size_where(fr_mat, sizes, thresh):
        """Largest size with fraction > thresh at any link in section."""
        if fr_mat.size == 0:
            return float("nan")
        present = (fr_mat > thresh).any(axis=0)
        return float(sizes[np.where(present)[0].max()]) if present.any() else float("nan")

    with open(out / "gsd_init_diagnostics.txt","w", encoding="utf-8") as fdiag:
        fdiag.write("# Initialized surface GSD (t=0) diagnostics\n")
        fdiag.write(f"# any>0, >FRAC_MIN ({FRAC_MIN:g}), >PRES_STRICT ({PRES_STRICT:g})\n")

        section_lids = {s: np.where(link_lab == s)[0] for s in SECTS}
        for sec in SECTS:
            lids = section_lids[sec]
            if lids.size == 0:
                fdiag.write(f"{sec}: no links\n"); continue
            fr_here = gsd0_link[lids]
            max_any   = _max_size_where(fr_here, SIZE_MM, 0.0)
            max_fmin  = _max_size_where(fr_here, SIZE_MM, FRAC_MIN)
            max_strct = _max_size_where(fr_here, SIZE_MM, PRES_STRICT)
            fdiag.write(
                f"{sec}: any>0 {max_any:.1f} mm, >FRAC_MIN {max_fmin:.1f} mm, "
                f">PRES_STRICT {max_strct:.1f} mm\n"
            )

    # Eligibility mask at t=0 using STRICT presence
    init_mask = (gsd0_link > PRES_STRICT)

    # For reporting, capture which sizes are present anywhere within each section at t=0
    init_present_reach = {}
    for sec in SECTS:
        lids = np.where(link_lab == sec)[0]
        present = (gsd0_link[lids] > PRES_STRICT).any(axis=0) if lids.size else np.zeros(gsd0_link.shape[1], bool)
        init_present_reach[sec] = present

    # ---- Main time loop ------------------------------------------------------
    t = 0.0
    while t < t_end:
        if t == 0.0:
            of.dt = max_dt
            t += of.dt
            continue

        # Update rainfall and hydraulics
        ridx, curr_r, rval, _ = update_rainfall_intensity(t, rtime, ridx, rint, g, curr_r, rval, False)
        of._rainfall_intensity = curr_r
        of.overland_flow(dt=max_dt)
        if not np.isfinite(of.dt) or of.dt <= 0:
            of.dt = max_dt * 1e-3

        # Dry-tail fast-forward
        if (g.at_node["surface_water__depth"].max() < 1e-9) and (abs(curr_r) < 1e-12):
            t = min(t + max(10.0, of.dt * 20.0), t_end)
            continue

        # Only run RBD when there’s enough water around
        if g.at_node["surface_water__depth"].max() >= DEPTH_MIN:
            depth_link = g.map_mean_of_link_nodes_to_link(g.at_node["surface_water__depth"])
            depth_link = np.maximum(depth_link, 1e-9)
            g["link"]["surface_water__velocity"] = g["link"]["surface_water__discharge"] / depth_link

            rbd._surface_water__velocity_prev_time_link = g["link"]["surface_water__velocity"]
            rbd._grid._dt = of.dt
            try:
                rbd.run_one_step()
            except ValueError as e:
                if "need at least one array to concatenate" in str(e):
                    warnings.warn("No stratigraphy layers this step; skipping RBD for this dt.")
                else:
                    raise

            # Pull updated fields from RBD and OF
            gsd_prev = link_frac().copy()
            qb_link  = np.abs(rbd._sed_transp__bedload_rate_link)
            tau_link = np.abs(rbd._surface_water__shear_stress_link)
            n_now    = gsd_prev.shape[1]
            sizes_now= SIZE_MM[:n_now]
            wet_link = (depth_link >= DEPTH_MIN)

            # Per-step section summaries
            sec_max_this_step = {s: 0.0 for s in SECTS}
            sec_union_this_step = {s: [] for s in SECTS}

            # Per-step counters and mass arrays for my two reporting criteria
            counts_2links = {s: np.zeros(n_now, dtype=int) for s in SECTS}
            mass_tmp_2    = {s: np.zeros(n_now)           for s in SECTS}
            counts_export = {s: np.zeros(n_now, dtype=int) for s in SECTS}
            mass_tmp_exp  = {s: np.zeros(n_now)            for s in SECTS}

            # Allocate run-long accumulators once sizes are known
            if mass[SECTS[0]] is None:
                mass = {s: np.zeros(n_now) for s in SECTS}
                moving_steps = np.zeros((g.number_of_links, n_now), dtype=int)
                two_links_mass = {s: np.zeros(n_now) for s in SECTS}
                export_mass    = {s: np.zeros(n_now) for s in SECTS}

            # Append placeholders for this step
            for s in SECTS:
                ts[s].append(0.0)
                ts_moved[s].append(np.array([], dtype=float))

            # Loop over links to accumulate stats and classify motion by size
            for lid in range(g.number_of_links):
                sec = link_lab[lid]

                # Aggregate τ and qb stats
                tau_stats[sec].add(tau_link[lid])
                qb_stats [sec].add(qb_link[lid])
                if wet_link[lid]:
                    tau_stats_wet[sec].add(tau_link[lid])
                    qb_stats_wet [sec].add(qb_link[lid])

                # Per-link stats and transport integral
                link_tau_stats[lid].add(tau_link[lid])
                link_qb_stats [lid].add(qb_link[lid])
                link_total_qb [lid] += qb_link[lid] * of.dt  # m^2

                # Log time series at selected channel links
                if lid in chan_ts.get(sec, {}):
                    D50_here = _D50_from_frac(gsd_prev[lid, :n_now], sizes_now)
                    rec = chan_ts[sec][int(lid)]
                    rec["t"].append(float(t))
                    rec["tau"].append(float(tau_link[lid]))
                    rec["qb"].append(float(qb_link[lid]))
                    rec["depth"].append(float(depth_link[lid]))
                    rec["D50"].append(float(D50_here))

                # Size-resolved motion with sustain criterion:
                frac      = gsd_prev[lid, :n_now]
                allowed   = init_mask[lid][:n_now]
                qb_by_size= qb_link[lid] * frac
                moving_now= allowed & (qb_by_size > QB_EPS)
                moving_steps[lid, moving_now] += 1
                moving_steps[lid, ~moving_now]  = 0

                moved_idx = (np.where(moving_steps[lid, :] >= int(SUSTAIN_STEPS))[0]
                             if SUSTAIN_STEPS > 1 else np.where(moving_now)[0])
                if moved_idx.size == 0:
                    continue

                # Integrate transported mass by size at the section level
                inc = qb_by_size[moved_idx] * of.dt
                mass[sec][moved_idx] += inc

                # Update per-step tallies for the two criteria
                counts_2links[sec][moved_idx] += 1
                mass_tmp_2[sec][moved_idx]    += inc
                if export_link[lid]:
                    counts_export[sec][moved_idx] += 1
                    mass_tmp_exp[sec][moved_idx]  += inc

                # After warm-up, update section-wide max size and union of moved sizes
                if t > float(WARMUP_S):
                    moved_sizes = sizes_now[moved_idx]
                    sec_max_this_step[sec] = max(sec_max_this_step[sec], float(moved_sizes.max()))
                    sec_union_this_step[sec].append(moved_sizes)
                    for sz in moved_sizes:
                        dur[sec][float(sz)] = dur[sec].get(float(sz), 0.0) + of.dt

            # Fill per-step section series once warm-up is over
            if t > float(WARMUP_S):
                for s in SECTS:
                    if sec_max_this_step[s] > 0.0:
                        ts[s][-1] = sec_max_this_step[s]
                    if sec_union_this_step[s]:
                        ts_moved[s][-1] = np.unique(np.concatenate(sec_union_this_step[s]))

            # Commit per-step masses for the two criteria into the run-long arrays
            for s in SECTS:
                sel2 = counts_2links[s] >= 2
                if np.any(sel2):
                    two_links_mass[s][:len(mass_tmp_2[s])] += mass_tmp_2[s] * sel2
                sele = counts_export[s] >= 1
                if np.any(sele):
                    export_mass[s][:len(mass_tmp_exp[s])] += mass_tmp_exp[s] * sele

            # Update per-link “largest size moved” map after warm-up
            if t > float(WARMUP_S):
                for s in SECTS:
                    if sec_max_this_step[s] > 0.0:
                        lids = np.where(link_lab==s)[0]
                        for lid in lids:
                            g.at_link["max_size_moved_mm"][lid] = max(
                                g.at_link["max_size_moved_mm"][lid], sec_max_this_step[s]
                            )

        else:
            # too shallow—append placeholders
            for s in SECTS:
                ts[s].append(0.0)
                ts_moved[s].append(np.array([], dtype=float))

        # Advance time robustly
        t_next = update_time(t, of, t_end, False, False, 0)[0]
        t = float(t_next) if np.isfinite(t_next) else t + float(of.dt or 0.0)

    # End of main loop
    gsd1_link = gsd_prev

    # ── POST-PROCESS RASTERS & STATS ─────────────────────────────────────────
    _ensure_dir(out)  # ← EDIT: ensure run_* directory still exists
    dz = g.at_node["topographic__elevation"] - z0
    g.add_field("elev_change", dz, at="node", overwrite=True)
    write_esri_ascii(out/"elev_change.asc", g, names="elev_change", clobber=True)

    def _Dx_map(before, after, pct, name):
        """Compute ΔDx at links, map to nodes by mean of incident links, then write ESRI ASCII."""
        D0 = np.array([_Dx(b, SIZE_MM, pct) for b in before])
        D1 = np.array([_Dx(a, SIZE_MM, pct) for a in after])
        d  = D1 - D0
        node = np.zeros(g.number_of_nodes); cnt = np.zeros_like(node)
        for lid, val in enumerate(d):
            n1, n2 = g.nodes_at_link[lid]
            node[n1] += val; node[n2] += val
            cnt[n1]  += 1;   cnt[n2]  += 1
        cnt[cnt==0] = 1
        node /= cnt
        fn = f"{name}_change_mm_node"
        g.add_field(fn, node, at="node", overwrite=True)
        write_esri_ascii(out/f"{name}_change.asc", g, names=fn, clobber=True)

    # Export ΔD10, ΔD50, ΔD90, and ΔDmax rasters
    for pct, nm in [(0.10,"D10"),(0.50,"D50"),(0.90,"D90"),(1.00,"Dmax")]:
        _Dx_map(gsd0_link, gsd1_link, pct, nm)

    # --- AFTER/BETWEEN: snapshots & deltas (text files) -----------------------
    def _norm_frac(row):
        f = np.clip(np.array(row, float), 0.0, None)
        s = f.sum()
        return f / s if s > 0 else f

    # 1) Write AFTER snapshots for every selected channel link
    for sec in SECTS:
        for lid in channel_links[sec]:
            _write_link_gsd_snapshot(out / f"channel_probe_{sec}_link{int(lid)}_after.txt",
                                     SIZE_MM, gsd1_link[int(lid), :])

    # 2) Per-link DELTA GSDs in both signs
    for sec in SECTS:
        for lid in channel_links[sec]:
            lid = int(lid)
            b = _norm_frac(gsd0_link_for_probe[lid, :])
            a = _norm_frac(gsd1_link[lid, :])
            n = min(len(b), len(a), len(SIZE_MM))
            sizes = SIZE_MM[:n]
            diff_a_minus_b = a[:n] - b[:n]
            diff_b_minus_a = b[:n] - a[:n]
            with open(out / f"channel_probe_{sec}_link{lid}_diff.txt", "w", encoding="utf-8") as f:
                f.write("# size_mm  Δfrac(after - before)\n")
                for sz, df in zip(sizes, diff_a_minus_b):
                    f.write(f"{sz:8.3f} {df:10.6f}\n")
            with open(out / f"channel_probe_{sec}_link{lid}_diff_b_minus_a.txt", "w", encoding="utf-8") as f:
                f.write("# size_mm  Δfrac(before - after)\n")
                for sz, df in zip(sizes, diff_b_minus_a):
                    f.write(f"{sz:8.3f} {df:10.6f}\n")

    # 3) Section-average BEFORE/AFTER and DELTAs (both signs)
    sec_mask = {s: (link_lab == s) for s in SECTS}
    def _section_mean_frac(frac_matrix, mask):
        if np.sum(mask) == 0:
            return np.zeros(frac_matrix.shape[1])
        m = np.clip(frac_matrix[mask], 0.0, None)
        m = np.nan_to_num(m, nan=0.0, posinf=0.0, neginf=0.0)
        f = m.mean(axis=0)
        s = f.sum()
        return f / s if s > 0 else f

    sec_before = {s: _section_mean_frac(gsd0_link, sec_mask[s]) for s in SECTS}
    sec_after  = {s: _section_mean_frac(gsd1_link, sec_mask[s]) for s in SECTS}

    for s in SECTS:
        f0 = _norm_frac(sec_before[s])
        f1 = _norm_frac(sec_after[s])
        n = min(len(f0), len(f1), len(SIZE_MM))
        sizes = SIZE_MM[:n]
        # write section before & after as 3-column (size, frac, cdf)
        def _write_sec_snapshot(path, frac_row):
            frac_row = _norm_frac(frac_row)[:n]
            cdf = np.cumsum(frac_row)
            with open(path, "w", encoding="utf-8") as fh:
                fh.write("# size_mm  frac  cdf\n")
                for sz, ff, cc in zip(sizes, frac_row, cdf):
                    fh.write(f"{sz:8.3f} {ff:9.6f} {cc:9.6f}\n")
        _write_sec_snapshot(out / f"section_{s}_gsd_before.txt", f0)
        _write_sec_snapshot(out / f"section_{s}_gsd_after.txt",  f1)

        # deltas
        fd_a_minus_b = f1[:n] - f0[:n]
        fd_b_minus_a = f0[:n] - f1[:n]
        with open(out / f"section_{s}_diff.txt", "w", encoding="utf-8") as fh:
            fh.write("# size_mm  Δfrac(after - before)\n")
            for sz, df in zip(sizes, fd_a_minus_b):
                fh.write(f"{sz:8.3f} {df:10.6f}\n")
        with open(out / f"section_{s}_diff_b_minus_a.txt", "w", encoding="utf-8") as fh:
            fh.write("# size_mm  Δfrac(before - after)\n")
            for sz, df in zip(sizes, fd_b_minus_a):
                fh.write(f"{sz:8.3f} {df:10.6f}\n")

    # 4) Channel-probe summary TXT (means) and per-link time series TXTs
    rows = []
    for sec in SECTS:
        for lid, rec in chan_ts[sec].items():
            ts_path = out / f"channel_probe_{sec}_link{int(lid)}_timeseries.txt"
            if len(rec["t"]):
                with open(ts_path, "w", encoding="utf-8") as fh:
                    fh.write("t_s depth_m tau_Pa qb_m2s D50_mm\n")
                    for T, D, TAU, QB, D50 in zip(rec["t"], rec["depth"], rec["tau"], rec["qb"], rec["D50"]):
                        fh.write(f"{T:.6f} {D:.9g} {TAU:.9g} {QB:.9g} {D50:.9g}\n")
            else:
                with open(ts_path, "w", encoding="utf-8") as fh:
                    fh.write("t_s depth_m tau_Pa qb_m2s D50_mm\n")

            def _nanmean(x):
                x = np.asarray(x, float)
                return float(np.nanmean(x)) if x.size else float("nan")
            rows.append(dict(
                section=sec,
                link_id=int(lid),
                t_mean_s=_nanmean(rec["t"]),
                depth_mean_m=_nanmean(rec["depth"]),
                tau_mean_Pa=_nanmean(rec["tau"]),
                qb_mean_m2s=_nanmean(rec["qb"]),
                D50_mean_mm=_nanmean(rec["D50"]),
            ))

    def _fmt(x, fp=":.9g"):
        try:
            xv = float(x)
        except Exception:
            return "nan"
        if not np.isfinite(xv):
            return "nan"
        spec = fp[1:] if (isinstance(fp, str) and fp.startswith(":")) else fp
        try:
            return format(xv, spec)
        except Exception:
            return f"{xv:.9g}"

    with open(out / "channel_probe_summary.txt", "w", encoding="utf-8") as fh:
        fh.write("section link_id t_mean_s depth_mean_m tau_mean_Pa qb_mean_m2s D50_mean_mm\n")
        for r in rows:
            fh.write(
                f"{r['section']} {int(r['link_id'])} "
                f"{_fmt(r['t_mean_s'], ':.6f')} "
                f"{_fmt(r['depth_mean_m'])} "
                f"{_fmt(r['tau_mean_Pa'])} "
                f"{_fmt(r['qb_mean_m2s'])} "
                f"{_fmt(r['D50_mean_mm'])}\n"
            )

    # 5) Largest moved grain summaries per section (timeseries + final) → TXT
    for s in SECTS:
        ts_path = out / f"section_{s}_largest_moved_mm_timeseries.txt"
        with open(ts_path, "w", encoding="utf-8") as fh:
            fh.write("step_index max_size_mm\n")
            for k, val in enumerate(ts[s]):
                fh.write(f"{k} {float(val):.6f}\n")
        final_max = float(np.nanmax(ts[s])) if len(ts[s]) else float("nan")
        with open(out / f"section_{s}_largest_moved_mm_final.txt", "w", encoding="utf-8") as fh:
            fh.write(f"{final_max:.6f}\n")

    # ── Run-health JSON ───────────────────────────────────────────────────────
    health = {
        "storm_file": rain_stem,
        "t_end_s": float(t_end),
        "max_depth_m": float(np.max(g.at_node["surface_water__depth"])),
        "any_rbd_steps": bool(np.isfinite(np.sum(link_total_qb)) and np.sum(link_total_qb) > 0),
        "sections_with_links": {
            s: int(np.sum(g.at_link["bed_surf__gsd_loc_link"]==i)) for i,s in enumerate(SECTS)
        },
    }
    with open(out / "run_health.json", "w", encoding="utf-8") as fh:
        json.dump(health, fh, indent=2)
    print("[RUN HEALTH]", health)

    # Return output directory for downstream scripts
    return str(out)