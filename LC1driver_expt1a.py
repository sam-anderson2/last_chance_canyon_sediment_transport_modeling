# -*- coding: utf-8 -*-
"""
driver_expt1a  –  v6-patch (“PNG-safe”, streaming stats, no RAM blow-up)
========================================================================
Single-event driver that couples:

    • Landlab.OverlandFlow      – hydraulics on a raster grid
    • Landlab.RiverBedDynamics  – mixed-grain sediment transport
    • (optional) simple plotting / ASCII raster export

Outputs written per run
-----------------------
run_<tag>/  (or ./output/ if tag == "")
    elev_change.asc / .png               Δz node raster
    D10_change.asc  / .png               ΔD₁₀ node raster
    D50_change.asc  / .png               ΔD₅₀ node raster
    D90_change.asc  / .png               ΔD₉₀ node raster
    Dmax_change.asc / .png               ΔD_max node raster
    max_size_moved.asc                   largest clast moved (node)
    largest_size_timeseries.png          per-band timeline of max clast
    shear_stats_<band>.txt               streaming |τ| stats
    bedload_stats_<band>.txt             streaming |qb| stats
    elev_change_stats_<band>.txt         summary stats of Δz
    surf_gsd_change_<band>.txt           reach-avg surface GSD shift
    summary.txt                          concise human-readable recap

Elevation-band logic (LC-1 vs LC-3)
-----------------------------------
The helper `_bands()` inspects the DEM filename:

    LC-3 DEM name contains “lc3” →  (1560–∞, 1530–1560, <1530) m
    otherwise                     →  (1560–∞, 1540–1560, <1540) m
    
---------------------------------------------------------------------
Folder layout expected by default
-----------------------------------------------------------------------
project_root/
│
├─ driver_expt1a.py                ← THIS FILE
├─ run_many_storms.py
│
├─ filled_lc1_dem.asc              ← LC-1 example DEM
├─ filled_lc3_dem.asc              ← LC-3 example DEM
├─ lc1_gsd_locations.asc           ← raster: which GSD applies where
├─ lc3_gsd_locations.asc
│
├─ GSDs/
│   ├─ LC1_grain_size_dist_with_max.txt
│   └─ LC3_grain_size_dist_with_max.txt
│
└─ working_climate/
    ├─ 15min_1000yr.txt
    ├─ 6hr_25yr.txt
    └─ …


If you swap watersheds or change thresholds, that’s the **only** place to edit.

Author: Sam Anderson
------

"""

# driver_expt1a.py — OverlandFlow + RiverBedDynamics driver (Windows-safe)

import matplotlib, argparse, math, random, warnings
matplotlib.use("Agg")                     # PNGs inside worker procs

from pathlib import Path
from shutil  import rmtree
from os      import getcwd, mkdir
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image  as mpimg
from numpy.ma import masked_equal

from landlab.io        import read_esri_ascii, write_esri_ascii
from landlab.components import (
    OverlandFlow, RiverBedDynamics, FlowAccumulator, FlowDirectorSteepest,
)
from rainfall_manager import rainfall_data, update_rainfall_intensity
from update_time      import update_time

warnings.filterwarnings("ignore", category=UserWarning)

# ── constants & tiny helpers ───────────────────────────────────────────────
EPS_MASS, DEPTH_MIN = 1e-3, 0.002
SECTS = ("upstream", "intermediate", "downstream")

def _clean(root: Path, tag: str) -> Path:
    out = root / ("output" if tag == "" else f"run_{tag}")
    if out.exists(): rmtree(out)
    mkdir(out); return out

def _bands(dem):
    return ((1560,np.inf),(1530,1560),(0,1530)) if "lc3" in dem.lower() \
           else ((1560,np.inf),(1540,1560),(0,1540))

def _label_links(g,rngs):
    up,mid,dn=rngs
    lab=np.full(g.number_of_links,"none","<U12")
    em=0.5*(g.at_node["topographic__elevation"][g.nodes_at_link[:,0]]+
            g.at_node["topographic__elevation"][g.nodes_at_link[:,1]])
    lab[(up [0]<=em)&(em<up [1])]="upstream"
    lab[(mid[0]<=em)&(em<mid[1])]="intermediate"
    lab[(dn [0]<=em)&(em<dn [1])]="downstream"
    return lab

def _label_nodes(g,rngs):
    up,mid,dn=rngs
    elev=g.at_node["topographic__elevation"]
    lab=np.full(g.number_of_nodes,"none","<U12")
    lab[(up [0]<=elev)&(elev<up [1])]="upstream"
    lab[(mid[0]<=elev)&(elev<mid[1])]="intermediate"
    lab[(dn [0]<=elev)&(elev<dn [1])]="downstream"
    return lab

def _Dx(vec, sizes, p):
    n=min(len(vec),len(sizes))
    v=np.asarray(vec[:n]); s=np.asarray(sizes[:n])
    if p>=1:           return float(s[-1])
    if v.sum()==0.0:   return float(s[0])
    cf=np.cumsum(v)/v.sum()
    return float(np.interp(p,cf,s))

def _safe_imsave(a,f,vcenter=0,cmap="RdBu"):
    m=np.nanmax(np.abs(a-vcenter))
    mpimg.imsave(f,a,cmap=cmap,vmin=vcenter-m,vmax=vcenter+m,format="png")

class OnlineStats:
    def __init__(self,sf=0.01):
        self.n=0; self.mean=0.; self.M2=0.; self.vmin=math.inf; self.vmax=-math.inf
        self.sf=sf; self.sample=[]
    def add(self,x):
        self.n+=1; d=x-self.mean
        self.mean+=d/self.n; self.M2+=d*(x-self.mean)
        self.vmin=min(self.vmin,x); self.vmax=max(self.vmax,x)
        if random.random()<self.sf: self.sample.append(x)
    def summary(self):
        var=self.M2/(self.n-1) if self.n>1 else 0.
        p5,p50,p95=(np.percentile(self.sample,[5,50,95])
                    if self.sample else (0,0,0))
        return dict(min=self.vmin,max=self.vmax,mean=self.mean,
                    median=p50,p5=p5,p95=p95,std=math.sqrt(var))

# ── grid builder ───────────────────────────────────────────────────────────
def _grid(dem,nodata,gsd_loc):
    g,z=read_esri_ascii(dem,name="topographic__elevation")
    g.at_node["topographic__elevation"]=masked_equal(z,nodata).filled(z.mean())
    g.add_zeros("surface_water__depth",at="node")
    g.add_zeros("surface_water__velocity",at="link")
    g.at_node["surface_water__depth"]+=1e-3
    _,locs=read_esri_ascii(gsd_loc,name="bed_surf__gsd_loc_node")
    g.add_field("node","bed_surf__gsd_loc_node",locs)
    g.add_zeros("link","max_size_moved_mm")
    return g

# ── main driver ────────────────────────────────────────────────────────────
def run_simulation(z_dem="filled_lc1_dem.asc",
                   rainfall_file="working_climate/15min_1000yr.txt",
                   gsd_file="GSDs/LC1_grain_size_dist_with_max.txt",
                   gsd_loc_raster="lc1_gsd_locations.asc",
                   nodata=-9999,max_dt=0.5,sim_max_t=1200,
                   cwd=getcwd(),tag="",plots=True):

    cwd=Path(cwd); out=_clean(cwd,tag)
    g=_grid(cwd/z_dem,nodata,cwd/gsd_loc_raster)

    of=OverlandFlow(g,mannings_n=0.025,alpha=0.25,
                    rainfall_intensity=0.,steep_slopes=True,theta=1.0)
    FlowAccumulator(g,flow_director=FlowDirectorSteepest(g)).run_one_step()

    fixed=np.zeros(g.number_of_nodes,int); fixed[[300,301,302]]=1

    # ── load raw cumulative GSD table (descending), init component ──
    raw_gsd=np.loadtxt(cwd/gsd_file)
    rbd=RiverBedDynamics(g,gsd=raw_gsd,
         outlet_boundary_condition="fixedValue",
         bed_surf__elev_fix_node=fixed,
         bedload_equation="WilcockAndCrowe")

    # ── convert to ascending size + per-class fractions for diagnostics ──
    sizes_raw=raw_gsd[:,0]
    descending=np.all(np.diff(sizes_raw)<0)
    cumlike =np.all(np.diff(raw_gsd[:,1])>=0) and raw_gsd[0,1]>90

    if descending and cumlike:
        SIZE_MM=sizes_raw[::-1]                     # ascending
        def link_frac():
            cum=rbd._bed_surf__gsd_link[:,::-1]     # flip cols
            return np.diff(np.concatenate([np.zeros((cum.shape[0],1)),cum],1),1)
    else:
        SIZE_MM=sizes_raw
        def link_frac(): return rbd._bed_surf__gsd_link

    gsd0_link=link_frac()               # “before” fractions

    rain=np.loadtxt(cwd/rainfall_file)
    rtime,rint,ridx,curr_r,rval=rainfall_data(rain,g)

    link_lab=_label_links(g,_bands(z_dem))
    node_lab=_label_nodes(g,_bands(z_dem))

    z0=g.at_node["topographic__elevation"].copy()

    ts={s:[] for s in SECTS}; lag={s:None for s in SECTS}; tlag={s:None for s in SECTS}
    mass={s:None for s in SECTS}
    tau_stats={s:OnlineStats(0.01) for s in SECTS}
    qb_stats ={s:OnlineStats(0.01) for s in SECTS}

    rbd._surface_water__velocity_prev_time_link=np.zeros(g.number_of_links)

    t=0.0
    while t<sim_max_t:
        if t==0.0: of.dt=max_dt; t+=of.dt; continue
        ridx,curr_r,rval,_=update_rainfall_intensity(
            t,rtime,ridx,rint,g,curr_r,rval,False)
        of._rainfall_intensity=curr_r
        of.overland_flow(dt=max_dt)

        if g.at_node["surface_water__depth"].max()>=DEPTH_MIN:
            g["link"]["surface_water__velocity"]=(
                g["link"]["surface_water__discharge"]/
                g["link"]["surface_water__depth"])
            rbd._surface_water__velocity_prev_time_link= \
                g["link"]["surface_water__velocity"]
            rbd._grid._dt=of.dt; rbd.run_one_step()

            qb_link =np.abs(rbd._sed_transp__bedload_rate_link)
            tau_link=np.abs(rbd._surface_water__shear_stress_link)
            tr      =link_frac()*qb_link[:,None]

            n_now=tr.shape[1]; sizes_now=SIZE_MM[:n_now]
            if mass[SECTS[0]] is None:
                mass={s:np.zeros(n_now) for s in SECTS}
            for s in SECTS: ts[s].append(0)

            for lid in range(g.number_of_links):
                sec=link_lab[lid]
                if sec=="none": continue
                tau_stats[sec].add(tau_link[lid]); qb_stats[sec].add(qb_link[lid])
                inc=tr[lid]*of.dt
                if inc.sum()<=EPS_MASS or curr_r.max()==0: continue
                idx=np.nonzero(inc>EPS_MASS)[0]
                if idx.size==0: continue
                mz=sizes_now[idx].max()
                ts[sec][-1]=mz; mass[sec]+=inc
                g.at_link["max_size_moved_mm"][lid]=max(g.at_link["max_size_moved_mm"][lid],mz)
                if lag[sec] is None or mz>lag[sec]:
                    lag[sec]=mz; tlag[sec]=t
        else:
            for s in SECTS: ts[s].append(0)

        of.dt=of.dt or max_dt*1e-3
        t,*_=update_time(t,of,sim_max_t,False,False,0)

    # rasters & PNGs ------------------------------------------------------
    dz=g.at_node["topographic__elevation"]-z0
    g.add_field("node","elev_change",dz,overwrite=True)
    write_esri_ascii(out/"elev_change.asc",g,names="elev_change",clobber=True)
    if plots: _safe_imsave(dz.reshape(g.shape),out/"elev_change.png")

    def _Dx_map(before,after,pct,name):
        D0=np.array([_Dx(b,SIZE_MM,pct) for b in before])
        D1=np.array([_Dx(a,SIZE_MM,pct) for a in after])
        d=D1-D0
        node=np.zeros(g.number_of_nodes); c=np.zeros_like(node)
        for lid,val in enumerate(d):
            n1,n2=g.nodes_at_link[lid]
            node[n1]+=val; node[n2]+=val; c[n1]+=1; c[n2]+=1
        c[c==0]=1; node/=c
        fn=f"{name}_change_mm_node"
        g.add_field("node",fn,node,overwrite=True)
        write_esri_ascii(out/f"{name}_change.asc",g,names=fn,clobber=True)
        if plots: _safe_imsave(node.reshape(g.shape),out/f"{name}_change.png")

    after=link_frac()
    _Dx_map(gsd0_link,after,0.10,"D10")
    _Dx_map(gsd0_link,after,0.50,"D50")
    _Dx_map(gsd0_link,after,0.90,"D90")
    _Dx_map(gsd0_link,after,1.00,"Dmax")

    for sec in SECTS:
        lids=np.where(link_lab==sec)[0]; nids=np.where(node_lab==sec)[0]
        with open(out/f"shear_stats_{sec}.txt","w",encoding="utf-8") as f:
            f.write("# |tau| Pa\n")
            for k,v in tau_stats[sec].summary().items(): f.write(f"{k} {v:.5e}\n")
        with open(out/f"bedload_stats_{sec}.txt","w",encoding="utf-8") as f:
            f.write("# |qb| kg m^-1 s^-1\n")
            for k,v in qb_stats[sec].summary().items(): f.write(f"{k} {v:.5e}\n")
        with open(out/f"elev_change_stats_{sec}.txt","w",encoding="utf-8") as f:
            arr=dz[nids]; p=np.percentile(arr,[5,50,95])
            f.write("# dz m\n")
            f.write(f"max {arr.max():.5e}\nmin {arr.min():.5e}\nmean {arr.mean():.5e}\n")
            f.write(f"median {p[1]:.5e}\np5 {p[0]:.5e}\np95 {p[2]:.5e}\n")
        gs0,gs1=gsd0_link[lids].mean(axis=0),after[lids].mean(axis=0)
        d10=_Dx(gs1,SIZE_MM,0.10)-_Dx(gs0,SIZE_MM,0.10)
        d50=_Dx(gs1,SIZE_MM,0.50)-_Dx(gs0,SIZE_MM,0.50)
        d90=_Dx(gs1,SIZE_MM,0.90)-_Dx(gs0,SIZE_MM,0.90)
        trim=min(gs0.size,SIZE_MM.size)
        idx0=np.where(gs0[:trim]>0)[0][-1] if gs0.sum() else 0
        idx1=np.where(gs1[:trim]>0)[0][-1] if gs1.sum() else 0
        dmax=SIZE_MM[idx1]-SIZE_MM[idx0]
        with open(out/f"surf_gsd_change_{sec}.txt","w",encoding="utf-8") as f:
            f.write("# surf GSD change mm\n")
            f.write(f"dD10 {d10:.2f}\ndD50 {d50:.2f}\ndD90 {d90:.2f}\ndDmax {dmax:.2f}\n")

    if plots:
        plt.figure()
        for sec in SECTS:
            plt.plot(np.arange(len(ts[sec]))*max_dt,ts[sec],label=sec)
        plt.xlabel("time [s]"); plt.ylabel("largest grain moved [mm]")
        plt.legend(); plt.tight_layout()
        plt.savefig(out/"largest_size_timeseries.png",dpi=300); plt.close()

    node_max=np.zeros(g.number_of_nodes)
    for lid,sz in enumerate(g.at_link["max_size_moved_mm"]):
        n1,n2=g.nodes_at_link[lid]
        node_max[n1]=max(node_max[n1],sz)
        node_max[n2]=max(node_max[n2],sz)
    g.add_field("node","max_size_moved_mm_node",node_max,overwrite=True)
    write_esri_ascii(out/"max_size_moved.asc",g,
                     names="max_size_moved_mm_node",clobber=True)

    with open(out/"summary.txt","w",encoding="utf-8") as f:
        for sec in SECTS:
            v=mass[sec]
            if v is None or v.sum()<=EPS_MASS:
                f.write(f"{sec}: no transport\n"); continue
            cf=np.cumsum(v)/v.sum()
            D10=SIZE_MM[np.searchsorted(cf,0.10)]
            D50=SIZE_MM[np.searchsorted(cf,0.50)]
            D90=SIZE_MM[np.searchsorted(cf,0.90)]
            f.write(
              f"{sec}: largest={lag[sec]:.1f} mm at {tlag[sec]:.0f} s; "
              f"D10/D50/D90={D10:.1f}/{D50:.1f}/{D90:.1f} mm; "
              f"mass={v.sum():.2e} kg\n")
    return out

# CLI
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--rainfall",default="working_climate/15min_1000yr.txt")
    p.add_argument("--tag",default="")
    p.add_argument("--sim_max_t",type=float,default=22000)
    p.add_argument("--plots",action="store_true")
    a=p.parse_args()
    folder=run_simulation(rainfall_file=a.rainfall,tag=a.tag,
                          sim_max_t=a.sim_max_t,plots=a.plots)
    print("results →",folder)
