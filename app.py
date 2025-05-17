#!/usr/bin/env python
# =======================================================================
# app.py  –  Interactive Circuit Canvas  (2025‑04‑20  ·  v7 - Final Polish)
# =======================================================================
#
# Shortcuts (Full legend: L key)
#   S=Solve, M=Info Pages, C=Edit Vals, F=Set Freq, W=Del Wires
#   Drag Palette=Place, Click Pin-Pin=Wire, R-Click=Delete Part
# -----------------------------------------------------------------------
from __future__ import annotations
import sys, cmath, math, itertools, re
import pygame as pg
import numpy  as np
import scipy.linalg

# ───────── UI SCALING FACTORS ─────────
SIZE_SCALE = 1.25
FONT_SCALE = 1.5

# ───────── UI CONSTANTS (Scaled) ─────────
WIN_W, WIN_H, PAL_W, INFO_H = int(1220*SIZE_SCALE), int(820*SIZE_SCALE), int(220*SIZE_SCALE), int(120*SIZE_SCALE)
CANVAS_W, CANVAS_H          = WIN_W-PAL_W, WIN_H-INFO_H
FPS, GRID, PIN_R            = 60, int(20*SIZE_SCALE), int(4*SIZE_SCALE)
GHOST_A                     = 120

FONT_PRIMARY_NAME = "DejaVu Sans Mono, Segoe UI, Arial, sans-serif"
FONT_MONO_NAME = "DejaVu Sans Mono, Consolas, Courier New, monospace"
FONT_SYMBOL_NAME = "DejaVu Sans Mono, Segoe UI Symbol, Arial Unicode MS, sans-serif"

# ───────── PART TYPES ───────────
class CType: R, C, L, VDC, VAC, GND = range(6)
LBL = {CType.R:"R", CType.C:"C", CType.L:"L", CType.VDC:"DC", CType.VAC:"AC", CType.GND:"GND"}
COL = {CType.R:(215,215,215), CType.C:(120,170,255), CType.L:(130,255,170),
       CType.VDC:(120,200,255), CType.VAC:(255,170,120), CType.GND:(150,150,255)}

# ───────── BASIC CLASSES ─────────
class Pin:
    def __init__(self, comp:'Comp', dx:int, dy:int, name:str=""):
        self.c, self.dx, self.dy = comp, int(dx*SIZE_SCALE), int(dy*SIZE_SCALE)
        self.net = -1; self.name = name
    @property
    def pos(self): return self.c.x+self.dx, self.c.y+self.dy

class Wire:
    def __init__(self, a:Pin, b:Pin): self.a=a; self.b=b

class Comp:
    comp_counters = {} # Initialized in App.__init__
    BW,BH = int(80*SIZE_SCALE), int(40*SIZE_SCALE)
    def __init__(self, ct:CType, x:int, y:int):
        self.ct,self.x,self.y,self.sel = ct,x,y,False
        self.val  = {CType.R:1e3, CType.C:1e-6, CType.L:1e-3, CType.VDC:5.0, CType.VAC:5.0}.get(ct,0.0)
        self.phase= 0.0
        
        Comp.comp_counters[ct] = Comp.comp_counters.get(ct, 0) + 1
        self.name = f"{LBL[ct]}{Comp.comp_counters[ct]}"

        pin_offset = int(30*SIZE_SCALE)
        if ct==CType.GND: self.pins = [Pin(self,0,0, name="ref")]
        else: self.pins = [Pin(self,-pin_offset,0, name="1"), Pin(self,pin_offset,0, name="2")]

    def rect(self):
        gnd_w, gnd_h = int(10*SIZE_SCALE), int(12*SIZE_SCALE)
        return pg.Rect(self.x-gnd_w, self.y-gnd_h, gnd_w*2, gnd_h*2) if self.ct==CType.GND \
               else pg.Rect(self.x-Comp.BW//2,self.y-Comp.BH//2,Comp.BW,Comp.BH)

    def label(self):
        if self.ct == CType.R: unit = "\u03A9"
        elif self.ct == CType.C: unit = "F"
        elif self.ct == CType.L: unit = "H"
        else: unit = ""
        if unit: return f"{value_to_str(self.val)}{unit}"
        if self.ct==CType.VDC: return f"{self.val:.3g}V"
        if self.ct==CType.VAC: return f"{self.val:.3g}V\u2220{self.phase:.0f}\u00B0"
        return ""
    def __repr__(self): return self.name

# ───────── VALUE PARSING / FORMAT ─────────
_SI = {"p":1e-12,"n":1e-9,"u":1e-6,"µ":1e-6, "m":1e-3,"":1,"k":1e3,"meg":1e6,"M":1e6,"g":1e9,"G":1e9}
def parse_value(txt:str) -> float|None:
    txt=txt.strip().lower()
    txt = re.sub(r"\s*(v|a|hz|f|h|ohm|\u03A9)\s*$", "", txt, flags=re.IGNORECASE)
    if "@" in txt: txt=txt.split("@")[0].strip()
    if "∠" in txt or "\u2220" in txt : txt=txt.split("∠")[0].split("\u2220")[0].strip()
    m=re.match(r"([-+]?[0-9.]+)\s*(p|n|u|µ|m|k|meg|M|g|G)?",txt)
    if not m: return None
    try: return float(m.group(1))*_SI.get(m.group(2) or "",1)
    except ValueError: return None

def value_to_str(v:float) -> str:
    if v == 0: return "0"
    prefixes = [("G",1e9),("M",1e6),("k",1e3),("",1),("m",1e-3),("u",1e-6),("n",1e-9),("p",1e-12)]
    for suf, fac in prefixes:
        if abs(v) >= fac:
            sv = v / fac
            fmt = "%.0f" if abs(sv) >= 100 else "%.1f" if abs(sv) >= 10 else "%.2f" if abs(sv) >= 1 else "%.3f"
            str_v = (fmt % sv).rstrip('0').rstrip('.') if '.' in (fmt % sv) else (fmt % sv)
            return f"{str_v}{suf}"
    return f"{v:.2e}"

# ───────── CIRCUIT / SOLVER ─────────
class Circuit:
    def __init__(self): # ... (same as v6)
        self.comps, self.wires = [], []
        self.lu_and_piv = None 
        self.invalidate()

    def invalidate(self): # ... (same as v6)
        self.solution=None; self.Y=self.b=None; self.lu_and_piv = None
        self.node_cnt=0; self.nmap={}; self.vsrc_idx={}; self.idx_to_nodename = {}
        self.net_id_to_comp_pins = {}

    def add_comp(self,c): self.comps.append(c); self.invalidate()
    def add_wire(self,a:Pin,b:Pin): # ... (same as v6)
        if a.c is b.c: return
        if any((w.a is a and w.b is b) or (w.a is b and w.b is a) for w in self.wires): return
        self.wires.append(Wire(a,b)); self.invalidate()
    def delete_comp(self,c_to_delete): # ... (same as v6)
        if c_to_delete in self.comps: self.comps.remove(c_to_delete)
        self.wires=[w for w in self.wires if w.a.c is not c_to_delete and w.b.c is not c_to_delete]
        self.invalidate()
    def delete_wire(self,w_to_delete):  # ... (same as v6)
        if w_to_delete in self.wires: self.wires.remove(w_to_delete)
        self.invalidate()

    def _nets(self)->tuple[bool, str | None]: # ... (same as v6)
        self.net_id_to_comp_pins.clear()
        for p in (p_ for c_ in self.comps for p_ in c_.pins): p.net=-1
        def bfs(start_pin:Pin, net_id:int):
            if start_pin.net != -1: return
            queue=[start_pin]; start_pin.net=net_id
            if net_id not in self.net_id_to_comp_pins: self.net_id_to_comp_pins[net_id] = []
            self.net_id_to_comp_pins[net_id].append(f"{start_pin.c.name}.{start_pin.name}")
            head = 0
            while head < len(queue):
                current_pin = queue[head]; head += 1
                for w in self.wires:
                    neighbor = None
                    if w.a is current_pin and w.b.net == -1: neighbor = w.b
                    elif w.b is current_pin and w.a.net == -1: neighbor = w.a
                    if neighbor:
                        neighbor.net = net_id
                        self.net_id_to_comp_pins[net_id].append(f"{neighbor.c.name}.{neighbor.name}")
                        queue.append(neighbor)
        nid=0 
        gnd_comps_found = [c for c in self.comps if c.ct == CType.GND]
        if gnd_comps_found: # Prioritize BFS from a GND if one exists
             if gnd_comps_found[0].pins[0].net == -1: bfs(gnd_comps_found[0].pins[0], nid); nid+=1
        
        for p in (p_ for c_ in self.comps for p_ in c_.pins):
            if p.net==-1: bfs(p,nid); nid+=1
        
        gnd_comps=[c for c in self.comps if c.ct==CType.GND]
        if not gnd_comps: return False,"No GND component. Circuit cannot be solved."
        if len(gnd_comps)>1: return False,"Multiple GNDs. Only one allowed."

        ground_net_id = gnd_comps[0].pins[0].net 
        if ground_net_id != 0: 
            for p in (p_ for c_ in self.comps for p_ in c_.pins):
                if p.net == ground_net_id: p.net = 0
                elif p.net == 0: p.net = ground_net_id
            if 0 in self.net_id_to_comp_pins and ground_net_id in self.net_id_to_comp_pins:
                 self.net_id_to_comp_pins[0], self.net_id_to_comp_pins[ground_net_id] = \
                 self.net_id_to_comp_pins[ground_net_id], self.net_id_to_comp_pins[0]
            elif ground_net_id in self.net_id_to_comp_pins : 
                 self.net_id_to_comp_pins[0] = self.net_id_to_comp_pins.pop(ground_net_id)
        return True, None

    def solve(self, freq:float): # Updated C, L DC behavior
        ok, msg = self._nets()
        if not ok: self.invalidate(); return False, msg

        omega = 2 * math.pi * freq
        self.nmap.clear(); self.idx_to_nodename.clear(); k=0
        
        unique_net_ids = sorted(list(set(p.net for c_ in self.comps for p in c_.pins if p.net > 0)))
        for net_id_val in unique_net_ids:
            self.nmap[net_id_val] = k
            connected_comps_str = ", ".join(list(set(cp_str.split('.')[0] for cp_str in self.net_id_to_comp_pins.get(net_id_val, []))))
            self.idx_to_nodename[k] = f"V(N{net_id_val}[{connected_comps_str[:15].strip(', ')}])" if connected_comps_str else f"V(N{net_id_val})"
            k += 1
        self.node_cnt = k

        v_sources = [c for c in self.comps if c.ct in (CType.VDC, CType.VAC)]
        m = len(v_sources)
        dim = self.node_cnt + m

        if dim == 0: # ... (same as v6)
            self.Y = np.array([], dtype=complex).reshape(0,0); self.b = np.array([], dtype=complex)
            self.solution = np.array([]); self.lu_and_piv = None
            return True, "Solved (empty circuit)"

        Y = np.zeros((dim, dim), dtype=complex)
        b = np.zeros(dim, dtype=complex)
        def idx(net_id_val): return -1 if net_id_val == 0 else self.nmap.get(net_id_val, -2)

        for c in self.comps:
            admittance = 0j # Use admittance for clarity
            if c.ct == CType.R: admittance = 1 / max(c.val, 1e-12)
            elif c.ct == CType.C:
                if omega == 0: admittance = 0j # Ideal capacitor is open at DC
                elif c.val == 0: admittance = 0j # 0F capacitor is open
                else: admittance = 1j * omega * c.val
            elif c.ct == CType.L:
                if omega == 0: admittance = 1e12 # Ideal inductor is short at DC (large admittance)
                elif c.val == 0: admittance = 1e12 # 0H inductor is short
                else: admittance = 1 / (1j * omega * c.val) # Standard AC

            if c.ct in (CType.R, CType.C, CType.L):
                if len(c.pins) != 2: continue
                n_a, n_b = c.pins[0].net, c.pins[1].net 
                i, j = idx(n_a), idx(n_b)
                if i >= 0: Y[i, i] += admittance
                if j >= 0: Y[j, j] += admittance
                if i >= 0 and j >= 0: Y[i, j] -= admittance; Y[j, i] -= admittance
        
        self.vsrc_idx.clear()
        for s_idx, c_vs in enumerate(v_sources): # ... (same as v6, but ensure idx_to_nodename uses Comp.name for source)
            if len(c_vs.pins) != 2: continue
            mna_row_idx = self.node_cnt + s_idx
            self.vsrc_idx[c_vs] = s_idx
            self.idx_to_nodename[mna_row_idx] = f"I({c_vs.name}:{c_vs.pins[0].name}\u2192{c_vs.pins[1].name})"

            node_pos_idx, node_neg_idx = idx(c_vs.pins[1].net), idx(c_vs.pins[0].net) 
            if node_pos_idx >= 0: Y[node_pos_idx, mna_row_idx] += 1; Y[mna_row_idx, node_pos_idx] += 1
            if node_neg_idx >= 0: Y[node_neg_idx, mna_row_idx] -= 1; Y[mna_row_idx, node_neg_idx] -= 1
            phasor = cmath.rect(c_vs.val, math.radians(c_vs.phase if c_vs.ct == CType.VAC else 0))
            b[mna_row_idx] = phasor
        
        self.Y, self.b = Y, b
        try: # ... (same as v6)
            if np.all(np.abs(Y) < 1e-18) and np.all(np.abs(b) < 1e-18) and dim > 0 :
                 self.solution = np.zeros_like(b) 
                 self.lu_and_piv = scipy.linalg.lu_factor(np.eye(dim) if dim > 0 else np.array([]).reshape(0,0), check_finite=False) if dim >0 else None # Factor identity for placeholder
            else:
                 self.lu_and_piv = scipy.linalg.lu_factor(Y, check_finite=False)
                 self.solution = scipy.linalg.lu_solve(self.lu_and_piv, b, check_finite=False)
            return True, "Solved (LU)"
        except (np.linalg.LinAlgError, ValueError, scipy.linalg.LinAlgError) as e:
            self.solution = None; self.lu_and_piv = None
            if "singular matrix" in str(e).lower() or isinstance(e, (np.linalg.LinAlgError, scipy.linalg.LinAlgError)):
                 return False, f"Singular matrix (check connections, GND)"
            return False, f"Solver error: {type(e).__name__}"

# ───────── APP ────────────
class App:
    def __init__(self): # ... (same as v6 for font loading) ...
        pg.init() 
        try: self.font = pg.font.SysFont(FONT_PRIMARY_NAME, int(15 * FONT_SCALE))
        except: self.font = pg.font.Font(None, int(15 * FONT_SCALE)) 
        try: self.font_small = pg.font.SysFont(FONT_MONO_NAME, int(12 * FONT_SCALE))
        except: self.font_small = pg.font.Font(None, int(12 * FONT_SCALE))
        try: self.font_tiny = pg.font.SysFont(FONT_MONO_NAME, int(10 * FONT_SCALE)) 
        except: self.font_tiny = pg.font.Font(None, int(10 * FONT_SCALE))
        try: self.font_symbol = pg.font.SysFont(FONT_SYMBOL_NAME, int(15 * FONT_SCALE)) 
        except: self.font_symbol = self.font 

        self.scr=pg.display.set_mode((WIN_W,WIN_H)); pg.display.set_caption("Circuit Canvas v7 (Final Polish)")
        self.clock=pg.time.Clock(); self.circ=Circuit()
        self.drag_type:CType|None = None; self.drag_inst:Comp|None = None
        self.offx=self.offy=0
        self.wire_start:Pin|None = None
        self.mode_edit = False; self.mode_wire_del = False
        self.show_page = 0; self.show_help = False
        self.status = "Drag parts to start. (L = Legend)"
        self.freq = 0.0
        self.input_buffer = ""; self.prompt = ""; self.active_edit:Comp|None = None

        self.pal_rects={}; pal_item_x = CANVAS_W + int(24*SIZE_SCALE)
        pal_item_y_start, pal_item_y_step = int(28*SIZE_SCALE), int(58*SIZE_SCALE)
        pal_item_w, pal_item_h = PAL_W - int(48*SIZE_SCALE), int(44*SIZE_SCALE)
        self.pal_tooltips = {
            CType.R: "Resistor (R)", CType.C: "Capacitor (C)", CType.L: "Inductor (L)",
            CType.VDC: "DC Voltage Source", CType.VAC: "AC Voltage Source", CType.GND: "Ground (0V)"
        }
        for i, t in enumerate((CType.R, CType.C, CType.L, CType.VDC, CType.VAC, CType.GND)):
            self.pal_rects[t] = pg.Rect(pal_item_x, pal_item_y_start + i*pal_item_y_step, pal_item_w, pal_item_h)
        Comp.comp_counters = {ct_val: 0 for ct_val in CType.__dict__.values() if isinstance(ct_val, int)}


    # ... (get_comp_at, get_pin_at, get_wire_at are fine)
    def get_comp_at(self,x,y): # unchanged
        for c in reversed(self.circ.comps):
            if c.rect().collidepoint(x,y): return c
        return None
    def get_pin_at(self,x,y): # unchanged
        click_r_sq = (PIN_R + int(3*SIZE_SCALE))**2
        for c in self.circ.comps:
            for p in c.pins:
                if (p.pos[0]-x)**2+(p.pos[1]-y)**2 <= click_r_sq: return p
        return None
    def get_wire_at(self,x,y): # unchanged
        tol_sq = int(6*SIZE_SCALE)**2
        for w in self.circ.wires:
            ax,ay=w.a.pos; bx,by=w.b.pos; dx,dy=bx-ax,by-ay
            len_sq = dx*dx + dy*dy
            if len_sq == 0: t = 0
            else: t = max(0, min(1, ((x-ax)*dx + (y-ay)*dy)/len_sq))
            cx, cy = ax + t*dx, ay + t*dy
            if (x-cx)**2+(y-cy)**2 < tol_sq: return w
        return None

    def run(self): # ... (same)
        while True:
            self.clock.tick(FPS)
            self.handle_events()
            self.draw()

    def handle_events(self): # ... (same as v6)
        mx, my = pg.mouse.get_pos()
        ctrl_pressed = pg.key.get_mods() & pg.KMOD_CTRL
        shift_pressed = pg.key.get_mods() & pg.KMOD_SHIFT

        for e in pg.event.get():
            if e.type == pg.QUIT: pg.quit(); sys.exit()

            if self.prompt and e.type == pg.KEYDOWN:
                if e.key == pg.K_RETURN or e.key == pg.K_KP_ENTER: self.process_prompt_input()
                elif e.key == pg.K_BACKSPACE: self.input_buffer = self.input_buffer[:-1]
                elif e.key == pg.K_ESCAPE: self.reset_action(); self.status = "Input cancelled."
                elif e.unicode.isprintable(): self.input_buffer += e.unicode
                continue

            if e.type == pg.KEYDOWN:
                if e.key == pg.K_ESCAPE:
                    if self.show_help: self.show_help = False; self.status = "Legend closed."
                    else: self.reset_action(); self.status = "Action cancelled."
                elif e.key in (pg.K_DELETE, pg.K_BACKSPACE):
                    if self.delete_selection(): self.status = "Deleted selected."
                elif e.key == pg.K_s:
                    ok, msg = self.circ.solve(self.freq)
                    self.status = msg + (" (M=View)" if ok else "")
                    if ok: self.show_page = 0 
                elif e.key == pg.K_m:
                    if self.circ.Y is not None and self.circ.Y.shape[0] > 0:
                        self.show_page = (self.show_page + 1) % 4 
                        pages = ["Solution Overlay", "Heatmap", "Matrix Text", "Circuit Info"] 
                        self.status = f"Info: {pages[self.show_page]}"
                    else: self.status = "Solve first (S) or add components."
                elif e.key == pg.K_e: 
                    if self.circ.Y is not None and self.circ.Y.shape[0] > 0:
                        self.show_page = 3; self.status = "Info: Circuit Info" 
                    else: self.status = "Solve first (S) or add components."
                elif e.key == pg.K_c:
                    self.mode_edit = not self.mode_edit
                    if self.mode_edit: self.mode_wire_del=False; self.reset_action(keep_edit=True)
                    self.status = "Component Edit "+("ON (Click part)" if self.mode_edit else "OFF")
                elif e.key == pg.K_f:
                    self.reset_action(keep_freq=True)
                    self.prompt = f"Frequency (now={value_to_str(self.freq)}Hz): "; self.input_buffer = ""
                    self.status = self.prompt
                elif e.key == pg.K_w:
                    self.mode_wire_del = not self.mode_wire_del
                    if self.mode_wire_del: self.mode_edit=False; self.reset_action(keep_wire_del=True)
                    self.status = "Wire Delete "+("ON (Click wire)" if self.mode_wire_del else "OFF")
                elif e.key == pg.K_a and ctrl_pressed:
                    for c_ in self.circ.comps: c_.sel = True
                    self.status = "Selected all."
                elif e.key in (pg.K_l, pg.K_h, pg.K_QUESTION):
                    self.show_help = not self.show_help
                    self.status = "Legend " + ("shown" if self.show_help else "hidden") + " (L/Esc)"


            if e.type == pg.MOUSEBUTTONDOWN: 
                if self.show_help and not self.prompt: 
                    self.show_help = False; self.status = "Legend closed."; continue
                if e.button == 1: 
                    if mx > CANVAS_W: 
                        self.reset_action(keep_freq=True)
                        for t, r_pal in self.pal_rects.items():
                            if r_pal.collidepoint(mx, my):
                                self.drag_type = t; self.status = f"Dragging {LBL[t]}..."
                                break
                    else: 
                        if self.mode_wire_del:
                            w = self.get_wire_at(mx, my)
                            if w: self.circ.delete_wire(w); self.status="Wire deleted."; self.circ.solve(self.freq)
                        elif self.mode_edit:
                            c_edit = self.get_comp_at(mx, my)
                            if c_edit and c_edit.ct!=CType.GND:
                                self.active_edit = c_edit
                                self.prompt = f"New val for {c_edit!r} ({c_edit.label()}): "; self.input_buffer=""
                                self.status = self.prompt
                        else: 
                            p_clk = self.get_pin_at(mx, my)
                            c_clk = self.get_comp_at(mx, my) if not p_clk else None
                            if p_clk:
                                self.drag_type = None 
                                if self.wire_start and p_clk is not self.wire_start and p_clk.c is not self.wire_start.c:
                                    self.circ.add_wire(self.wire_start, p_clk); self.status = "Wire added."
                                    self.wire_start = None; self.circ.solve(self.freq)
                                elif self.wire_start is p_clk: self.wire_start=None; self.status="Wiring cancelled."
                                else: self.wire_start = p_clk; self.status = "Click another pin..."
                            elif c_clk:
                                self.drag_type=None; self.wire_start=None
                                if not ctrl_pressed and not c_clk.sel: 
                                    for c_oth in self.circ.comps: c_oth.sel=False
                                c_clk.sel = not c_clk.sel if ctrl_pressed else True 
                                self.status = f"{c_clk!r} selected."
                                self.drag_inst = c_clk 
                                self.offx = c_clk.x - mx; self.offy = c_clk.y - my
                            else: 
                                self.drag_type=None; self.wire_start=None
                                if not ctrl_pressed: 
                                    for c_oth in self.circ.comps: c_oth.sel=False
                                self.status = ""
                elif e.button == 3: 
                    if not self.prompt:
                        c_del = self.get_comp_at(mx, my)
                        if c_del: self.circ.delete_comp(c_del); self.status=f"Deleted {c_del!r}."; self.circ.solve(self.freq)
            
            if e.type == pg.MOUSEBUTTONUP: 
                if e.button == 1:
                    if self.drag_type is not None and mx < CANVAS_W:
                        self.place(mx,my) 
                    self.drag_type = None; self.drag_inst = None 
            
            if e.type == pg.MOUSEMOTION: 
                if self.drag_inst is not None and e.buttons[0]: 
                    mdx = mx + self.offx - self.drag_inst.x
                    mdy = my + self.offy - self.drag_inst.y
                    moved_any = False
                    for c_mv in self.circ.comps:
                        if c_mv.sel:
                            new_x, new_y = c_mv.x + mdx, c_mv.y + mdy
                            hbw, hbh = (Comp.BW//2, Comp.BH//2) if c_mv.ct != CType.GND else (int(10*SIZE_SCALE), int(12*SIZE_SCALE))
                            c_mv.x = max(hbw, min(CANVAS_W - hbw, new_x))
                            gnd_y_visual_offset = 0 
                            c_mv.y = max(hbh + gnd_y_visual_offset, min(CANVAS_H - hbh - gnd_y_visual_offset, new_y))
                            moved_any = True
                    if moved_any: self.circ.invalidate() 
                
                if mx > CANVAS_W and self.drag_type is None and not self.prompt:
                    current_tooltip = ""
                    for t, r_pal in self.pal_rects.items():
                        if r_pal.collidepoint(mx, my):
                            current_tooltip = self.pal_tooltips.get(t, "")
                            break
                    if current_tooltip != self.status : self.status = current_tooltip 
            
            if e.type == pg.MOUSEWHEEL: 
                 if not self.prompt:
                    c_adj = self.get_comp_at(mx, my)
                    if c_adj and c_adj.ct in (CType.R,CType.C,CType.L,CType.VDC,CType.VAC):
                        factor = 10.0 if e.y > 0 else 0.1
                        if shift_pressed: factor = (10**0.5) if e.y > 0 else (1/(10**0.5)) 
                        c_adj.val *= factor
                        if c_adj.ct in (CType.R,CType.C,CType.L): c_adj.val = max(c_adj.val, 1e-15)
                        self.status = f"Set {c_adj!r} to {c_adj.label()}"; self.circ.invalidate(); self.circ.solve(self.freq)

    def process_prompt_input(self): # ... (same as v6)
        input_str = self.input_buffer.strip()
        if not input_str: self.reset_action(); self.status = "Edit cancelled."; return

        if self.prompt.startswith("New val"):
            c = self.active_edit
            if c is None: self.reset_action(); return
            new_val = parse_value(input_str)
            new_phase = None
            if c.ct == CType.VAC:
                 match_phase = re.search(r"(?:∠|\u2220)\s*([-+]?\s*\d+\.?\d*)", self.input_buffer) 
                 if match_phase:
                     try: new_phase = float(match_phase.group(1).replace(" ",""))
                     except ValueError: new_phase = None
            
            changed = False
            if new_val is not None:
                if c.ct in (CType.R,CType.C,CType.L):
                    if new_val > 1e-15: c.val = new_val; changed = True
                else: c.val = new_val; changed = True 
            if new_phase is not None and c.ct == CType.VAC: c.phase = new_phase; changed = True
            
            if changed:
                self.status=f"Updated {c!r} to {c.label()}"; self.circ.invalidate(); self.circ.solve(self.freq)
            else: self.status = f"Invalid input: '{self.input_buffer}'"
        elif self.prompt.startswith("Frequency"):
            new_freq = parse_value(input_str)
            if new_freq is not None and new_freq >= 0:
                self.freq = new_freq; self.status = f"Freq set to {value_to_str(self.freq)}Hz."
                self.circ.invalidate(); self.circ.solve(self.freq) 
            else: self.status = f"Invalid freq: '{self.input_buffer}'"
        self.prompt = ""; self.input_buffer = ""; self.active_edit = None

    def reset_action(self, keep_freq=False, keep_edit=False, keep_wire_del=False): # ... (same as v6)
        self.drag_type = None; self.drag_inst = None; self.wire_start = None
        if not keep_edit: self.mode_edit = False
        if not keep_wire_del: self.mode_wire_del = False
        self.prompt = ""; self.input_buffer = ""; self.active_edit = None
        if not (keep_edit or keep_wire_del or self.show_help or self.show_page > 0): 
            for c_ in self.circ.comps: c_.sel = False # Renamed loop var
        
        mx, my = pg.mouse.get_pos()
        is_tooltip_active = mx > CANVAS_W and self.status.startswith(tuple(self.pal_tooltips.values()))
        is_info_page_status = self.status.startswith("Info:")
        
        if not (keep_freq or is_tooltip_active or is_info_page_status or self.show_help):
            self.status = "Cancelled."

    def delete_selection(self): # ... (same as v6)
        sel = [c for c in self.circ.comps if c.sel]
        if not sel: return False
        for c_del in sel: self.circ.delete_comp(c_del) # Renamed loop var
        self.circ.solve(self.freq) 
        return True

    def place(self, x: int, y: int): # ... (same as v6)
        if self.drag_type is None: self.status = "Error: No type selected."; return
        ctype = self.drag_type
        gx, gy = round(x/GRID)*GRID, round(y/GRID)*GRID
        hbw,hbh = (Comp.BW//2, Comp.BH//2) if ctype!=CType.GND else (int(10*SIZE_SCALE),int(12*SIZE_SCALE))
        gx = max(hbw, min(CANVAS_W-hbw, gx))
        gy = max(hbh, min(CANVAS_H-hbh, gy)) 
        
        new_c = Comp(ctype, gx, gy)
        self.circ.add_comp(new_c)
        self.status = f'Placed {new_c!r}'
        
        SNAP_SQ = (PIN_R * 3.5)**2; snapped_wire = False
        other_pins = [p for c_ in self.circ.comps if c_ is not new_c for p in c_.pins]
        for p1 in new_c.pins:
            for p2 in other_pins:
                if (p1.pos[0]-p2.pos[0])**2 + (p1.pos[1]-p2.pos[1])**2 <= SNAP_SQ:
                    self.circ.add_wire(p1, p2); snapped_wire = True
        if snapped_wire: self.status += ' (auto-wired)'
        self.circ.solve(self.freq)
    
    def draw(self): # Added draw_mini_legend
        self.scr.fill((40, 45, 50))
        self.draw_grid()
        self.draw_wires()
        self.draw_comps() 
        self.draw_palette()
        self.draw_info_bar() 
        if self.show_help: self.draw_legend()
        else: self.draw_mini_legend() # Draw mini-legend if full is not shown
        pg.display.flip()

    def draw_grid(self): # ... (same as v6)
        grid_c = (60,65,70)
        for x_ in range(0, CANVAS_W, GRID): pg.draw.line(self.scr, grid_c, (x_,0), (x_,CANVAS_H))
        for y_ in range(0, CANVAS_H, GRID): pg.draw.line(self.scr, grid_c, (0,y_), (CANVAS_W,y_))

    def draw_wires(self): # ... (same as v6)
        wire_c, active_c, del_c = (220,220,130), (255,255,0), (255,100,100)
        thick, active_thick = max(1,int(2*SIZE_SCALE)), max(1,int(1*SIZE_SCALE))
        mx,my = pg.mouse.get_pos()
        hover_w = self.get_wire_at(mx,my) if self.mode_wire_del else None
        for w in self.circ.wires:
            col = del_c if w is hover_w else wire_c # Renamed c to col
            pg.draw.line(self.scr, col, w.a.pos, w.b.pos, thick)
        if self.wire_start:
            pg.draw.line(self.scr, active_c, self.wire_start.pos, (mx,my), active_thick)

    def draw_comps(self): # ... (same as v6)
        for c in self.circ.comps: self.draw_single_comp(c) 
        if self.drag_type is not None and pg.mouse.get_pos()[0] < CANVAS_W:
            mx,my = pg.mouse.get_pos()
            gx,gy = round(mx/GRID)*GRID, round(my/GRID)*GRID
            hbw,hbh = (Comp.BW//2, Comp.BH//2) if self.drag_type!=CType.GND else (int(10*SIZE_SCALE),int(12*SIZE_SCALE))
            gx = max(hbw, min(CANVAS_W-hbw, gx))
            gy = max(hbh, min(CANVAS_H-hbh, gy))
            ghost_color_rgb = COL[self.drag_type]
            rect_w = Comp.BW if self.drag_type != CType.GND else hbw*2
            rect_h = Comp.BH if self.drag_type != CType.GND else hbh*2
            ghost_s = pg.Surface((rect_w, rect_h), pg.SRCALPHA)
            ghost_rect_on_surf = pg.Rect(0,0,rect_w, rect_h) 
            pg.draw.rect(ghost_s, (*ghost_color_rgb, GHOST_A), ghost_rect_on_surf, border_radius=int(6*SIZE_SCALE))
            label_s = self.font.render(LBL[self.drag_type], True, (255,255,255, GHOST_A+80))
            ghost_s.blit(label_s, label_s.get_rect(center=ghost_rect_on_surf.center))
            self.scr.blit(ghost_s, (gx - rect_w//2, gy - rect_h//2))

    def draw_single_comp(self, c: Comp): # ... (same as v6)
        r = c.rect(); base_col = COL[c.ct]; border_c=(60,60,60); text_c=(0,0,0)
        pin_c=(240,240,120); sel_c=(255,255,0); edit_c=(255,100,0)
        br = int(6*SIZE_SCALE); thick=max(1,int(SIZE_SCALE)); sel_inf = int(6*SIZE_SCALE)
        if c.ct == CType.GND: 
            x,y = c.x,c.y; stem_h=int(12*SIZE_SCALE); bar_y0=int(4*SIZE_SCALE)
            bar_ystep=int(4*SIZE_SCALE); bar_ws=[int(w*SIZE_SCALE) for w in (14,10,6)]
            gnd_thick = max(1, int(2*SIZE_SCALE))
            pg.draw.line(self.scr,base_col,(x,y-stem_h),(x,y),gnd_thick)
            for i,w_bar in enumerate(bar_ws): pg.draw.line(self.scr,base_col,(x-w_bar//2,y+bar_y0+i*bar_ystep),(x+w_bar//2,y+bar_y0+i*bar_ystep),gnd_thick)
        else:
            pg.draw.rect(self.scr,base_col,r,border_radius=br)
            pg.draw.rect(self.scr,border_c,r,thick,border_radius=br)
            name_s = self.font.render(c.name, True, text_c) 
            self.scr.blit(name_s, name_s.get_rect(center=(c.x, r.top + int(10*SIZE_SCALE))))
            val_s = self.font_symbol.render(c.label(), True, text_c) 
            self.scr.blit(val_s, val_s.get_rect(center=(c.x, r.bottom - int(10*SIZE_SCALE))))
        for p in c.pins: 
            pg.draw.circle(self.scr,pin_c,p.pos,PIN_R)
            pg.draw.circle(self.scr,border_c,p.pos,PIN_R,thick)
        if c.sel: pg.draw.rect(self.scr,sel_c,r.inflate(sel_inf,sel_inf),max(1,int(2*SIZE_SCALE)),border_radius=br+int(2*SIZE_SCALE))
        if self.mode_edit and c.ct!=CType.GND and r.collidepoint(pg.mouse.get_pos()):
             pg.draw.rect(self.scr,edit_c,r.inflate(sel_inf,sel_inf),max(1,int(2*SIZE_SCALE)),border_radius=br+int(2*SIZE_SCALE))
        if self.circ.solution is not None and self.show_page == 0 and len(self.circ.solution) > 0 :
            self.draw_solution_overlay(c) 

    def draw_solution_overlay(self, c: Comp): # ... (same as v6)
        if not self.circ.solution.size: return
        volt_col=(120,255,120); curr_col=(255,200,120); font_to_use=self.font_symbol
        node_voltages=self.circ.solution[:self.circ.node_cnt]; vsrc_currents=self.circ.solution[self.circ.node_cnt:]
        net_map=self.circ.nmap
        for p_idx, p in enumerate(c.pins):
            volt_phasor=0j
            if p.net > 0 and p.net in net_map:
                node_idx=net_map[p.net]
                if 0 <= node_idx < len(node_voltages): volt_phasor=node_voltages[node_idx]
            mag, ph_rad = cmath.polar(volt_phasor); ph_deg = math.degrees(ph_rad)
            v_str = "0V" if mag < 1e-9 else f"{value_to_str(mag)}V\u2220{ph_deg:.0f}\u00B0"
            v_surf = font_to_use.render(v_str, True, volt_col)
            v_off_x = PIN_R + int(3*SIZE_SCALE)
            text_anchor_x = p.pos[0] + v_off_x if p.dx >= 0 else p.pos[0] - v_off_x
            anchor_point_key = 'midleft' if p.dx >= 0 else 'midright'
            if p.dx == 0 and p.dy == 0: 
                v_rect = v_surf.get_rect(center=(p.pos[0], p.pos[1] - PIN_R - int(5*SIZE_SCALE)))
            else: v_rect = v_surf.get_rect(**{anchor_point_key: (text_anchor_x, p.pos[1])})
            self.scr.blit(v_surf, v_rect)
        current_phasor = None 
        if c.ct in (CType.R,CType.C,CType.L) and len(c.pins)==2:
            p1,p2=c.pins[0],c.pins[1]; v_at_p1,v_at_p2=0j,0j
            if p1.net>0 and p1.net in net_map and net_map[p1.net]<len(node_voltages): v_at_p1=node_voltages[net_map[p1.net]]
            if p2.net>0 and p2.net in net_map and net_map[p2.net]<len(node_voltages): v_at_p2=node_voltages[net_map[p2.net]]
            omega=2*math.pi*self.freq; imp=float('inf')
            if c.ct==CType.R: imp=max(c.val,1e-12)
            elif c.ct==CType.C: imp=1/(1j*omega*c.val) if omega!=0 and c.val!=0 else float('inf')
            elif c.ct==CType.L: imp=(1j*omega*c.val) if omega!=0 else 1e-12 # Use impedance for L
            if imp!=float('inf') and abs(imp)>1e-12: current_phasor=(v_at_p1-v_at_p2)/imp # I = (V1-V2)/Z
            else: current_phasor=0j 
        elif c.ct in (CType.VDC,CType.VAC):
            vs_idx=self.circ.vsrc_idx.get(c)
            if vs_idx is not None and 0<=vs_idx<len(vsrc_currents): current_phasor=vsrc_currents[vs_idx]
        if current_phasor is not None:
            mag,ph_rad=cmath.polar(current_phasor); ph_deg=math.degrees(ph_rad)
            curr_str="0A" if mag<1e-9 else f"{value_to_str(mag)}A\u2220{ph_deg:.0f}\u00B0"
            curr_surf=font_to_use.render(curr_str,True,curr_col)
            self.scr.blit(curr_surf,curr_surf.get_rect(midtop=(c.x,c.rect().bottom+int(3*SIZE_SCALE))))

    def draw_palette(self): # ... (same as v6)
        pg.draw.rect(self.scr, (45,50,55), (CANVAS_W,0,PAL_W,WIN_H))
        br = int(6*SIZE_SCALE); thick = max(1,int(SIZE_SCALE))
        mx,my = pg.mouse.get_pos()
        for t, r_pal in self.pal_rects.items():
            is_hover = r_pal.collidepoint(mx,my) and not self.drag_type
            is_active_drag = self.drag_type == t
            base_c = COL[t]
            draw_c = [min(255,x+30) for x in base_c] if is_hover or is_active_drag else base_c
            pg.draw.rect(self.scr, draw_c, r_pal, border_radius=br)
            pg.draw.rect(self.scr, (80,80,80), r_pal, thick, border_radius=br)
            lbl_s = self.font.render(LBL[t],True,(0,0,0))
            self.scr.blit(lbl_s, lbl_s.get_rect(center=r_pal.center))

    def draw_info_bar(self): # ... (same as v6)
        bar_r = pg.Rect(0,CANVAS_H,WIN_W,INFO_H)
        pg.draw.rect(self.scr, (28,30,32), bar_r)
        pg.draw.line(self.scr, (60,65,70), (0,CANVAS_H), (WIN_W,CANVAS_H), max(1,int(SIZE_SCALE)))
        disp_txt = self.status; txt_c = (220,220,255); cur = ""
        if self.prompt:
            disp_txt = self.prompt + self.input_buffer; txt_c = (255,255,180)
            if int(pg.time.get_ticks()/400)%2==0: cur="_"
        elif self.freq > 0 and self.circ.solution is not None and len(self.circ.solution)>0: disp_txt += f"    |    AC: {value_to_str(self.freq)}Hz"
        elif self.freq == 0 and self.circ.solution is not None and len(self.circ.solution)>0: disp_txt += f"    |    DC Analysis"
        status_s = self.font.render(disp_txt+cur, True, txt_c)
        self.scr.blit(status_s, (int(15*SIZE_SCALE), CANVAS_H + int(5*SIZE_SCALE))) 
        if self.show_page > 0 and self.circ.Y is not None and self.circ.Y.shape[0] > 0:
            page_names = {1:"Heatmap", 2:"Matrix Text", 3:"Circuit Info"} 
            page_title = f"Page {self.show_page}/{len(page_names)}: {page_names[self.show_page]}"
            title_s = self.font_small.render(page_title, True, (150,150,150))
            self.scr.blit(title_s, title_s.get_rect(topright=(WIN_W-int(10*SIZE_SCALE), CANVAS_H + int(25*SIZE_SCALE))))
            draw_func = {1:self.draw_heat_map, 2:self.draw_matrix_text, 3:self.draw_stamps_info}.get(self.show_page)
            if draw_func: draw_func()
    
    def draw_legend(self): # ... (same as v6, with minor content updates if needed)
        # Use the sections dict from v6, which was already improved.
        # The drawing logic for the full legend from v6 can be used here.
        # For brevity, I'll skip pasting the full legend drawing code again if it's unchanged from v6.
        # Assume the v6 legend drawing code is here.
        sections = { # From v6
            "Mouse": [("Drag Palette","Place part"), ("Pin→Pin","Wire"), ("Drag Part","Move"),
                      ("R-Click Part","Del"), ("Wheel","Val ×/÷10"), ("Sh+Wheel","Val ×/÷√10"),("Ctrl+Click","Multi-Sel")],
            "Keyboard": [("S","Solve"), ("C","Edit Vals"), ("F","Set Freq"),("W","Del Wires"),
                         ("M","Info Page"),("E","Circuit Info Page"), ("Del","Del Sel"),("Ctrl+A","Sel All"),
                         ("Esc","Cancel/Deselect"),("L/H/?","This Legend")],
            "Info Pages": [("0","Solution Overlay"),("1","Matrix Heatmap"),("2","Matrix Text"),("3","Circuit Info Details")]
        }
        legend_bg=(25,30,35,235); hdr_c=(180,220,255); txt_c=(220,220,220); key_c=(255,255,180)
        bord_c=(100,120,140,220); pad=int(20*SIZE_SCALE); txt_sp=int(5*SIZE_SCALE)
        sect_sp=int(15*SIZE_SCALE); br=int(12*SIZE_SCALE); bord_th=max(1,int(2*SIZE_SCALE))
        f_leg=self.font_small; f_hdr=self.font; f_title=pg.font.SysFont(FONT_PRIMARY_NAME,int(18*FONT_SCALE), bold=True)
        max_kw=0; max_dw=0; sect_hs=[]; line_h_leg=f_leg.get_linesize()+txt_sp; hdr_h_leg=f_hdr.get_linesize()+txt_sp
        for sn, itms in sections.items(): # Calculate dimensions
            sh=hdr_h_leg; current_max_dw_section = 0
            for k,d_orig in itms:
                max_kw=max(max_kw, f_leg.render(k+":",1,key_c).get_width())
                words=d_orig.split(' '); current_line_desc=""; desc_lines_count=1; temp_max_desc_w_for_item=0
                for word in words:
                    test_line_desc=current_line_desc+word+" "
                    if f_leg.render(test_line_desc,1,txt_c).get_width() > (WIN_W*0.35): # Wrap width
                        temp_max_desc_w_for_item=max(temp_max_desc_w_for_item, f_leg.render(current_line_desc,1,txt_c).get_width())
                        current_line_desc=word+" "; desc_lines_count+=1
                    else: current_line_desc=test_line_desc
                temp_max_desc_w_for_item=max(temp_max_desc_w_for_item, f_leg.render(current_line_desc.strip(),1,txt_c).get_width())
                current_max_dw_section=max(current_max_dw_section, temp_max_desc_w_for_item)
                sh+=line_h_leg*desc_lines_count
            max_dw=max(max_dw, current_max_dw_section); sect_hs.append(sh)
        total_th=sum(sect_hs)+sect_sp*(len(sect_hs)-1) if sect_hs else 0
        col_w=max_kw+max_dw+int(10*SIZE_SCALE)+int(10*SIZE_SCALE) 
        leg_w=min(WIN_W-pad*2, col_w+pad*2); title_h_leg=f_title.get_linesize()+txt_sp
        leg_h=total_th+title_h_leg+pad*2 
        leg_r=pg.Rect((WIN_W-leg_w)//2, max(int(10*SIZE_SCALE),(WIN_H-leg_h)//2), leg_w, min(WIN_H-int(20*SIZE_SCALE),leg_h))
        leg_s=pg.Surface(leg_r.size, pg.SRCALPHA)
        pg.draw.rect(leg_s,leg_bg,leg_s.get_rect(),border_radius=br); pg.draw.rect(leg_s,bord_c,leg_s.get_rect(),bord_th,border_radius=br)
        title_s=f_title.render("CONTROLS & REFERENCE",1,(255,255,255)); leg_s.blit(title_s,title_s.get_rect(midtop=(leg_w//2,pad)))
        cur_y=title_s.get_rect(midtop=(leg_w//2,pad)).bottom+pad//2; x_st=pad; itm_ind=int(10*SIZE_SCALE); k_d_gap=int(10*SIZE_SCALE)
        for idx,(sn,itms) in enumerate(sections.items()): # Draw content
            hdr_s=f_hdr.render(sn,1,hdr_c); leg_s.blit(hdr_s,(x_st,cur_y)); cur_y+=hdr_h_leg
            for k,d_orig in itms:
                key_s=f_leg.render(k+":",1,key_c); leg_s.blit(key_s,(x_st+itm_ind,cur_y))
                words=d_orig.split(' '); current_line_desc=""; line_start_x_desc=x_st+itm_ind+max_kw+k_d_gap
                temp_y_desc=cur_y; initial_cur_y_for_item=cur_y
                for word_idx,word in enumerate(words):
                    test_line_desc=current_line_desc+word+" "
                    if line_start_x_desc+f_leg.render(test_line_desc,1,txt_c).get_width() > leg_w-pad or word_idx==len(words)-1:
                        if word_idx==len(words)-1 and line_start_x_desc+f_leg.render(test_line_desc,1,txt_c).get_width()<=leg_w-pad: current_line_desc=test_line_desc
                        desc_s=f_leg.render(current_line_desc.strip(),1,txt_c); leg_s.blit(desc_s,(line_start_x_desc,temp_y_desc))
                        if not (word_idx==len(words)-1 and current_line_desc==test_line_desc): temp_y_desc+=line_h_leg; current_line_desc=word+" "
                        else: current_line_desc = "" # Ensure it doesn't carry over if it was the last word of the last line
                cur_y=max(initial_cur_y_for_item+line_h_leg, temp_y_desc+(line_h_leg if current_line_desc.strip() else 0) )
            cur_y+=sect_sp//2 
        self.scr.blit(leg_s,leg_r)
        if leg_r.bottom < WIN_H-int(10*SIZE_SCALE):
            close_s=self.font_tiny.render("Press L or Esc to close",1,(180,180,180)); self.scr.blit(close_s,close_s.get_rect(midbottom=(WIN_W//2,WIN_H-int(5*SIZE_SCALE))))


    def draw_mini_legend(self):
        """Draws a small, persistent legend of essential shortcuts."""
        if self.show_help: return # Don't draw if full legend is shown

        lines = [
            "S: Solve",
            "M: Info Pages",
            "L: Full Legend",
            "Esc: Cancel"
        ]
        line_height = self.font_tiny.get_linesize()
        padding = int(5 * SIZE_SCALE)
        
        y_pos = WIN_H - padding - (len(lines) * line_height) - ( (len(lines)-1) * int(2*SIZE_SCALE) )
        
        for i, text in enumerate(lines):
            surf = self.font_tiny.render(text, True, (180, 180, 180))
            rect = surf.get_rect(bottomright=(WIN_W - padding, y_pos + (i+1)*line_height + i*int(2*SIZE_SCALE) ))
            self.scr.blit(surf, rect)

    def draw_heat_map(self): # ... (same as v6)
        Y = self.circ.Y; b_vec = self.circ.b 
        if Y is None or b_vec is None or Y.shape[0] == 0: return
        dim = Y.shape[0]
        info_bar_content_y_start = CANVAS_H + int(40 * SIZE_SCALE) 
        available_h = WIN_H - info_bar_content_y_start - int(5 * SIZE_SCALE) 
        available_w = WIN_W - int(40 * SIZE_SCALE)
        max_cell_s = int(12*SIZE_SCALE) 
        cell_s = max(1, min(available_w // (dim + 1.5) if (dim+1.5)>0 else available_w, 
                             available_h // dim if dim>0 else available_h, 
                             max_cell_s))
        spacing = max(1, int(cell_s / 8)) 
        total_h = dim * (cell_s + spacing) - spacing if dim > 0 else 0
        total_w_y = dim * (cell_s + spacing) - spacing if dim > 0 else 0
        total_w_b = cell_s
        b_gap = spacing * 2 
        ox = (WIN_W - total_w_y - b_gap - total_w_b) // 2
        oy = info_bar_content_y_start + (available_h - total_h) // 2
        mY = np.max(np.abs(Y)) if np.any(Y) else 1.0; mY = max(mY, 1e-12)
        mb = np.max(np.abs(b_vec)) if np.any(b_vec) else 1.0; mb = max(mb, 1e-12)
        border_th = max(1, int(0.5*SIZE_SCALE)) 
        for r_iter, c_iter in itertools.product(range(dim), repeat=2): 
            v = np.abs(Y[r_iter,c_iter]); br_val = min(1, v/mY if mY!=0 else 0) 
            col_val = int(220*br_val) 
            cell_r = pg.Rect(ox+c_iter*(cell_s+spacing), oy+r_iter*(cell_s+spacing), cell_s, cell_s) 
            pg.draw.rect(self.scr, (col_val,col_val,col_val), cell_r)
            if cell_s > 3: pg.draw.rect(self.scr, (80,80,80), cell_r, border_th)
        bx_pos = ox + total_w_y + b_gap
        for r_idx, v_val_b in enumerate(b_vec): 
            v_abs = np.abs(v_val_b); br_val_b = min(1, v_abs/mb if mb!=0 else 0) 
            col_val_b = int(220*br_val_b) 
            cell_r_b = pg.Rect(bx_pos, oy+r_idx*(cell_s+spacing), cell_s, cell_s) 
            pg.draw.rect(self.scr, (int(col_val_b*0.6),col_val_b,int(col_val_b*0.6)), cell_r_b)
            if cell_s > 3: pg.draw.rect(self.scr, (80,80,80), cell_r_b, border_th)


    def draw_matrix_text(self): # Using the improved version from your previous response
        # ... (Paste the improved draw_matrix_text from the previous good response here) ...
        # ... (It was quite long, so assuming it's correctly pasted) ...
        Y, b_vec = self.circ.Y, self.circ.b 
        if Y is None or b_vec is None or Y.shape[0] == 0:
            no_data_surf = self.font_small.render("No matrix data. Solve circuit (S).", True, (200,200,200))
            self.scr.blit(no_data_surf, no_data_surf.get_rect(centerx=WIN_W//2, top=CANVAS_H + int(50*SIZE_SCALE)))
            return
        dim = Y.shape[0]
        def fmt_c(val):
            mag = abs(val)
            if mag < 1e-12: return "0"
            if abs(val.imag) < 1e-9 * max(abs(val.real), 1e-12): return f"{val.real:.2g}"
            if abs(val.real) < 1e-9 * max(abs(val.imag), 1e-12): return f"{val.imag:.2g}j"
            return f"{mag:.2g}\u2220{math.degrees(cmath.phase(val)):.0f}\u00B0"
        unknown_labels = [self.circ.idx_to_nodename.get(i, f"x{i}") for i in range(dim)]
        Y_str_elements = [[fmt_c(Y[r, c_idx]) for c_idx in range(dim)] for r in range(dim)]
        b_str_elements = [fmt_c(b_val) for b_val in b_vec]
        max_row_label_w = 0
        if unknown_labels: max_row_label_w = max(self.font_tiny.render(lbl + ": ", True, (0,0,0)).get_width() for lbl in unknown_labels)
        Y_col_widths = [0] * dim
        for c_idx in range(dim):
            header_w = self.font_tiny.render(unknown_labels[c_idx], True, (0,0,0)).get_width()
            Y_col_widths[c_idx] = header_w
            for r_idx in range(dim): Y_col_widths[c_idx] = max(Y_col_widths[c_idx], self.font_tiny.render(Y_str_elements[r_idx][c_idx], True, (0,0,0)).get_width())
        x_vec_col_width = 0
        if unknown_labels: x_vec_col_width = max(self.font_tiny.render(lbl, True, (0,0,0)).get_width() for lbl in unknown_labels)
        b_vec_col_width = 0
        if b_str_elements: b_vec_col_width = max(self.font_tiny.render(s, True, (0,0,0)).get_width() for s in b_str_elements)
        col_spacing = int(10 * SIZE_SCALE); op_spacing = int(8 * SIZE_SCALE); bracket_spacing = int(3*SIZE_SCALE)
        line_h = self.font_tiny.get_linesize() + int(2 * SIZE_SCALE) 
        info_bar_content_y_start = CANVAS_H + int(40 * SIZE_SCALE) 
        start_x = int(10 * SIZE_SCALE); current_y = info_bar_content_y_start
        current_x_col_header = start_x + max_row_label_w + bracket_spacing 
        for c_idx in range(dim):
            header_surf = self.font_tiny.render(unknown_labels[c_idx], True, (180, 180, 220)) 
            header_rect = header_surf.get_rect(centerx=current_x_col_header + Y_col_widths[c_idx] // 2, top=current_y)
            self.scr.blit(header_surf, header_rect)
            current_x_col_header += Y_col_widths[c_idx] + col_spacing
        current_y += line_h
        for r_idx in range(dim):
            if current_y + line_h > WIN_H - int(5*SIZE_SCALE): break 
            current_x_row = start_x
            row_label_surf = self.font_tiny.render(unknown_labels[r_idx] + ":", True, (180, 180, 220))
            row_label_rect = row_label_surf.get_rect(topleft=(current_x_row, current_y))
            self.scr.blit(row_label_surf, row_label_rect)
            current_x_row += max_row_label_w
            bracket_surf = self.font_tiny.render("|", True, (200,210,200))
            self.scr.blit(bracket_surf, (current_x_row, current_y))
            current_x_row += bracket_surf.get_width() 
            for c_idx in range(dim):
                elem_surf = self.font_tiny.render(Y_str_elements[r_idx][c_idx], True, (200,210,200))
                elem_rect = elem_surf.get_rect(right=current_x_row + Y_col_widths[c_idx], centery=current_y + line_h // 2 -1)
                self.scr.blit(elem_surf, elem_rect)
                current_x_row += Y_col_widths[c_idx] + col_spacing
            current_x_row -= col_spacing 
            self.scr.blit(bracket_surf, (current_x_row, current_y))
            current_x_row += bracket_surf.get_width() + op_spacing
            mul_surf = self.font_tiny.render("*", True, (220,220,100))
            self.scr.blit(mul_surf, (current_x_row, current_y))
            current_x_row += mul_surf.get_width() + op_spacing
            self.scr.blit(bracket_surf, (current_x_row, current_y))
            current_x_row += bracket_surf.get_width()
            x_elem_surf = self.font_tiny.render(unknown_labels[r_idx], True, (200,210,200))
            x_elem_rect = x_elem_surf.get_rect(centerx=current_x_row + x_vec_col_width // 2, centery=current_y + line_h // 2 -1)
            self.scr.blit(x_elem_surf, x_elem_rect)
            current_x_row += x_vec_col_width
            self.scr.blit(bracket_surf, (current_x_row, current_y))
            current_x_row += bracket_surf.get_width() + op_spacing
            eq_surf = self.font_tiny.render("=", True, (220,220,100))
            self.scr.blit(eq_surf, (current_x_row, current_y))
            current_x_row += eq_surf.get_width() + op_spacing
            self.scr.blit(bracket_surf, (current_x_row, current_y))
            current_x_row += bracket_surf.get_width()
            b_elem_surf = self.font_tiny.render(b_str_elements[r_idx], True, (180,220,180))
            b_elem_rect = b_elem_surf.get_rect(centerx=current_x_row + b_vec_col_width // 2, centery=current_y + line_h // 2 -1)
            self.scr.blit(b_elem_surf, b_elem_rect)
            current_x_row += b_vec_col_width
            self.scr.blit(bracket_surf, (current_x_row, current_y))
            current_y += line_h

    def draw_stamps_info(self): # ... (same as v6, ensure it's now called "Circuit Info") ...
        if self.circ.Y is None or self.circ.Y.shape[0] == 0: return
        lines = []
        lines.append(f"Nodes (Excl.GND): {self.circ.node_cnt}, V-Sources: {len(self.circ.vsrc_idx)}, Matrix Dim: {self.circ.Y.shape[0]}")
        lines.append(f"Frequency: {value_to_str(self.freq)}Hz")
        lines.append("-" * 40)
        lines.append("Solved Variables (MNA Unknowns):")
        for i in range(self.circ.Y.shape[0]):
            var_name = self.circ.idx_to_nodename.get(i, f"Var{i}")
            val_str = ""
            if self.circ.solution is not None and i < len(self.circ.solution):
                phasor = self.circ.solution[i]
                mag, ph_rad = cmath.polar(phasor); ph_deg = math.degrees(ph_rad)
                unit = "(V)" if var_name.startswith("V(") else "(A)" 
                val_str = " = 0" if mag < 1e-12 else f" = {value_to_str(mag)}\u2220{ph_deg:.0f}\u00B0 {unit}"
            lines.append(f"  x[{i}]: {var_name}{val_str}")
        lines.append("-" * 40)
        lines.append("Component Details & Connections:")
        omega = 2 * math.pi * self.freq
        for c_comp in self.circ.comps: # Renamed loop var
            pin_info_list = []
            for p_idx, p_pin in enumerate(c_comp.pins): # Renamed loop vars
                net_display = "GND" if p_pin.net == 0 else f"N{p_pin.net}"
                pin_info_list.append(f"{p_pin.name}:{net_display}")
            nets_str = ", ".join(pin_info_list)
            details_str = ""
            try:
                if c_comp.ct == CType.GND: details_str = "Ground Reference (Net 0)"
                elif c_comp.ct == CType.R: details_str = f"R={c_comp.label()}, G={1/max(1e-12,c_comp.val):.2g}S"
                elif c_comp.ct == CType.C:
                    if omega==0: details_str=f"C={c_comp.label()}, Yc=0S (DC Open)"
                    else: details_str=f"C={c_comp.label()}, Yc={1j*omega*c_comp.val:.2g}S"
                elif c_comp.ct == CType.L:
                    if omega==0: details_str=f"L={c_comp.label()}, Yl=inf S (DC Short)"
                    elif c_comp.val==0: details_str=f"L={c_comp.label()}, Yl=inf S (0H Short)"
                    else: details_str=f"L={c_comp.label()}, Yl={1/(1j*omega*c_comp.val):.2g}S"
                elif c_comp.ct in (CType.VDC,CType.VAC): details_str=f"V={c_comp.label()}"
            except (ZeroDivisionError, OverflowError): details_str="Math Error in calc"
            lines.append(f"  {c_comp.name} [{nets_str}]: {details_str}")
        line_h_info = self.font_tiny.get_linesize() 
        info_bar_content_y_start = CANVAS_H + int(40 * SIZE_SCALE)
        start_y_info = info_bar_content_y_start + int(2*SIZE_SCALE) 
        left_pad_info = int(10*SIZE_SCALE) 
        max_fit_lines = (WIN_H - start_y_info - int(5*SIZE_SCALE)) // line_h_info if line_h_info > 0 else 0
        for i, line_txt in enumerate(lines[:max_fit_lines]):
            surf = self.font_tiny.render(line_txt, True, (180, 220, 180))
            self.scr.blit(surf, (left_pad_info, start_y_info + i * line_h_info))


# ───────── main ─────────
if __name__=="__main__":
    try: App().run()
    except KeyboardInterrupt: pass
    finally: pg.quit()
