# Run once if needed:
# !pip install ipympl ipywidgets numpy matplotlib

%matplotlib widget

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import ipywidgets as widgets
from IPython.display import display
import csv
import os
from itertools import combinations_with_replacement

# ─────────────────────────────────────────────
#  PASTE YOUR MONOISOTOPIC MASSES HERE
# ─────────────────────────────────────────────

x_raw = """

"""

# ─────────────────────────────────────────────
#  PASTE THEORETICAL MASSES FROM PROSIGHT HERE (optional)
#  These are theoretical c/z-ion monoisotopic masses - H.
#  This is NOT part of the de novo sequencing method.
#  This tests mass error in experimental fragment masses.
#  Leave blank if not needed.
#  1.007825 will be added to each value automatically.
# ─────────────────────────────────────────────

theo_raw = """

"""

# ─────────────────────────────────────────────
#  Amino acid residue masses (monoisotopic, Da)
# ─────────────────────────────────────────────
AA_RESIDUES = {
    'G':   57.021464,
    'A':   71.037114,
    'S':   87.032028,
    'P':   97.052764,
    'V':   99.068414,
    'T':  101.047679,
    'C':  103.009185,
    'L/I':113.084064,
    'N':  114.042927,
    'D':  115.026943,
    'Q':  128.058578,
    'K':  128.094963,
    'E':  129.042593,
    'M':  131.040485,
    'H':  137.058912,
    'F':  147.068414,
    'R':  156.101111,
    'Y':  163.063329,
    'W':  186.079313,
    'C-H': 102.001360,   # Cysteine - H (103.009185 - 1.007825) — disulphide-bridged
}

AA_PAIRS = []
aa_names  = list(AA_RESIDUES.keys())
aa_masses = list(AA_RESIDUES.values())
for i, j in combinations_with_replacement(range(len(aa_names)), 2):
    AA_PAIRS.append((f'{aa_names[i]}+{aa_names[j]}', aa_masses[i] + aa_masses[j]))

AA_TRIPLES = []
for i, j, k in combinations_with_replacement(range(len(aa_names)), 3):
    AA_TRIPLES.append((f'{aa_names[i]}+{aa_names[j]}+{aa_names[k]}', aa_masses[i] + aa_masses[j] + aa_masses[k]))

PPM_THRESHOLD = 10.0
H2O = 18.010565

# ─────────────────────────────────────────────
#  Parse & calculate Y
# ─────────────────────────────────────────────
def parse_data(raw):
    values = []
    for token in raw.replace('\t', '\n').split('\n'):
        token = token.strip()
        if token:
            try:
                values.append(float(token))
            except ValueError:
                pass
    return np.array(values)

x = parse_data(x_raw)
averagine_scaled = 0.999493892 * x
nominal_mass     = np.round(averagine_scaled)
y                = nominal_mass - averagine_scaled

# Parse theoretical masses and add 1.007825
_theo_raw_parsed = parse_data(theo_raw)
theo_masses = _theo_raw_parsed + 1.007825 if len(_theo_raw_parsed) > 0 else np.array([])

# ─────────────────────────────────────────────
#  Constants
# ─────────────────────────────────────────────
x_range = float(x.max() - x.min()) or 1.0
y_range = float(y.max() - y.min()) or 1.0

X_PAD  = x_range * 0.05
X_LEFT = float(x.min() - X_PAD)
X_RIGHT= float(x.max() + X_PAD)
Y_PAD  = y_range * 0.12

TOL_X  = x_range * 0.015
TOL_Y  = y_range * 0.04

sections = []
drag     = {'active': False, 'type': None, 'index': None}

# ─────────────────────────────────────────────
#  Figure
# ─────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 6))
fig.subplots_adjust(left=0.12, right=0.97, top=0.97, bottom=0.12)

sc = ax.scatter(x, y, c='grey', marker='s', s=60, zorder=3)

ax.set_xlabel('Monoisotopic Mass (Da)', fontsize=12, labelpad=8)
ax.set_ylabel('Averagine Scaled Mass Defect (Da)', fontsize=12, labelpad=8)
ax.grid(True, linestyle=':', alpha=0.5)
ax.set_xlim(X_LEFT, X_RIGHT)
ax.set_ylim(y.min() - Y_PAD, y.max() + Y_PAD)

red_patch   = mpatches.Patch(color='red',   label='c-ions')
black_patch = mpatches.Patch(color='black', label='z-ions')
grey_patch  = mpatches.Patch(color='grey',  label='unassigned')
ax.legend(handles=[red_patch, black_patch, grey_patch],
          loc='upper right', framealpha=0.85, fontsize=10)

# ─────────────────────────────────────────────
#  Section helpers
# ─────────────────────────────────────────────
def left_of(i):
    return X_LEFT if i == 0 else sections[i - 1]['x_right']

def update_hline_spans():
    for i, sec in enumerate(sections):
        sec['hline'].set_xdata([left_of(i), sec['x_right']])
        sec['hline'].set_ydata([sec['h_pos'], sec['h_pos']])
    fig.canvas.draw_idle()

def update_colors():
    colors = ['grey'] * len(x)
    for idx, (xi, yi) in enumerate(zip(x, y)):
        for i, sec in enumerate(sections):
            if left_of(i) <= xi <= sec['x_right']:
                colors[idx] = 'red' if yi >= sec['h_pos'] else 'black'
                break
    sc.set_facecolor(colors)
    fig.canvas.draw_idle()

# ─────────────────────────────────────────────
#  Add Section
# ─────────────────────────────────────────────
def add_section(_):
    n      = len(sections)
    x_left = sections[-1]['x_right'] if sections else X_LEFT
    new_x  = (x_left + X_RIGHT) / 2
    h_new  = float(np.mean(y))

    vline = ax.axvline(x=new_x, color='steelblue', lw=1.8, linestyle='-', zorder=5)
    hline, = ax.plot([left_of(n), new_x], [h_new, h_new],
                     color='royalblue', lw=1.8, linestyle='--', zorder=5)

    sections.append({'x_right': new_x, 'h_pos': h_new,
                     'vline': vline, 'hline': hline})
    update_colors()

# ─────────────────────────────────────────────
#  Mouse interaction
# ─────────────────────────────────────────────
def on_press(event):
    if event.inaxes != ax or event.button != 1:
        return
    ex, ey = event.xdata, event.ydata
    for i, sec in enumerate(sections):
        if abs(ex - sec['x_right']) < TOL_X:
            drag.update(active=True, type='vline', index=i)
            return
    for i, sec in enumerate(sections):
        if left_of(i) <= ex <= sec['x_right'] and abs(ey - sec['h_pos']) < TOL_Y:
            drag.update(active=True, type='hline', index=i)
            return

def on_motion(event):
    if not drag['active'] or event.inaxes != ax:
        return
    ex, ey = event.xdata, event.ydata
    i = drag['index']

    if drag['type'] == 'vline':
        lo = (sections[i - 1]['x_right'] if i > 0 else X_LEFT) + TOL_X * 2
        hi = (sections[i + 1]['x_right'] if i + 1 < len(sections) else X_RIGHT) - TOL_X * 2
        new_x = float(np.clip(ex, lo, hi))
        sections[i]['x_right'] = new_x
        sections[i]['vline'].set_xdata([new_x, new_x])
        update_hline_spans()
        update_colors()

    elif drag['type'] == 'hline':
        sections[i]['h_pos'] = float(ey)
        update_hline_spans()
        update_colors()

def on_release(event):
    drag.update(active=False, type=None, index=None)

fig.canvas.mpl_connect('button_press_event',   on_press)
fig.canvas.mpl_connect('motion_notify_event',  on_motion)
fig.canvas.mpl_connect('button_release_event', on_release)

# ─────────────────────────────────────────────
#  Shared analysis function (used for both ion types)
# ─────────────────────────────────────────────
def analyse_ions(masses, defects, labels):
    """
    For a sorted list of ion masses, compute single-residue matches
    and best 2- and 3-residue combination matches for each ion.
    Returns (rows, max_single_matches).
    """
    rows = []
    max_single = 0

    for i in range(len(masses)):
        # Single residue matches
        single_matches = []
        for j in range(i + 1, len(masses)):
            exp_diff = masses[j] - masses[i]
            for aa, theo_mass in AA_RESIDUES.items():
                ppm = ((exp_diff - theo_mass) / theo_mass) * 1e6
                if abs(ppm) <= PPM_THRESHOLD:
                    single_matches.append({
                        'fragment': labels[j],
                        'aa':       aa,
                        'ppm':      round(ppm, 4)
                    })
        single_matches.sort(key=lambda m: abs(m['ppm']))
        single_matches = single_matches[:3]
        max_single = max(max_single, len(single_matches))

        # 2-residue combination matches – keep only best
        pair_candidates = []
        for j in range(i + 1, len(masses)):
            exp_diff = masses[j] - masses[i]
            for pair_name, pair_mass in AA_PAIRS:
                ppm = ((exp_diff - pair_mass) / pair_mass) * 1e6
                if abs(ppm) <= PPM_THRESHOLD:
                    pair_candidates.append({
                        'fragment': labels[j],
                        'pair':     pair_name,
                        'ppm':      round(ppm, 4)
                    })
        pair_candidates.sort(key=lambda m: abs(m['ppm']))

        # 3-residue combination matches – keep only best
        triple_candidates = []
        for j in range(i + 1, len(masses)):
            exp_diff = masses[j] - masses[i]
            for triple_name, triple_mass in AA_TRIPLES:
                ppm = ((exp_diff - triple_mass) / triple_mass) * 1e6
                if abs(ppm) <= PPM_THRESHOLD:
                    triple_candidates.append({
                        'fragment': labels[j],
                        'triple':   triple_name,
                        'ppm':      round(ppm, 4)
                    })
        triple_candidates.sort(key=lambda m: abs(m['ppm']))

        # Theoretical mass match – check bare, +H2O, +2H2O
        theo_match = None
        if len(theo_masses) > 0:
            candidates = []
            for adduct_label, offset in [('', 0), ('(+H2O)', H2O), ('(+2H2O)', 2 * H2O)]:
                adjusted = theo_masses + offset
                ppms = ((adjusted - masses[i]) / masses[i]) * 1e6
                for k in range(len(theo_masses)):
                    if abs(ppms[k]) <= PPM_THRESHOLD:
                        candidates.append({
                            'theo_mass': round(theo_masses[k], 6),
                            'adduct':    adduct_label,
                            'ppm':       round(float(ppms[k]), 4)
                        })
            if candidates:
                theo_match = min(candidates, key=lambda t: abs(t['ppm']))

        rows.append({
            'label':        labels[i],
            'mass':         round(masses[i], 6),
            'defect':       round(defects[i], 6),
            'singles':      single_matches,
            'best_pair':    pair_candidates[0] if pair_candidates else None,
            'best_triple':  triple_candidates[0] if triple_candidates else None,
            'theo_match':   theo_match,
        })

    return rows, max_single

# ─────────────────────────────────────────────
#  Calculate & export to CSV
# ─────────────────────────────────────────────
def calculate(_):
    # ── Gather c-ions and z-ions ──
    c_indices, z_indices = [], []
    for idx, (xi, yi) in enumerate(zip(x, y)):
        for i, sec in enumerate(sections):
            if left_of(i) <= xi <= sec['x_right']:
                if yi >= sec['h_pos']:
                    c_indices.append(idx)
                else:
                    z_indices.append(idx)
                break

    if len(c_indices) < 2 and len(z_indices) < 2:
        print("Need at least 2 assigned ions before calculating.")
        return

    c_indices.sort(key=lambda i: x[i])
    z_indices.sort(key=lambda i: x[i])

    c_masses  = [x[i] for i in c_indices]
    c_defects = [y[i] for i in c_indices]
    c_labels  = [f'{n+1}' for n in range(len(c_indices))]

    z_masses  = [x[i] for i in z_indices]
    z_defects = [y[i] for i in z_indices]
    z_labels  = [f'{n+1}' for n in range(len(z_indices))]

    # ── Analyse both ion series ──
    c_rows, c_max_single = analyse_ions(c_masses, c_defects, c_labels)
    z_rows, z_max_single = analyse_ions(z_masses, z_defects, z_labels)
    max_single = max(c_max_single, z_max_single)

    # ── Build header ──
    header = ['Ion', 'Monoisotopic Mass (Da)', 'Mass Defect (Da)']
    for m in range(1, max_single + 1):
        header += [f'Match {m} Fragment', f'Match {m} Residue', f'Match {m} PPM Error']
    header += ['Best 2-Residue Fragment', 'Best 2-Residue Pair', 'Best 2-Residue PPM Error']
    header += ['Best 3-Residue Fragment', 'Best 3-Residue Combination', 'Best 3-Residue PPM Error']
    if len(theo_masses) > 0:
        header += ['Theoretical Mass Match (Da)', 'Adduct', 'Theoretical Mass PPM Error']

    def build_row(row):
        line = [row['label'], row['mass'], row['defect']]
        for match in row['singles']:
            line += [match['fragment'], match['aa'], match['ppm']]
        line += [''] * (max_single - len(row['singles'])) * 3
        if row['best_pair']:
            line += [row['best_pair']['fragment'],
                     row['best_pair']['pair'],
                     row['best_pair']['ppm']]
        else:
            line += ['', '', '']
        if row['best_triple']:
            line += [row['best_triple']['fragment'],
                     row['best_triple']['triple'],
                     row['best_triple']['ppm']]
        else:
            line += ['', '', '']
        if len(theo_masses) > 0:
            if row['theo_match']:
                line += [row['theo_match']['theo_mass'],
                         row['theo_match']['adduct'],
                         row['theo_match']['ppm']]
            else:
                line += ['', '', '']
        return line

    # ── Write CSV ──
    out_path = os.path.join(os.getcwd(), 'de_novo_sequencing.csv')
    with open(out_path, 'w', newline='', encoding='utf-8-sig') as f:
        writer = csv.writer(f)
        writer.writerow(header)

        if c_rows:
            writer.writerow(['--- c-ions ---'] + [''] * (len(header) - 1))
            for row in c_rows:
                writer.writerow(build_row(row))

        writer.writerow([])

        if z_rows:
            writer.writerow(['--- z-ions ---'] + [''] * (len(header) - 1))
            for row in z_rows:
                writer.writerow(build_row(row))

    print(f"Saved: {out_path}")
    print(f"c-ions: {len(c_labels)} | z-ions: {len(z_labels)}")

# ─────────────────────────────────────────────
#  Buttons
# ─────────────────────────────────────────────
btn_add = widgets.Button(
    description  = '+ Add Section',
    button_style = 'info',
    layout       = widgets.Layout(width='160px', height='36px')
)
btn_calc = widgets.Button(
    description  = 'Calculate',
    button_style = 'success',
    layout       = widgets.Layout(width='160px', height='36px')
)

btn_add.on_click(add_section)
btn_calc.on_click(calculate)

display(widgets.HBox([btn_add, btn_calc]))
plt.show()
