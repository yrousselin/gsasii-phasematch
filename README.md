# Phase search (search-match) for GSAS-II

*[Version française : LISEZMOI.md](LISEZMOI.md)*

A plugin for **GSAS-II** that helps identify the **crystalline phases present in a powder diffraction pattern**, including mixtures. It compares your pattern with a local index of reference structures built from the CSD, the COD and/or your own CIF files, and shows the best candidates as overlays on the pattern. **You** decide which phase to keep. The next search then looks at what remains to be explained.

## What it does

- Adds the menu entry **Export → Entire project as → (user) Phase search (search-match)** to GSAS-II.
- Opens an interactive window. **The interface is in French.** It shows:
  - the measured pattern after background subtraction;
  - the phases you have kept;
  - the selected candidate;
  - the line positions;
  - what remains to be explained.
- **How a search works.** First, a coarse screening of the whole index checks whether the strong lines of each reference are present in the pattern, allowing an isotropic expansion of the cell. Then a fine screening fits the cell to the observed peaks. Finally, the best candidates are validated with the full structure: cell, scale and peak width are fitted to the profile.
- When you keep a phase (**Retenir cette phase**), all kept phases are fitted together with non-negative scale factors, and the search starts again on the residual.
- Optional filter by **required elements**, for example `C H N O`, with a strict mode (**aucun autre**).
- **Voir dans Mercury** opens a candidate in Mercury (macOS).
- The kept phases are exported, with a summary saved as `<project>_phases.json`.

## Requirements

| What | Why | Where to get it |
|---|---|---|
| **GSAS-II** | hosts the plugin; provides numpy, scipy, matplotlib and wxPython | [GSAS-II](https://advancedphotonsource.github.io/GSAS-II-tutorials/install.html) |
| **A reference index** (built once, see below) | the structures to compare with | CSD, COD and/or your own CIF folders |
| CCDC software with the CSD Python API *(optional)* | adding the CSD to the index | [CCDC](https://www.ccdc.cam.ac.uk/) (requires a CSD licence) |
| gemmi *(optional)* | symmetry of COD files that have no symmetry operators | `pip install gemmi` |
| rsync *(optional)* | fast COD download and updates; built into macOS | otherwise use `--archive` |

Only **constant-wavelength** powder data (X-rays or neutrons) are supported.

## Installation

### Step 1: copy the files into GSAS-II

**macOS (and Linux)**
1. Download this repository: green **Code** button → **Download ZIP**, then double-click the ZIP to unpack it.
2. Open the **Terminal**. Type `bash `, with a space after it, drag `install_phasematch_plugin.sh` into the window, and press Enter. The files are copied into `~/.GSASII/exports/` and gemmi is installed.

**Windows**
1. Download the ZIP as above, then right-click it → **Extract All**.
2. Double-click **`install_windows.bat`**. The files are copied into `C:\Users\<your name>\.GSASII\exports\`.
3. Optional: install gemmi with GSAS-II's Python:
   `"<GSAS-II folder>\python.exe" -m pip install gemmi`

### Step 2: build the reference index (once)

The index is a single file, `~/.GSASII/pm_index.sqlite`. It stores only the strongest lines and a short description of each structure. The full structure is read again from its source when needed.

Run these commands in a terminal, with GSAS-II's Python. Below, `PY` stands for that Python, for example `~/g2main/bin/python` on macOS, and `EXP` for the folder `~/.GSASII/exports` (`%USERPROFILE%\.GSASII\exports` on Windows). Choose one or more sources:

```
# Your own CIF files (any folder, for example ICSD exports)
PY EXP/pm_index.py build --cif-dir /path/to/my_CIFs --source MyCIFs

# The CSD (requires the CCDC software)
PY EXP/pm_index.py build --csd

# The COD: first a local copy of the CIF files (long the first time), then the index
PY EXP/pm_index.py cod-sync              # macOS / Linux (rsync)
PY EXP/pm_index.py cod-sync --archive    # Windows, or if rsync is missing
PY EXP/pm_index.py build --cod

# Check what the index contains
PY EXP/pm_index.py stats
```

`build` is **incremental**. You can stop it with Ctrl-C and run it again later. After a CSD or COD update, run the same commands again: only the new entries are added.

If the CCDC's Python is not found, add `--csd-python` followed by its path. Typical locations:
- **Windows:** `"C:\Program Files\CCDC\ccdc-software\csd-python-api\miniconda\python.exe"` (to find it: `dir /s /b "C:\Program Files\CCDC\python.exe"`)
- **macOS:** `~/CCDC/ccdc-software/csd-python-api/miniconda/bin/python`

#### Full example on Windows

Replace `<name>` with your Windows user name, and the `python.exe` path with your GSAS-II one (here `C:\ProgramData\gsas2main`). Open the **Command Prompt**:

```
"C:\ProgramData\gsas2main\python.exe" C:\Users\<name>\.GSASII\exports\pm_index.py build --csd
"C:\ProgramData\gsas2main\python.exe" C:\Users\<name>\.GSASII\exports\pm_index.py cod-sync --archive
"C:\ProgramData\gsas2main\python.exe" C:\Users\<name>\.GSASII\exports\pm_index.py build --cod
"C:\ProgramData\gsas2main\python.exe" C:\Users\<name>\.GSASII\exports\pm_index.py stats
```

What to expect (COD, about 520,000 files):
- `cod-sync --archive` downloads an archive of more than 10 GB, then extracts the CIF files into `C:\Users\<name>\COD\cif` (about 35 min of extraction). Allow several tens of GB of free space. If the download is interrupted, run the same command again: it resumes where it stopped.
- `build --cod` shows its progress every 30 s, with an estimated time remaining. On a 12-core PC: about 120 structures per second, roughly 1 h 15 min.
- Prevent the PC from going to sleep during these steps (Settings → System → Power).
- Ctrl-C prints a long series of "KeyboardInterrupt" messages: this is normal, the work already done is kept.

**To check an installation**, first run a small test into a separate index (a few seconds):
```
"C:\ProgramData\gsas2main\python.exe" C:\Users\<name>\.GSASII\exports\pm_index.py build --cod --workers 1 --limit 300 --db C:\Users\<name>\test_cod.sqlite
```
It should end with "0 en erreur" (0 errors). Then delete `test_cod.sqlite`.

**If the build is unusually slow**, or a line "ATTENTION : … processus de calcul ont planté" (worker processes crashed) appears, run the test above: it shows exactly where any crash happens. On Windows, a `0xc06d007f` error in numpy came from the numerical DLLs of the GSAS-II Python, which are not found when `python.exe` is run directly. `pm_index.py` now fixes this itself: make sure you have the latest version.

### Step 3: restart GSAS-II

## Usage

1. Open a GSAS-II project containing a powder pattern.
2. Open **Export → Entire project as → (user) Phase search (search-match)**.
3. Look at the candidates and their overlays. Click **Retenir cette phase** (keep this phase) for a phase that matches. The search then continues on the remaining signal.
4. Repeat until everything is explained, then click **Exporter les phases retenues…** (export kept phases).

**From the command line (without the window):**

```
PY EXP/pm_search.py project.gpx --elements "C H N O" --max-phases 4
PY EXP/pm_search.py --help
```

## Limitations

- Constant-wavelength data only (no time-of-flight or energy-dispersive data).
- The result is a **proposal**. The final decision, and the confirmation by Rietveld refinement, are up to the user.
- The quality of the search depends on what is in the index.
- Developed and used on macOS. Windows should work, but has not been fully tested yet. Feedback is welcome.

## Credits and references

- **GSAS-II:** B. H. Toby, R. B. Von Dreele, *J. Appl. Cryst.* **2013**, 46, 544–549. The X-ray scattering factor coefficients in `pm_f0.py` are taken from GSAS-II.
- **COD:** Crystallography Open Database, <https://www.crystallography.net/>. **CSD:** Cambridge Crystallographic Data Centre. **gemmi:** <https://gemmi.readthedocs.io/>.
- Written by Yoann Rousselin with the help of Claude (Anthropic) as a coding assistant.

## Licence

MIT. See [LICENSE](LICENSE).

## Author

**Yoann Rousselin**, CNRS research engineer (X-ray crystallography), ICMUB, Dijon (France).
Questions or bug reports: please open an *Issue* on this page.
