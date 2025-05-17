# Interactive Circuit Canvas (v7)

## Description

The Interactive Circuit Canvas is a Python-based application that allows users to design, simulate, and analyze electronic circuits. It provides an intuitive graphical interface for placing components, wiring them together, and viewing simulation results, including DC and AC analysis. The simulation is powered by Modified Nodal Analysis (MNA) with LU decomposition for solving the circuit equations.

This project is built using Pygame for the graphical user interface and event handling, NumPy for numerical operations, and SciPy for matrix solving.

**Current Version:** v7 - Final Polish (as of 2025-04-20)

![Screenshot of Circuit Canvas in action (placeholder - replace with actual screenshot)](./screenshot.png)
*(Suggestion: Add a good screenshot of your application here named `screenshot.png` in the same directory as the README)*

## Features

*   **Interactive Canvas:**
    *   Drag-and-drop components from a palette onto a gridded canvas.
    *   Click-to-connect wiring between component pins.
    *   Move and delete components.
    *   Multi-select components (Ctrl+Click) for group actions (drag, delete).
*   **Component Support:**
    *   Resistors (R)
    *   Capacitors (C)
    *   Inductors (L)
    *   DC Voltage Sources (VDC)
    *   AC Voltage Sources (VAC) - with magnitude and phase.
    *   Ground (GND) - essential for simulation.
*   **Circuit Simulation:**
    *   **Modified Nodal Analysis (MNA):** Core algorithm for setting up circuit equations.
    *   **LU Decomposition:** Uses `scipy.linalg.lu_factor` and `scipy.linalg.lu_solve` for solving the system `Y*x = b`.
    *   **DC Analysis:** Solves the circuit with frequency set to 0 Hz.
    *   **AC Analysis:** Solves the circuit at a user-specified frequency, handling complex impedances and phasors.
*   **Real-time Feedback & Results:**
    *   **Solution Overlay:** Displays calculated node voltages (magnitude and phase) and component currents directly on the schematic after a successful solve.
    *   **Status Bar:** Provides continuous feedback on user actions, errors, and simulation status.
    *   **Info Pages (Cycle with 'M' key):**
        1.  **Heatmap:** Visual representation of the MNA matrix `Y` and knowns vector `b`.
        2.  **Matrix Text:** Detailed textual display of the `Y` matrix, unknowns vector `x` (labels), and `b` vector, showing the full `Y*x = b` equation structure.
        3.  **Circuit Info:** Lists solved variables (node voltages, source currents with descriptive names), component details (value, calculated stamp characteristic like G or Yc/Yl), and pin-to-net connections.
*   **User Interface Enhancements:**
    *   Scalable UI elements and fonts.
    *   Component value editing via keyboard input (with SI prefix parsing).
    *   Frequency setting for AC analysis via keyboard input.
    *   Wire deletion mode.
    *   Automatic snapping of wires when placing components near existing pins.
    *   Comprehensive pop-up Legend for controls and shortcuts.
    *   Persistent mini-legend for essential shortcuts visible at the bottom-right.
    *   Improved font rendering for special characters (Ω, ∠, °, µ) using font fallbacks.

## Requirements

*   Python 3.x
*   Pygame: `pip install pygame`
*   NumPy: `pip install numpy`
*   SciPy: `pip install scipy`

*(Optional, but recommended for best font rendering if system defaults are poor: DejaVu Sans Mono font installed, or other fonts listed in `FONT_*_NAME` constants in the script).*

## How to Run

1.  Ensure all required libraries are installed.
2.  Save the script as `app.py` (or your preferred name).
3.  Run from the command line:
    ```bash
    python app.py
    ```

## Keyboard Shortcuts

*(A brief summary - press 'L' in the app for the full, detailed legend)*

*   **Mouse:**
    *   **Drag Palette Item:** Place component.
    *   **Click Pin → Pin:** Create wire.
    *   **Drag Component:** Move selected component(s).
    *   **Right-Click Part:** Delete component under cursor.
    *   **Mouse Wheel over Part:** Adjust value (×10 / ÷10).
    *   **Shift + Mouse Wheel:** Fine adjust value (×√10 / ÷√10).
    *   **Ctrl + Click Part:** Toggle selection.
*   **General Keys:**
    *   **S:** Solve circuit.
    *   **Esc:** Cancel current action / Deselect all / Close prompts & legend.
    *   **L / H / ?:** Toggle full help Legend.
    *   **Delete / Backspace:** Delete selected component(s).
    *   **Ctrl + A:** Select all components.
*   **Modes & Settings:**
    *   **C:** Toggle Component Edit Mode (type new value).
    *   **F:** Set global AC frequency for analysis.
    *   **W:** Toggle Wire Delete Mode (click wire to delete).
*   **Information View:**
    *   **M:** Cycle through Info Pages (Solution Overlay, Heatmap, Matrix Text, Circuit Info).
    *   **E:** Show Circuit Info Page directly.

## Code Structure Overview

*   **UI Scaling Factors & Constants:** Define base sizes and scaling for UI elements.
*   **Part Types (`CType`, `LBL`, `COL`):** Enumeration and metadata for circuit components.
*   **Basic Classes (`Pin`, `Wire`, `Comp`):** Core data structures for representing circuit elements.
    *   `Comp` includes automatic naming (e.g., R1, C2).
*   **Value Parsing/Format (`parse_value`, `value_to_str`):** Utilities for handling numbers with SI prefixes.
*   **`Circuit` Class:**
    *   Manages components and wires.
    *   `_nets()`: Performs netlisting using Breadth-First Search (BFS) and identifies ground.
    *   `solve()`: Implements Modified Nodal Analysis (MNA).
        *   Calculates admittances for R, L, C components (handles DC and AC cases correctly).
        *   Builds the `Y` (MNA) matrix and `b` (knowns) vector.
        *   Uses `scipy.linalg.lu_factor` and `scipy.linalg.lu_solve` for solving the system.
        *   Populates `idx_to_nodename` for descriptive variable names.
*   **`App` Class:**
    *   Initializes Pygame, fonts, and UI elements.
    *   `run()`: Main application loop.
    *   `handle_events()`: Processes user input (keyboard, mouse).
    *   `draw*()` methods: Handle all rendering to the screen:
        *   `draw_grid()`, `draw_wires()`, `draw_comps()`, `draw_single_comp()`, `draw_solution_overlay()`
        *   `draw_palette()`, `draw_info_bar()`, `draw_legend()`, `draw_mini_legend()`
        *   `draw_heat_map()`, `draw_matrix_text()`, `draw_stamps_info()` (Circuit Info)
    *   Manages application state (modes, selections, prompts).

## Potential Future Enhancements

*   Save/Load circuit designs to/from file (e.g., JSON, custom format).
*   Undo/Redo functionality.
*   Zoom and Pan for the canvas.
*   Component rotation and flipping.
*   More advanced components (e.g., dependent sources, diodes, transistors, op-amps).
*   Transient analysis (time-domain simulation).
*   Plotting of results (e.g., frequency response, transient waveforms).
*   More robust error reporting for singular matrices (identifying floating nodes, etc.).
*   Bundling a `.ttf` font file for consistent cross-platform rendering of special characters.
*   In-place editing of component values directly on the canvas or via a properties panel.

## Known Issues / Considerations

*   Font rendering for special symbols (Ω, ∠, µ, °) relies on system-installed fonts specified in `FONT_*_NAME`. If suitable fonts are not found, Pygame's default font will be used, which may not render these characters correctly. Bundling a font would be a more robust solution.
*   Performance for very large circuits might degrade due to linear searches in event handling and O(N^3) complexity of matrix solving.
*   The "Matrix Text" view can still become very wide for larger matrices, potentially exceeding screen width even with `font_tiny`.

## License

Or: This project is open source and free to use.
