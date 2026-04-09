#!/usr/bin/env python3
"""
Interactive Plotly visualizations for the Poisson-Laplace solver output.

Reads JSON files exported by the C++ solver and generates:
  1. Potential heatmap with equipotential contours
  2. Electric field vector plot (quiver)
  3. E-field magnitude heatmap
  4. Solver convergence comparison
  5. Line profile cuts through the potential

Usage:
    python visualize.py <field_data.json> [--convergence conv.json] [--output-dir plots/]
    python visualize.py --benchmark benchmark.json
"""

import json
import argparse
import numpy as np
import os

# Try plotly, fall back to instructions
try:
    import plotly.graph_objects as go
    import plotly.subplots as sp
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    print("Plotly not found. Install with: pip install plotly")
    print("Generating static fallback with matplotlib instead.")

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def load_field_data(path):
    """Load field JSON exported by the C++ solver."""
    with open(path) as f:
        data = json.load(f)

    phi = np.array(data["phi"])
    Ex = np.array(data["Ex"])
    Ey = np.array(data["Ey"])

    nx, ny = data["nx"], data["ny"]
    x = np.linspace(data["x0"], data["x1"], nx)
    y = np.linspace(data["y0"], data["y1"], ny)

    return x, y, phi, Ex, Ey


def load_convergence(path):
    """Load convergence JSON."""
    with open(path) as f:
        data = json.load(f)
    return data["residuals"]


def load_benchmark(path):
    """Load benchmark JSON with multiple solver results."""
    with open(path) as f:
        data = json.load(f)
    return data["entries"]


# ── Plotly Visualizations ─────────────────────────────────────

def plot_potential_plotly(x, y, phi, title="Scalar Potential φ"):
    """Interactive heatmap of the scalar potential with contours."""
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("Potential Heatmap", "Equipotential Contours"),
                        horizontal_spacing=0.12)

    # Heatmap
    fig.add_trace(go.Heatmap(
        z=phi.T, x=x * 1000, y=y * 1000,
        colorscale="RdBu_r",
        colorbar=dict(title="φ (V)", x=0.45),
        zmid=0
    ), row=1, col=1)

    # Contour
    fig.add_trace(go.Contour(
        z=phi.T, x=x * 1000, y=y * 1000,
        colorscale="RdBu_r",
        contours=dict(showlabels=True,
                      labelfont=dict(size=10, color="black")),
        ncontours=20,
        showscale=False
    ), row=1, col=2)

    fig.update_xaxes(title_text="x (mm)", scaleanchor="y", scaleratio=1)
    fig.update_yaxes(title_text="y (mm)")
    fig.update_layout(title=title, height=500, width=1100)
    return fig


def plot_efield_plotly(x, y, Ex, Ey, title="Electric Field"):
    """Quiver plot of E-field + magnitude heatmap."""
    E_mag = np.sqrt(Ex**2 + Ey**2)

    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("|E| Magnitude", "E-field Vectors"),
                        horizontal_spacing=0.12)

    # Magnitude heatmap
    fig.add_trace(go.Heatmap(
        z=E_mag.T, x=x * 1000, y=y * 1000,
        colorscale="Inferno",
        colorbar=dict(title="|E| (V/m)", x=0.45)
    ), row=1, col=1)

    # Quiver: subsample for clarity
    step = max(1, len(x) // 20)
    xs = x[::step]
    ys = y[::step]
    Exs = Ex[::step, ::step]
    Eys = Ey[::step, ::step]
    Es = np.sqrt(Exs**2 + Eys**2)
    Es[Es == 0] = 1  # Avoid division by zero

    # Normalize arrows
    scale = (x[-1] - x[0]) / 20 * 1000
    Uxn = Exs / Es * scale
    Uyn = Eys / Es * scale

    # Draw arrows as lines
    arrow_x = []
    arrow_y = []
    for i in range(len(xs)):
        for j in range(len(ys)):
            x0 = xs[i] * 1000
            y0 = ys[j] * 1000
            arrow_x.extend([x0, x0 + Uxn[i, j], None])
            arrow_y.extend([y0, y0 + Uyn[i, j], None])

    fig.add_trace(go.Scatter(
        x=arrow_x, y=arrow_y,
        mode="lines",
        line=dict(color="white", width=1),
        showlegend=False
    ), row=1, col=2)

    # Background: magnitude
    fig.add_trace(go.Heatmap(
        z=E_mag.T, x=x * 1000, y=y * 1000,
        colorscale="Inferno",
        showscale=False,
        opacity=0.7
    ), row=1, col=2)

    fig.update_xaxes(title_text="x (mm)", scaleanchor="y", scaleratio=1)
    fig.update_yaxes(title_text="y (mm)")
    fig.update_layout(title=title, height=500, width=1100)
    return fig


def plot_convergence_plotly(residuals_dict, title="Solver Convergence"):
    """Plot convergence curves for one or more solvers."""
    fig = go.Figure()

    colors = ["#636EFA", "#EF553B", "#00CC96", "#AB63FA"]

    for idx, (name, residuals) in enumerate(residuals_dict.items()):
        iters = [i * 10 for i in range(len(residuals))]  # Sampled every 10
        fig.add_trace(go.Scatter(
            x=iters, y=residuals,
            mode="lines",
            name=name,
            line=dict(color=colors[idx % len(colors)], width=2)
        ))

    fig.update_layout(
        title=title,
        xaxis_title="Iteration",
        yaxis_title="L2 Residual Norm",
        yaxis_type="log",
        height=450, width=700,
        legend=dict(x=0.7, y=0.95)
    )
    return fig


def plot_benchmark_plotly(entries):
    """Bar chart comparing solver performance from benchmark data."""
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("Iterations to Converge",
                                        "Wall-clock Time (ms)"))

    solvers = sorted(set(e["solver"] for e in entries))
    grid_sizes = sorted(set(e["n"] for e in entries))
    colors = {"Jacobi": "#636EFA", "Gauss-Seidel": "#EF553B", "SOR": "#00CC96"}

    for solver in solvers:
        solver_entries = [e for e in entries if e["solver"] == solver]
        solver_entries.sort(key=lambda e: e["n"])

        x_labels = [f"N={e['n']}" for e in solver_entries]
        iterations = [e["iterations"] for e in solver_entries]
        times = [e["elapsed_ms"] for e in solver_entries]

        fig.add_trace(go.Bar(
            x=x_labels, y=iterations,
            name=solver, marker_color=colors.get(solver, "#999"),
            showlegend=True
        ), row=1, col=1)

        fig.add_trace(go.Bar(
            x=x_labels, y=times,
            name=solver, marker_color=colors.get(solver, "#999"),
            showlegend=False
        ), row=1, col=2)

    fig.update_layout(
        title="Solver Benchmark: Quadrupole Problem",
        barmode="group",
        height=450, width=1000
    )
    fig.update_yaxes(title_text="Iterations", row=1, col=1)
    fig.update_yaxes(title_text="Time (ms)", row=1, col=2)

    # Also plot convergence curves
    fig_conv = go.Figure()
    for entry in entries:
        if entry["n"] == max(grid_sizes):
            iters = [i * 10 for i in range(len(entry["residuals"]))]
            fig_conv.add_trace(go.Scatter(
                x=iters, y=entry["residuals"],
                mode="lines", name=f"{entry['solver']} (N={entry['n']})",
                line=dict(color=colors.get(entry["solver"], "#999"), width=2)
            ))

    fig_conv.update_layout(
        title=f"Convergence Comparison (N={max(grid_sizes)})",
        xaxis_title="Iteration",
        yaxis_title="L2 Residual",
        yaxis_type="log",
        height=450, width=700
    )

    return fig, fig_conv


def plot_line_cut_plotly(x, y, phi, title="Potential Line Cuts"):
    """Horizontal and vertical line cuts through the center."""
    fig = make_subplots(rows=1, cols=2,
                        subplot_titles=("Horizontal cut (y=0)",
                                        "Vertical cut (x=0)"))
    mid_j = len(y) // 2
    mid_i = len(x) // 2

    fig.add_trace(go.Scatter(
        x=x * 1000, y=phi[:, mid_j],
        mode="lines", name="φ(x, y=0)",
        line=dict(color="#636EFA", width=2)
    ), row=1, col=1)

    fig.add_trace(go.Scatter(
        x=y * 1000, y=phi[mid_i, :],
        mode="lines", name="φ(x=0, y)",
        line=dict(color="#EF553B", width=2)
    ), row=1, col=2)

    fig.update_xaxes(title_text="x (mm)", row=1, col=1)
    fig.update_xaxes(title_text="y (mm)", row=1, col=2)
    fig.update_yaxes(title_text="φ (V)")
    fig.update_layout(title=title, height=400, width=1000)
    return fig


# ── Matplotlib fallback ───────────────────────────────────────

def plot_potential_mpl(x, y, phi, output_dir, prefix):
    """Matplotlib fallback for potential visualization."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    im = ax1.pcolormesh(x * 1000, y * 1000, phi.T, cmap="RdBu_r", shading="auto")
    ax1.set_xlabel("x (mm)")
    ax1.set_ylabel("y (mm)")
    ax1.set_title("Potential Heatmap")
    ax1.set_aspect("equal")
    plt.colorbar(im, ax=ax1, label="φ (V)")

    cs = ax2.contour(x * 1000, y * 1000, phi.T, levels=20, cmap="RdBu_r")
    ax2.clabel(cs, fontsize=8)
    ax2.set_xlabel("x (mm)")
    ax2.set_ylabel("y (mm)")
    ax2.set_title("Equipotential Contours")
    ax2.set_aspect("equal")

    plt.tight_layout()
    path = os.path.join(output_dir, f"{prefix}_potential.png")
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"Saved: {path}")


def plot_efield_mpl(x, y, Ex, Ey, output_dir, prefix):
    """Matplotlib fallback for E-field."""
    E_mag = np.sqrt(Ex**2 + Ey**2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    im = ax1.pcolormesh(x * 1000, y * 1000, E_mag.T, cmap="inferno", shading="auto")
    ax1.set_xlabel("x (mm)")
    ax1.set_ylabel("y (mm)")
    ax1.set_title("|E| Magnitude")
    ax1.set_aspect("equal")
    plt.colorbar(im, ax=ax1, label="|E| (V/m)")

    step = max(1, len(x) // 20)
    xs, ys = np.meshgrid(x[::step] * 1000, y[::step] * 1000, indexing="ij")
    Exs = Ex[::step, ::step]
    Eys = Ey[::step, ::step]

    ax2.quiver(xs, ys, Exs, Eys, pivot="mid", scale_units="xy")
    ax2.set_xlabel("x (mm)")
    ax2.set_ylabel("y (mm)")
    ax2.set_title("E-field Vectors")
    ax2.set_aspect("equal")

    plt.tight_layout()
    path = os.path.join(output_dir, f"{prefix}_efield.png")
    plt.savefig(path, dpi=150)
    plt.close()
    print(f"Saved: {path}")


# ── Main ──────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Visualize Poisson-Laplace solver output")
    parser.add_argument("field_data", nargs="?", help="Field data JSON file")
    parser.add_argument("--convergence", "-c", nargs="*", help="Convergence JSON file(s)")
    parser.add_argument("--benchmark", "-b", help="Benchmark JSON file")
    parser.add_argument("--output-dir", "-o", default="plots", help="Output directory")
    parser.add_argument("--prefix", "-p", default="field", help="Output file prefix")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    if args.benchmark:
        entries = load_benchmark(args.benchmark)
        if HAS_PLOTLY:
            fig_bar, fig_conv = plot_benchmark_plotly(entries)
            path1 = os.path.join(args.output_dir, "benchmark_bars.html")
            path2 = os.path.join(args.output_dir, "benchmark_convergence.html")
            fig_bar.write_html(path1)
            fig_conv.write_html(path2)
            print(f"Saved: {path1}")
            print(f"Saved: {path2}")
        return

    if not args.field_data:
        parser.error("field_data is required unless --benchmark is used")

    x, y, phi, Ex, Ey = load_field_data(args.field_data)
    prefix = args.prefix

    if HAS_PLOTLY:
        # Interactive HTML plots
        fig = plot_potential_plotly(x, y, phi)
        path = os.path.join(args.output_dir, f"{prefix}_potential.html")
        fig.write_html(path)
        print(f"Saved: {path}")

        fig = plot_efield_plotly(x, y, Ex, Ey)
        path = os.path.join(args.output_dir, f"{prefix}_efield.html")
        fig.write_html(path)
        print(f"Saved: {path}")

        fig = plot_line_cut_plotly(x, y, phi)
        path = os.path.join(args.output_dir, f"{prefix}_linecuts.html")
        fig.write_html(path)
        print(f"Saved: {path}")

    elif HAS_MPL:
        plot_potential_mpl(x, y, phi, args.output_dir, prefix)
        plot_efield_mpl(x, y, Ex, Ey, args.output_dir, prefix)

    # Convergence
    if args.convergence:
        residuals_dict = {}
        for cpath in args.convergence:
            name = os.path.splitext(os.path.basename(cpath))[0]
            residuals_dict[name] = load_convergence(cpath)

        if HAS_PLOTLY:
            fig = plot_convergence_plotly(residuals_dict)
            path = os.path.join(args.output_dir, f"{prefix}_convergence.html")
            fig.write_html(path)
            print(f"Saved: {path}")


if __name__ == "__main__":
    main()
