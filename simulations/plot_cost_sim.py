
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter


# ============================================================
# PREFIXES
# ============================================================

prefixes = {
    "church": "C",
    "wukong": "W",
    "goldman": "G",
    "gcplus": "GC",
    "hedges": "H",
    "max_density": "MD",
    "no_homopolymer": "NH",
}


# ============================================================
# COLOUR FOR EACH ENCODING
# ============================================================

encoding_colors = {
    "church": "tab:blue",
    "wukong": "tab:orange",
    "goldman": "tab:green",
    "gcplus": "tab:red",
    "hedges": "tab:purple",
    "max_density": "tab:brown",
    "no_homopolymer": "tab:pink",
}


# ============================================================
# CREATE LABEL FROM CSV CONFIGURATION STRING
# ============================================================

def make_configuration_label(encoding, configuration):

    prefix = prefixes.get(encoding, encoding)

    if encoding in ["church", "wukong"]:

        # Example:
        # add_redundancy=True_rs_num=3
        match = re.search(
            r"add_redundancy=(True|False)_rs_num=(\d+)",
            configuration
        )

        if match:
            redundancy = match.group(1) == "True"
            rs_num = match.group(2)

            return f"{prefix}{int(redundancy)}-{rs_num}"

    elif encoding == "goldman":

        # Example:
        # mean=1
        match = re.search(
            r"mean=([^_]+)",
            configuration
        )

        if match:
            return f"G{match.group(1)}"

    elif encoding == "gcplus":

        # Example:
        # gcplus_c1=3
        match = re.search(
            r"gcplus_c1=([^_]+)",
            configuration
        )

        if match:
            return f"GC{match.group(1)}"

    elif encoding == "hedges":

        # Example:
        # hedges_coderate=3
        match = re.search(
            r"hedges_coderate=([^_]+)",
            configuration
        )

        if match:
            return f"H{match.group(1)}"

    elif encoding in ["max_density", "no_homopolymer"]:

        # Example:
        # reed_solo_percentage=0.7
        match = re.search(
            r"reed_solo_percentage=([0-9.]+)",
            configuration
        )

        if match:
            percentage = float(match.group(1))

            return f"{prefix}{percentage:g}"

    return prefix


# ============================================================
# PLOT COST VS ERROR TOLERANCE
# ============================================================

def plot_cost_vs_tolerance_csv(
    csv_file,
    output_dir,
    target_recovery
):

    # --------------------------------------------------------
    # Load CSV
    # --------------------------------------------------------

    df = pd.read_csv(csv_file)

    # --------------------------------------------------------
    # Make sure numeric columns are numeric
    # --------------------------------------------------------

    df["target_recovery"] = pd.to_numeric(
        df["target_recovery"]
    )

    df["dna_cost_nt_per_bit"] = pd.to_numeric(
        df["dna_cost_nt_per_bit"]
    )

    df["error_tolerance"] = pd.to_numeric(
        df["error_tolerance"],
        errors="coerce"
    )

    # --------------------------------------------------------
    # Select requested recovery
    # --------------------------------------------------------

    target_df = df[
        df["target_recovery"] == target_recovery
    ].copy()

    if target_df.empty:
        print(
            f"No data found for "
            f"{target_recovery * 100:g}% recovery."
        )
        return

    # --------------------------------------------------------
    # Create figure
    # --------------------------------------------------------

    plt.figure(figsize=(14, 9))

    # --------------------------------------------------------
    # Plot each encoding separately
    # --------------------------------------------------------

    for encoding, encoding_df in target_df.groupby(
        "encoding",
        sort=False
    ):

        # ----------------------------------------------------
        # Remove rows without an error tolerance
        # ----------------------------------------------------

        encoding_df = encoding_df[
            encoding_df["error_tolerance"].notna()
        ].copy()

        if encoding_df.empty:
            continue

        # ----------------------------------------------------
        # X / Y values
        # ----------------------------------------------------

        x_values = encoding_df[
            "dna_cost_nt_per_bit"
        ].tolist()

        y_values = encoding_df[
            "error_tolerance"
        ].tolist()

        # ----------------------------------------------------
        # Labels
        # ----------------------------------------------------

        labels = [
            make_configuration_label(
                encoding,
                configuration
            )
            for configuration in encoding_df[
                "configuration"
            ]
        ]

        # ----------------------------------------------------
        # Colour
        # ----------------------------------------------------

        encoding_color = encoding_colors.get(
            encoding,
            "black"
        )

        # ----------------------------------------------------
        # Scatter
        # ----------------------------------------------------

        plt.scatter(
            x_values,
            y_values,
            s=70,
            alpha=0.8,
            color=encoding_color,
            label=encoding
        )

        # ----------------------------------------------------
        # Label offsets
        # ----------------------------------------------------

        label_offsets = [
            (7, 7),
            (7, -12),
            (-7, 7),
            (-7, -12),
            (12, 0),
            (-12, 0),
            (0, 12),
            (0, -15),
        ]

        # ----------------------------------------------------
        # Add labels
        # ----------------------------------------------------

        for i, (
            x_value,
            y_value,
            label
        ) in enumerate(
            zip(
                x_values,
                y_values,
                labels
            )
        ):

            offset = label_offsets[
                i % len(label_offsets)
            ]

            plt.annotate(
                label,
                (x_value, y_value),
                xytext=offset,
                textcoords="offset points",
                fontsize=8,
                ha="center",
                va="center",
                color=encoding_color,
                bbox=dict(
                    boxstyle="round,pad=0.2",
                    facecolor="white",
                    edgecolor="none",
                    alpha=0.8
                )
            )

    # ========================================================
    # AXES / TITLE
    # ========================================================

    recovery_percent = int(
        target_recovery * 100
    )

    plt.xlabel(
        "DNA cost (nt/bit)"
    )

    plt.ylabel(
        "Tolerated IID substitution error rate"
    )

    plt.title(
        f"DNA storage cost vs error tolerance "
        f"({recovery_percent}% recovery)",
        fontsize=16,
        pad=15
    )

    plt.xscale("symlog", linthresh=0.1)
    plt.yscale(
        "symlog",
        linthresh=0.0001,
        linscale=1
    )


    # Normal decimal labels
    x_formatter = ScalarFormatter()
    x_formatter.set_scientific(False)
    x_formatter.set_useOffset(False)

    y_formatter = ScalarFormatter()
    y_formatter.set_scientific(False)
    y_formatter.set_useOffset(False)

    plt.gca().xaxis.set_major_formatter(x_formatter)
    plt.gca().yaxis.set_major_formatter(
        FuncFormatter(lambda x, pos: f"{x:.4f}")
    )

    # Grid including zero
    plt.grid(
        True,
        which="both",
        alpha=0.3
    )

    plt.legend(
        title="Encoding"
    )

    plt.tight_layout(
        rect=[0, 0, 1, 0.95]
    )

    # ========================================================
    # CREATE OUTPUT DIRECTORY
    # ========================================================

    os.makedirs(
        output_dir,
        exist_ok=True
    )

    # ========================================================
    # SAVE
    # ========================================================

    png_file = os.path.join(
        output_dir,
        f"cost_vs_error_tolerance_"
        f"{recovery_percent}pct.png"
    )

    pdf_file = os.path.join(
        output_dir,
        f"cost_vs_error_tolerance_"
        f"{recovery_percent}pct.pdf"
    )

    plt.savefig(
        png_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.savefig(
        pdf_file,
        bbox_inches="tight"
    )

    plt.savefig(png_file, dpi=300, bbox_inches="tight")
    plt.savefig(pdf_file, bbox_inches="tight")

    plt.show()
    plt.close()


# ============================================================
# USAGE
# ============================================================

csv_file = "simulation_results_costs/cost_vs_error_tolerance.csv"

output_dir = "plots"

# 100% recovery
plot_cost_vs_tolerance_csv(
    csv_file,
    output_dir,
    1.0
)

# 90% recovery
plot_cost_vs_tolerance_csv(
    csv_file,
    output_dir,
    0.90
)
