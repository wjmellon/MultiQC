import logging
import re
import numpy as np
from typing import Dict, Optional

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table, scatter

log = logging.getLogger(__name__)

# ------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------

def _clean_circtools_column_name(raw: str) -> str:
    return (
        raw.replace(".Chimeric.out.junction", "")
           .replace("_Chimeric.out.junction", "")
           .replace(".Chimeric", "")
           .replace("_Chimeric", "")
    )


# ------------------------------------------------------------------
# STAR parser (simplified, file-contents based)
# ------------------------------------------------------------------

def parse_star_report(contents: str) -> Optional[Dict[str, float]]:
    regexes = {
        "total_reads": r"Number of input reads \|\s+(\d+)",
        "uniquely_mapped": r"Uniquely mapped reads number \|\s+(\d+)",
        "multimapped": r"Number of reads mapped to multiple loci \|\s+(\d+)",
        "multimapped_toomany": r"Number of reads mapped to too many loci \|\s+(\d+)",
    }

    parsed: Dict[str, float] = {}
    for key, rgx in regexes.items():
        m = re.search(rgx, contents, re.MULTILINE)
        if m:
            parsed[key] = float(m.group(1))

    if not parsed:
        return None

    parsed["mapped"] = parsed.get("uniquely_mapped", 0) + parsed.get("multimapped", 0)
    return parsed


# ------------------------------------------------------------------
# MultiQC Module
# ------------------------------------------------------------------

class MultiqcModule(BaseMultiqcModule):
    """
    Metrics based on circtools quickcheck module
    """

    def __init__(self):
        super().__init__(
            name="circtools",
            anchor="circtools",
            href="https://github.com/jakobilab/circtools",
            info="Produced by circtools",
        )

        # --------------------------------------------------
        # Parse circRNA counts
        # --------------------------------------------------
        data_by_sample: Dict[str, Dict[str, float]] = {}

        for f in self.find_log_files("circtools/detect", filehandles=True):
            if f["fn"].endswith("CircRNACountName.txt"):
                continue

            parsed = self.parse_circrnacount(f)
            for s, metrics in parsed.items():
                data_by_sample[s] = metrics
                self.add_data_source(f, section="CircRNACount")

        # --------------------------------------------------
        # Parse linear counts
        # --------------------------------------------------
        linear_by_sample: Dict[str, int] = {}

        for f in self.find_log_files("circtools/detect/linear", filehandles=True):
            parsed = self.parse_linearcount(f)
            linear_by_sample.update(parsed)

        for s in data_by_sample:
            data_by_sample[s]["total_linear_reads"] = linear_by_sample.get(s, np.nan)

        data_by_sample = self.ignore_samples(data_by_sample)
        if not data_by_sample:
            raise ModuleNoSamplesFound

        # --------------------------------------------------
        # Parse STAR logs
        # --------------------------------------------------
        star_by_sample: Dict[str, Dict[str, float]] = {}

        for f in self.find_log_files("star"):
            parsed = parse_star_report(f["f"])
            if parsed:
                star_by_sample[f["s_name"]] = parsed

        # --------------------------------------------------
        # Derived metrics
        # --------------------------------------------------
        for s, d in data_by_sample.items():
            uniq = star_by_sample.get(s, {}).get("uniquely_mapped")
            if uniq and uniq > 0:
                d["circRNAs_per_million_unique"] = d["num_detected_circRNAs"] / (uniq / 1e6)
            else:
                d["circRNAs_per_million_unique"] = np.nan

        # --------------------------------------------------
        # Output
        # --------------------------------------------------
        self.write_data_file(data_by_sample, "multiqc_circtools")
        self.stats_tables(data_by_sample)

        self.add_section(
            name="circRNAs per Million Unique Reads",
            anchor="circtools_norm",
            plot=circtools_detection_plot(data_by_sample),
        )

        self.add_section(
            name="Circular vs Linear Reads",
            anchor="circtools_circ_vs_linear",
            plot=circ_vs_linear_plot(data_by_sample),
        )

        self.add_section(
            name="Unique Reads vs Detected circRNAs",
            anchor="circtools_unique_vs_circs",
            plot=unique_vs_circs_plot(data_by_sample, star_by_sample),
        )

    # ------------------------------------------------------------------
    # Tables
    # ------------------------------------------------------------------

    def stats_tables(self, data_by_sample: Dict[str, Dict[str, float]]) -> None:
        headers = {
            "num_detected_circRNAs": {"title": "circRNAs"},
            "total_circRNA_reads": {"title": "circRNA reads"},
            "total_linear_reads": {"title": "Linear reads"},
            "circRNAs_per_million_unique": {"title": "circRNAs / M uniq"},
        }

        self.general_stats_addcols(data_by_sample, headers, namespace="circtools")

        self.add_section(
            name="Summary Statistics",
            anchor="circtools_summary",
            plot=table.plot(
                data_by_sample,
                headers,
                pconfig={
                    "id": "circtools_summary_table",
                    "title": "circtools: Summary Statistics",
                },
            ),
        )

    # ------------------------------------------------------------------
    # Parsers (instance methods – REQUIRED)
    # ------------------------------------------------------------------

    def parse_linearcount(self, f) -> Dict[str, int]:
        header = None
        sample_names: list[str] = []
        tmp: Dict[str, list[int]] = {}

        for line in f["f"]:
            cols = line.rstrip().split("\t")
            if not cols:
                continue

            # Header
            if header is None:
                header = cols
                sample_names = [
                    self.clean_s_name(_clean_circtools_column_name(s), f)
                    for s in header[3:]   # ← linear starts after Chr/Start/End
                ]
                for s in sample_names:
                    tmp[s] = []
                continue

            # Data rows
            for i, s in enumerate(sample_names):
                try:
                    tmp[s].append(int(cols[3 + i]))
                except (ValueError, IndexError):
                    pass

        # Sum per sample
        return {s: int(np.sum(v)) for s, v in tmp.items() if v}




    def parse_circrnacount(self, f) -> Dict[str, Dict[str, float]]:
        header = None
        sample_names: list[str] = []
        tmp: Dict[str, list[int]] = {}

        for line in f["f"]:
            cols = line.rstrip().split("\t")
            if not cols:
                continue

            # Header
            if header is None:
                header = cols
                sample_names = [
                    self.clean_s_name(_clean_circtools_column_name(s), f)
                    for s in header[4:]
                ]
                for s in sample_names:
                    tmp[s] = []
                continue

            # Data
            for i, s in enumerate(sample_names):
                try:
                    tmp[s].append(int(cols[4 + i]))
                except Exception:
                    pass

        out: Dict[str, Dict[str, float]] = {}
        for s, vals in tmp.items():
            if not vals:
                continue

            arr = np.array(vals)
            out[s] = {
                "total_circRNA_reads": int(arr.sum()),
                "num_detected_circRNAs": int((arr > 0).sum()),
            }

        return out




# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------

def circtools_detection_plot(data_by_sample):
    keys = {
        "circRNAs_per_million_unique": {
            "color": "#3874c8",
            "name": "circRNAs per million uniquely mapped reads",
        }
    }

    pconfig = {
        "id": "circtools_norm_bar_v2",  
        "title": "Detected circRNAs per Million Unique Reads",
        "ylab": "circRNAs / million reads",
        "cpswitch_counts_label": "circRNAs per million unique reads",
    }

    return bargraph.plot(data_by_sample, keys, pconfig)




def circ_vs_linear_plot(data_by_sample):
    plot_data = {
        s: {"x": d["total_circRNA_reads"], "y": d["total_linear_reads"]}
        for s, d in data_by_sample.items()
        if d.get("total_circRNA_reads") is not None
        and d.get("total_linear_reads") is not None
    }

    return scatter.plot(
        plot_data,
        pconfig={
            "id": "circtools_circ_vs_linear_scatter_v2",  
            "title": "Circular vs Linear Read Counts",
            "xlab": "Circular reads (BSJ)",
            "ylab": "Linear reads",
            "xlog": True,
            "ylog": True,
        },
    )



def unique_vs_circs_plot(data_by_sample, star_by_sample):
    plot_data = {
        s: {"x": star_by_sample[s]["uniquely_mapped"], "y": d["num_detected_circRNAs"]}
        for s, d in data_by_sample.items()
        if s in star_by_sample and "uniquely_mapped" in star_by_sample[s]
    }

    return scatter.plot(
        plot_data,
        pconfig={
            "id": "circtools_unique_vs_circs_scatter_v2", 
            "title": "Unique Reads vs Detected circRNAs",
            "xlab": "Uniquely mapped reads",
            "ylab": "Detected circRNAs",
        },
    )


