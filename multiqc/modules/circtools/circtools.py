import logging
import numpy as np
from typing import Dict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)


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

        data_by_sample: Dict[str, Dict[str, float]] = {}

        for f in self.find_log_files("circtools/detect", filehandles=True):
            parsed_samples = parse_circrnacount(f)
            if not parsed_samples:
                continue

            for s_name, metrics in parsed_samples.items():
                if s_name in data_by_sample:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")

                data_by_sample[s_name] = metrics
                self.add_data_source(f, section="CircRNACount")
                self.add_software_version(None, sample=s_name)

        data_by_sample = self.ignore_samples(data_by_sample)
        if not data_by_sample:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(data_by_sample)} circtools samples")

        self.write_data_file(data_by_sample, "multiqc_circtools")
        self.stats_tables(data_by_sample)

        self.add_section(
            name="circRNA Detection",
            anchor="circtools_detection",
            plot=circtools_detection_plot(data_by_sample),
        )

    # ------------------------------------------------------------------
    # Tables
    # ------------------------------------------------------------------

    def stats_tables(self, data_by_sample: Dict[str, Dict[str, float]]) -> None:
        headers = {
            "num_detected_circRNAs": {
                "namespace": "circtools",
                "title": "circRNAs",
                "description": "Detected circRNAs (>0 BSJ reads)",
                "scale": "Blues",
                "hidden": False,
            },
            "total_circRNA_reads": {
                "namespace": "circtools",
                "title": "circRNA reads",
                "description": "Total backsplice junction reads (BSJ)",
                "scale": "PuRd",
                "format": "{:,.0f}",
                "hidden": False,
            },
            "mean_circRNA_reads": {
                "namespace": "circtools",
                "title": "Mean BSJ",
                "description": "Mean BSJ reads per detected circRNA",
                "format": "{:.2f}",
                "scale": "OrRd",
                "hidden": True,
            },
            "median_circRNA_reads": {
                "namespace": "circtools",
                "title": "Median BSJ",
                "description": "Median BSJ reads per detected circRNA",
                "format": "{:.1f}",
                "scale": "OrRd",
                "hidden": True,
            },
            "max_circRNA_reads": {
                "namespace": "circtools",
                "title": "Max BSJ",
                "description": "Maximum BSJ reads for a single circRNA",
                "scale": "Reds",
                "hidden": True,
            },
        }

        self.general_stats_addcols(data_by_sample, headers, namespace="circtools")

        self.add_section(
            name="Summary Statistics",
            anchor="circtools_summary",
            description="Summary statistics from circtools detect (CircRNACount).",
            plot=table.plot(
                data_by_sample,
                headers,
                pconfig={
                    "id": "circtools_summary_table",
                    "title": "circtools: Summary Statistics",
                    "namespace": "circtools",
                },
            ),
        )


# ------------------------------------------------------------------
# Parsers
# ------------------------------------------------------------------

def parse_circrnacount(f) -> Dict[str, Dict[str, float]]:
    header = None
    sample_names = []
    tmp: Dict[str, list] = {}

    for line in f["f"]:
        line = line.strip()
        if not line:
            continue

        cols = line.split("\t")

        # Header
        if header is None:
            header = cols
            sample_names = [
                s.replace(".Chimeric.out.junction", "")
                .replace("_Chimeric.out.junction", "")
                .replace(".Chimeric", "")
                .replace("_Chimeric", "")
                for s in header[4:]
            ]


            for s in sample_names:
                tmp[s] = []
            continue

        # Data
        for i, s in enumerate(sample_names):
            try:
                tmp[s].append(int(cols[4 + i]))
            except (ValueError, IndexError):
                continue

    out: Dict[str, Dict[str, float]] = {}
    for s, values in tmp.items():
        if not values:
            continue

        arr = np.array(values)
        out[s] = {
            "total_circRNA_reads": int(arr.sum()),
            "num_detected_circRNAs": int((arr > 0).sum()),
            "mean_circRNA_reads": float(arr.mean()),
            "median_circRNA_reads": float(np.median(arr)),
            "max_circRNA_reads": int(arr.max()),
        }

    return out


# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------

def circtools_detection_plot(data_by_sample):
    keys = {
        "num_detected_circRNAs": {
            "color": "#437bb1",
            "name": "Detected circRNAs"
        },
    }

    pconfig = {
        "id": "circtools_detection_plot",
        "title": "circtools: circRNA Detection",
        "ylab": "Count",
        "cpswitch_counts_label": "Number of circRNAs",
    }

    return bargraph.plot(data_by_sample, keys, pconfig)
