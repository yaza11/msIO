"""
Overview figure of measured features for the annotation GUI.

Three modes:
    "match"       a measured feature and its library match (library MS2 as mirror)
    "feature"     a single measured feature
    "comparison"  two measured features (second MS2 as mirror)

Layout: MS2 spectra on top, below the MS1 spectrum (measured and theoretical isotope pattern, peaks
labelled with their distance to the monoisotopic peak) and the intensities in the samples, next to
them the chemical structure.

Gathering the data (database queries) and building the figure don't touch tk, so both can run in a
worker thread. Peak annotations of the MS2 spectra can be changed by the user, the figure returns the
plotted peaks so that clicks can be mapped to them.
"""

import os
from dataclasses import dataclass, field

import numpy as np
from LipidCalculator import CompoundDict
from LipidCalculator.isotopes.isotope_pattern import IsotopePattern
from LipidCalculator.rdkit.plotting import mplt_mol
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from rdkit import Chem

from msIO import PeakList
from msIO.feature_managers.db import FeatureManagerDB, Library
from msIO.metrics import cosine_similarity_sym

# identifies a spectrum for user defined peak annotations: ("meas", feature id) or ("lib", "library:feature id"),
# feature ids of different libraries overlap (see library_spectrum_key)
SpectrumKey = tuple[str, int | str]

# peaks within this distance of the precursor m/z are hidden by the isotope filter
ISOTOPE_WINDOW_DA: float = 4.
# measured MS1 peaks are shown from precursor m/z - MS1_WINDOW_DA[0] to + MS1_WINDOW_DA[1]
MS1_WINDOW_DA: tuple[float, float] = (1.5, 5.5)
# MS1 peaks below this relative intensity are not labelled with their distance to the monoisotopic peak
MS1_LABEL_CUTOFF: float = 0.01
# known distances to the monoisotopic peak, used to label MS1 peaks
_C13: float = (CompoundDict('C[13]') - CompoundDict('C')).mass
KNOWN_ISOTOPE_DISTANCES: dict[str, float] = {
    '13C': _C13,
    '2x13C': 2 * _C13,
    '34S': (CompoundDict('S[34]') - CompoundDict('S')).mass,
    '18O': (CompoundDict('O[18]') - CompoundDict('O')).mass,
    '13C, z=2': _C13 / 2,
}
KNOWN_DISTANCE_TOLERANCE_DA: float = 2e-3
# labels of the known distances (mathtext for the superscripts)
KNOWN_ISOTOPE_LABELS: dict[str, str] = {
    '13C': r'$^{13}$C',
    '2x13C': r'2×$^{13}$C',
    '34S': r'$^{34}$S',
    '18O': r'$^{18}$O',
    '13C, z=2': r'$^{13}$C, z=2',
}
# height (relative intensity) of the lowest distance arrow in the MS1 plot and the step between stacked arrows
MS1_ARROW_BASE: float = 1.1
MS1_ARROW_STEP: float = 0.17

FEATURE_COLORS: tuple[str, str] = ("C2", "C1")  # same as the compared features in the scatter plot
LIBRARY_COLOR: str = "C3"
THEORETICAL_COLOR: str = "C7"
# plots of the overview, OverviewSettings.panels selects the shown ones
OVERVIEW_PANELS: tuple[str, ...] = ("ms2", "ms1", "intensities", "structure")


@dataclass
class FeatureData:
    """Properties and spectra of a measured feature."""
    f_id: int
    mz: float
    rt_min: float | None
    ccs: float | None
    ms1: PeakList | None
    ms2: PeakList | None
    intensities: dict[str, float]
    # manual annotation of the feature, if any (name, formula, adduct, smiles)
    annotation: dict | None = None


@dataclass
class LibraryData:
    f_id: int
    name: str | None
    formula: str | None
    adduct: str | None
    smiles: str | None
    mz: float | None
    ccs: float | None
    ms2: PeakList | None
    library: str = ""  # file name of the library

    @property
    def key(self) -> SpectrumKey:
        return library_spectrum_key(self.library, self.f_id)


def library_spectrum_key(library: str, f_id: int) -> SpectrumKey:
    """Key of the MS2 spectrum of a library feature."""
    return "lib", f"{library}:{f_id}" if library else f_id


@dataclass
class OverviewData:
    mode: str
    features: list[FeatureData]
    library: LibraryData | None = None
    match_result: dict | None = None


@dataclass
class OverviewSettings:
    hide_isotope_peaks: bool = False
    annotation_cutoff: float = 0.3
    ms2_tolerance_da: float = 0.01
    # mass accuracy of the theoretical isotope pattern
    ms1_tolerance_da: float = 0.003
    peak_annotations: dict[SpectrumKey, dict[float, str]] = field(default_factory=dict)
    # plots to draw, see OVERVIEW_PANELS
    panels: tuple[str, ...] = ("ms2", "ms1", "intensities", "structure")


@dataclass
class PlottedSpectrum:
    """Peaks of a spectrum as plotted (signed relative intensities), to find the peak under a click."""
    key: SpectrumKey
    ax: Axes
    mzs: np.ndarray
    heights: np.ndarray


def peak_key(mz: float) -> float:
    """Peaks are identified by their m/z rounded to 4 decimals for user defined annotations."""
    return round(float(mz), 4)


def gather_feature(
        meas_db: FeatureManagerDB,
        f_id: int,
        ms2_spectra: dict | None = None,
        annotation: dict | None = None
) -> FeatureData:
    rt = meas_db.retention_times_in_seconds.get(f_id)
    ccs = meas_db.collisional_cross_sections.get(f_id)
    ms2 = ms2_spectra.get(f_id) if ms2_spectra is not None else meas_db.get_ms_spectrum(f_id, level=2)
    return FeatureData(
        f_id=f_id,
        mz=meas_db.mzs[f_id],
        rt_min=None if rt is None else rt / 60,
        ccs=ccs if ccs is not None and ccs > 0 else None,
        ms1=meas_db.get_ms_spectrum(f_id, level=1),
        ms2=ms2,
        intensities=meas_db.get_intensities(f_id),
        annotation=annotation,
    )


def gather_library(lib_db: Library, f_id: int) -> LibraryData:
    props: dict = lib_db.get_values_for_feature(f_id)
    ccs = props.get('CCS')
    return LibraryData(
        f_id=f_id,
        name=props.get('name_sirius'),
        formula=props.get('formula_metaboscape'),
        adduct=props.get('adduct_metaboscape'),
        smiles=props.get('smiles'),
        mz=props.get('mz_meas'),
        ccs=ccs if ccs is not None and ccs > 0 else None,
        ms2=lib_db.get_ms_spectrum(f_id, level=2),
        library=os.path.basename(lib_db.path_file),
    )


def _filter_isotope_peaks(mzs: np.ndarray, intensities: np.ndarray, annotations: np.ndarray,
                          precursor_mz: float | None, hide: bool):
    if not hide or precursor_mz is None:
        return mzs, intensities, annotations
    keep = np.abs(mzs - precursor_mz) > ISOTOPE_WINDOW_DA
    return mzs[keep], intensities[keep], annotations[keep]


def _plot_ms2(
        ax: Axes,
        spectrum: PeakList | None,
        key: SpectrumKey,
        precursor_mz: float | None,
        color: str,
        as_mirror: bool,
        settings: OverviewSettings
) -> PlottedSpectrum | None:
    """
    Stem plot normalized to 1 (negative for mirror). Peaks annotated by the user (stored in a measured spectrum or set
    in this session, settings.peak_annotations) are always labelled. Other peaks above the cutoff are labelled with
    the annotation of the library spectrum (e.g. "-H2O", every peak of some libraries has one) or their m/z.
    """
    if spectrum is None or len(spectrum.mzs) == 0:
        ax.text(0.5, 0.25 if as_mirror else 0.75, f"no MS2 spectrum for {'library ' if key[0] == 'lib' else ''}"
                f"feature {key[1]}", transform=ax.transAxes, ha="center", color=color)
        return None
    mzs, intensities, stored_annotations = _filter_isotope_peaks(
        np.asarray(spectrum.mzs, dtype=float), np.asarray(spectrum.intensities, dtype=float),
        np.asarray(spectrum.annotations, dtype=object), precursor_mz, settings.hide_isotope_peaks
    )
    if len(mzs) == 0:
        return None
    heights = intensities / intensities.max() * (-1 if as_mirror else 1)
    ax.vlines(mzs, 0, heights, colors=color, linewidth=1)

    user_annotations: dict[float, str] = settings.peak_annotations.get(key, {})
    # annotations stored in measured spectra are set by the user (see AnnotationGUI.on_click_overview)
    stored_by_user: bool = key[0] == "meas"
    for mz, height, stored in zip(mzs, heights, stored_annotations):
        label = user_annotations.get(peak_key(mz)) or (stored if stored_by_user else None)
        if not label:
            if abs(height) < settings.annotation_cutoff:
                continue
            label = stored or f"{mz:.4f}"
        ax.annotate(label, xy=(mz, height), ha="center", va="top" if as_mirror else "bottom", rotation=90,
                    fontsize=7, color="dimgray" if label == f"{mz:.4f}" else "black")
    return PlottedSpectrum(key=key, ax=ax, mzs=mzs, heights=heights)


def _plot_ms1(ax: Axes, data: OverviewData, settings: OverviewSettings):
    """
    Measured MS1 spectra (and the theoretical isotope pattern if formula and adduct are known), all pointing up and
    normalized to 1. Measured peaks are labelled with their distance to the monoisotopic peak. Two compared features
    are shifted to their precursor m/z.
    """
    relative: bool = data.mode == "comparison"
    formula, adduct = _get_formula_adduct(data)
    handles: list = []
    # distance arrows above the spectra, two levels per feature: neighboring peaks and every second peak
    level: int = 0
    for feature, color in zip(data.features, FEATURE_COLORS):
        if feature.ms1 is None or len(feature.ms1.mzs) == 0:
            continue
        mzs, intensities = np.asarray(feature.ms1.mzs, dtype=float), np.asarray(feature.ms1.intensities, dtype=float)
        in_window = (mzs > feature.mz - MS1_WINDOW_DA[0]) & (mzs < feature.mz + MS1_WINDOW_DA[1])
        mzs, intensities = mzs[in_window], intensities[in_window]
        if len(mzs) == 0:
            continue
        heights = intensities / intensities.max()
        offset = feature.mz if relative else 0.
        ax.vlines(mzs - offset, 0, heights, colors=color, linewidth=2)
        handles.append(Line2D([], [], color=color, label=f"feature {feature.f_id}"))

        # labelled peaks sorted by m/z, the monoisotopic peak is the one closest to the precursor m/z
        order = np.argsort(mzs)
        peaks = [(mz, height) for mz, height in zip(mzs[order], heights[order]) if height >= MS1_LABEL_CUTOFF]
        if len(peaks) < 2:
            continue
        idx_mono: int = int(np.argmin([abs(mz - feature.mz) for mz, _ in peaks]))
        y_neighbors: float = MS1_ARROW_BASE + level * MS1_ARROW_STEP
        y_second: float = y_neighbors + MS1_ARROW_STEP
        level += 2
        # lower level: neighboring peaks, upper level: every second peak counted from the monoisotopic one
        # (M -> M+2, M+2 -> M+4, ...), these arrows don't overlap
        pairs = [(i, i + 1, y_neighbors) for i in range(len(peaks) - 1)]
        pairs += [(i, i + 2, y_second) for i in range(idx_mono % 2, len(peaks) - 2, 2)]
        # dotted guides from the peaks up to their highest arrow
        y_top = {i: y for _, i, y in sorted(pairs, key=lambda pair: pair[2])}
        y_top |= {i: max(y, y_top.get(i, 0)) for i, _, y in pairs}
        ax.vlines([peaks[i][0] - offset for i in y_top], [peaks[i][1] for i in y_top], list(y_top.values()),
                  colors=color, linewidth=0.5, linestyles="dotted")
        for i, j, y in pairs:
            (mz_a, _), (mz_b, _) = peaks[i], peaks[j]
            distance = mz_b - mz_a
            known = [KNOWN_ISOTOPE_LABELS[name] for name, d in KNOWN_ISOTOPE_DISTANCES.items()
                     if abs(distance - d) <= KNOWN_DISTANCE_TOLERANCE_DA]
            ax.annotate("", xy=(mz_b - offset, y), xytext=(mz_a - offset, y),
                        arrowprops=dict(arrowstyle="<->", color=color, linewidth=0.8, shrinkA=0, shrinkB=0))
            ax.text((mz_a + mz_b) / 2 - offset, y, f"{distance:.4f}" + (f" ({', '.join(known)})" if known else ""),
                    ha="center", va="bottom", fontsize=7, color=color)

    window = (data.features[0].mz - MS1_WINDOW_DA[0], data.features[0].mz + MS1_WINDOW_DA[1])
    if formula and adduct and not relative:
        try:
            pattern = IsotopePattern.from_formula(formula=formula, adduct=adduct, mass_accuracy=settings.ms1_tolerance_da,
                                                  merge_method='weighted_average')
            masses, intensities = np.asarray(pattern.masses), np.asarray(pattern.intensities)
            # short horizontal markers at the theoretical intensity, the measured peaks stay visible below
            ax.plot(masses, intensities / intensities.max(), linestyle="none", marker="_", markersize=18,
                    markeredgewidth=2.5, color=THEORETICAL_COLOR, zorder=3)
            handles.append(Line2D([], [], linestyle="none", marker="_", markersize=12, markeredgewidth=2.5,
                                  color=THEORETICAL_COLOR, label=f"theoretical {formula} {adduct}"))
            if not window[0] < masses[np.argmax(intensities)] < window[1]:
                ax.text(0.02, 0.95, f"theoretical pattern at m/z {masses[np.argmax(intensities)]:.4f} (outside)",
                        transform=ax.transAxes, fontsize=7, va="top")
        except Exception as error:
            ax.text(0.02, 0.95, f"no theoretical pattern: {error}", transform=ax.transAxes, fontsize=7, va="top")
    if not relative:  # keep the measured peaks in view
        ax.set_xlim(*window)

    if not handles:
        ax.text(0.5, 0.5, "no MS1 spectrum", transform=ax.transAxes, ha="center")
    else:
        # the distance arrows are at the top, the higher isotope peaks on the right are small
        ax.legend(handles=handles, loc="center right", fontsize=7)
    ax.set_ylim(0, MS1_ARROW_BASE + max(level, 1) * MS1_ARROW_STEP + 0.1)  # room for the distance arrows
    ax.set_xlabel("m/z - precursor m/z" if relative else "m/z")
    ax.set_ylabel("Relative intensity")
    ax.set_title("MS1", fontsize=9)


def _plot_intensities(ax: Axes, data: OverviewData):
    """Grouped bar plot of the intensities of the measured features in each sample."""
    samples: list[str] = sorted({sample for feature in data.features for sample in feature.intensities})
    if not samples:
        ax.text(0.5, 0.5, "no intensities", transform=ax.transAxes, ha="center")
        return
    width: float = 0.8 / len(data.features)
    positions = np.arange(len(samples))
    for i, (feature, color) in enumerate(zip(data.features, FEATURE_COLORS)):
        values = [feature.intensities.get(sample) or 0 for sample in samples]
        ax.bar(positions + (i - (len(data.features) - 1) / 2) * width, values, width=width, color=color,
               label=f"feature {feature.f_id}")
    ax.set_xticks(positions, [s if len(s) <= 25 else s[:22] + "..." for s in samples], rotation=90, fontsize=6)
    ax.set_ylabel("Intensity")
    ax.set_title("Intensities in samples", fontsize=9)
    if len(data.features) > 1:
        ax.legend(fontsize=7)


def _get_formula_adduct(data: OverviewData) -> tuple[str | None, str | None]:
    """Formula and adduct for the theoretical isotope pattern: library match, else manual annotation."""
    if data.library is not None:
        return data.library.formula, data.library.adduct
    annotation = data.features[0].annotation if data.features else None
    if annotation:
        return annotation.get('formula'), annotation.get('adduct')
    return None, None


def _plot_structure(ax: Axes, data: OverviewData):
    """Structure of the library match or of the manual annotation of the (first) feature, with a short summary."""
    ax.axis("off")
    if data.library is not None:
        smiles, name = data.library.smiles, data.library.name
    elif data.features and data.features[0].annotation:
        smiles, name = data.features[0].annotation.get('smiles'), data.features[0].annotation.get('name')
    else:
        smiles, name = None, None

    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is not None:
        mplt_mol(mol=mol, ax=ax)
    else:
        ax.text(0.5, 0.6, "no structure", transform=ax.transAxes, ha="center", color="gray")

    lines: list[str] = []
    if name:
        lines.append(name)
    for feature in data.features:
        parts = [f"feature {feature.f_id}: m/z {feature.mz:.4f}"]
        if feature.rt_min is not None:
            parts.append(f"RT {feature.rt_min:.2f} min")
        if feature.ccs is not None:
            parts.append(f"CCS {feature.ccs:.1f} Å²")
        lines.append(", ".join(parts))
        if feature.annotation and data.mode != "match":
            lines.append(f"  annotated: {feature.annotation.get('name')} {feature.annotation.get('adduct') or ''}")
    if (lib := data.library) is not None:
        lines.append(f"library {lib.f_id} ({lib.library}): {lib.formula or ''} {lib.adduct or ''}"
                     + (f", m/z {lib.mz:.4f}" if lib.mz is not None else ""))
        if lib.ccs is not None and data.features[0].ccs is not None:
            lines.append(f"CCS lib. {lib.ccs:.1f} Å² (Δ {(data.features[0].ccs - lib.ccs) / lib.ccs * 100:+.1f} %)")
    if data.mode == "comparison" and len(data.features) > 1:
        lines.append(_describe_differences(*data.features[:2]))
    ax.text(0.5, 0.0, "\n".join(lines), transform=ax.transAxes, ha="center", va="top", fontsize=8, wrap=True)


def _describe_differences(first: FeatureData, second: FeatureData) -> str:
    """
    Differences of two compared features (second - first): m/z (Da, ppm), RT (min) and CCS (Å², %), one per line
    (the summary column is narrow).
    """
    dmz = second.mz - first.mz
    lines = [f"feature {second.f_id} - {first.f_id}:", f"Δm/z {dmz:+.4f} Da ({dmz / first.mz * 1e6:+.1f} ppm)"]
    if first.rt_min is not None and second.rt_min is not None:
        lines.append(f"ΔRT {second.rt_min - first.rt_min:+.2f} min")
    if first.ccs is not None and second.ccs is not None:
        dccs = second.ccs - first.ccs
        lines.append(f"ΔCCS {dccs:+.1f} Å² ({dccs / first.ccs * 100:+.1f} %)")
    return "\n".join(lines)


def _plot_ms2_panel(ax: Axes, data: OverviewData, settings: OverviewSettings) -> list[PlottedSpectrum]:
    """MS2: first feature on top, library match or second feature as mirror. Returns the plotted spectra."""
    first = data.features[0]
    spectra: list[PlottedSpectrum | None] = [
        _plot_ms2(ax, first.ms2, ("meas", first.f_id), first.mz, FEATURE_COLORS[0], False, settings)
    ]
    legend: list[Line2D] = [Line2D([], [], color=FEATURE_COLORS[0], label=f"feature {first.f_id}")]
    mirror: PeakList | None = None
    if data.mode == "match" and data.library is not None:
        lib = data.library
        mirror = lib.ms2
        spectra.append(_plot_ms2(ax, lib.ms2, lib.key, lib.mz, LIBRARY_COLOR, True, settings))
        legend.append(Line2D([], [], color=LIBRARY_COLOR, label=f"library {lib.f_id}: {lib.name}"))
        result = data.match_result or {}
        score = result.get('ms2_score')
        ax.set_title(f"MS2 score: {score:.3f} ({result.get('n_hits_ms2')} peaks matched)"
                     if score is not None and not np.isnan(score) else "MS2 (no score)", fontsize=10)
    elif data.mode == "comparison" and len(data.features) > 1:
        second = data.features[1]
        mirror = second.ms2
        spectra.append(_plot_ms2(ax, second.ms2, ("meas", second.f_id), second.mz, FEATURE_COLORS[1], True,
                                 settings))
        legend.append(Line2D([], [], color=FEATURE_COLORS[1], label=f"feature {second.f_id}"))
        if first.ms2 is not None and second.ms2 is not None:
            score, n_hits = cosine_similarity_sym(first.ms2, second.ms2, settings.ms2_tolerance_da, return_nhits=True)
            ax.set_title(f"Symmetric cosine similarity: {score:.3f} ({n_hits} peaks matched)", fontsize=10)
    else:
        ax.set_title(f"MS2 of feature {first.f_id}", fontsize=10)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.set_ylim(-1.35 if mirror is not None else 0, 1.35)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Relative intensity")
    ax.legend(handles=legend, loc="upper left", fontsize=8)
    return [spectrum for spectrum in spectra if spectrum is not None]


def build_overview_figure(fig: Figure, data: OverviewData, settings: OverviewSettings) -> list[PlottedSpectrum]:
    """
    Draw the overview into fig, returns the plotted MS2 spectra (for clickable peaks).

    Only the plots in settings.panels are drawn: MS2 on top over the full width, below MS1 and intensities on the
    left and the structure on the right (each of them takes the space of hidden ones).
    """
    shown: list[str] = [panel for panel in OVERVIEW_PANELS if panel in settings.panels]
    if not shown:
        fig.text(0.5, 0.5, "all plots are hidden", ha="center", va="center", color="gray")
        return []
    with_ms2: bool = "ms2" in shown
    left: list[str] = [panel for panel in ("ms1", "intensities") if panel in shown]
    with_structure: bool = "structure" in shown
    n_lower_rows: int = max(len(left), int(with_structure))
    n_columns: int = 2 if left and with_structure else 1
    # the bottom margin leaves room for the sample names of the intensity plot and the text below the structure
    gs = GridSpec(
        int(with_ms2) + n_lower_rows, n_columns, figure=fig,
        width_ratios=[2, 1] if n_columns == 2 else None,
        height_ratios=([1.3] if with_ms2 else []) + [1] * n_lower_rows,
        hspace=0.6, wspace=0.15, left=0.07, right=0.98, top=0.95,
        bottom=0.17 if "intensities" in shown or with_structure else 0.08
    )

    spectra: list[PlottedSpectrum] = []
    if with_ms2:
        spectra = _plot_ms2_panel(fig.add_subplot(gs[0, :]), data, settings)
    first_lower_row: int = int(with_ms2)
    for row, panel in enumerate(left, start=first_lower_row):
        ax = fig.add_subplot(gs[row, 0])
        if panel == "ms1":
            _plot_ms1(ax, data, settings)
        else:
            _plot_intensities(ax, data)
    if with_structure:
        _plot_structure(fig.add_subplot(gs[first_lower_row:, 1 if left else 0]), data)
    return spectra
