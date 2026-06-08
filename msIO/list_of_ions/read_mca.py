import json
from dataclasses import dataclass, field
from typing import Any, Iterable, Literal


@dataclass
class Score:
    name: str
    value: float
    unit: str | None = None

    @staticmethod
    def from_dict(d: dict[str, Any]) -> 'Score':
        return Score(
            name=d.get("name"),
            value=float(d.get("value")) if d.get("value") is not None else None,
            unit=d.get("unit")
        )


@dataclass
class Source:
    tool_name: str | None = None
    db_name: str | None = None

    @staticmethod
    def from_dict(d: dict[str, Any] | None) -> 'Source':
        if d is None:
            return Source()
        return Source(
            tool_name=d.get("toolName"),
            db_name=d.get("dbName")
        )


@dataclass
class Annotation:
    name: str
    formula: str | None = None
    scores: list[Score] = field(default_factory=list)
    source: Source | None = None

    @staticmethod
    def from_dict(d: dict[str, Any]) -> 'Annotation':
        return Annotation(
            name=d.get("name"),
            formula=d.get("formula"),
            scores=[Score.from_dict(s) for s in d.get("scores", [])],
            source=Source.from_dict(d.get("source"))
        )


@dataclass
class Ion:
    mz: float
    rt_seconds: float
    mobility: float | None = None
    ccs: float | None = None
    notation: str | None = None

    @staticmethod
    def from_dict(d: dict[str, Any]) -> 'Ion':
        return Ion(
            mz=float(d.get("mz")),
            rt_seconds=float(d.get("rt")),  # JSON uses "rt", we store as seconds
            mobility=float(d.get("mobility")) if d.get("mobility") is not None else None,
            ccs=float(d.get("ccs")) if d.get("ccs") is not None else None,
            notation=d.get("notation")
        )


@dataclass
class MoleculeAnnotation:
    neutral_mass: float | None = None
    ions: list[Ion] = field(default_factory=list)
    annotations: list[Annotation] = field(default_factory=list)
    flag_names: list[str] = field(default_factory=list)

    @staticmethod
    def from_dict(d: dict[str, Any]) -> 'MoleculeAnnotation':
        return MoleculeAnnotation(
            neutral_mass=float(d.get("neutralMass")) if d.get("neutralMass") is not None else None,
            ions=[Ion.from_dict(i) for i in d.get("ions", [])],
            annotations=[Annotation.from_dict(a) for a in d.get("annotations", [])],
            flag_names=d.get("flagNames", [])
        )


@dataclass
class AnnotationMatch:
    molecule: MoleculeAnnotation
    ion: Ion
    annotations: list[Annotation]
    deltas: dict[str, float]  # keys: "mz_da", "mz_ppm", "rt_s", "ccs_A"


class McaImportManager:
    model_version: str | None = None
    title: str | None = None
    export_date: str | None = None
    molecule_annotations: list[MoleculeAnnotation] = []

    def __init__(self, path_mca: str):
        self.path_mca = path_mca
        self.model_version = None
        self.title = None
        self.export_date = None
        self.molecule_annotations = []

        # Build a flat index of ions for fast scanning
        self._ion_index: list[dict[str, Any]] = []

        self._load()

    def _load(self):
        with open(self.path_mca, "r", encoding="utf-8") as f:
            data = json.load(f)

        self.model_version = data.get("modelVersion")
        self.title = data.get("title")
        self.export_date = data.get("exportDate")
        self.molecule_annotations = [
            MoleculeAnnotation.from_dict(ma) for ma in data.get("moleculeAnnotations", [])
        ]

        # Build a flat index (list of dicts pointing to MoleculeAnnotation and Ion)
        for mol in self.molecule_annotations:
            for ion in mol.ions:
                self._ion_index.append({
                    "molecule": mol,
                    "ion": ion
                })

    @staticmethod
    def _ppm(delta_da: float, ref_mz: float) -> float:
        if ref_mz == 0:
            return float("inf")
        return abs(delta_da) * 1e6 / ref_mz

    def find_annotations(
        self,
        mz: float,
        rt_seconds: float,
        ccs: float,
        mz_tol_ppm: float | None = 10.0,
        mz_tol_da: float | None = None,
        rt_tol_seconds: float = 5.0,
        ccs_tol_A: float = 1.0,
        require_all_within_tolerance: bool = True,
        sort_by: Literal["ppm", "rt", "ccs", "combined"] = "combined",
        require_unique_match: bool = True
    ) -> list[AnnotationMatch] | AnnotationMatch:
        """
        Find annotations whose associated ion matches the given mz/rt/ccs
        within specified tolerances.

        Parameters:
        - mz_tol_ppm: mass tolerance in ppm (used if mz_tol_da is None)
        - mz_tol_da: mass tolerance in Da (overrides mz_tol_ppm if provided)
        - rt_tol_seconds: retention time tolerance in seconds
        - ccs_tol_A: CCS tolerance in Å
        - require_all_within_tolerance: if True, all three dimensions must be
          within tolerance
        - sort_by: how to sort matches; "combined" sorts by a simple sum of
          normalized deltas

        Returns:
        - A list of AnnotationMatch objects.
        """
        matches: list[AnnotationMatch] = []

        for rec in self._ion_index:
            ion: Ion = rec["ion"]
            mol: MoleculeAnnotation = rec["molecule"]

            if ion.mz is None or ion.rt_seconds is None or ion.ccs is None:
                continue

            dmz_da = ion.mz - mz
            dmz_abs_da = abs(dmz_da)
            dmz_ppm = self._ppm(dmz_da, mz)

            drt_s = abs(ion.rt_seconds - rt_seconds)
            dccs_A = abs(ion.ccs - ccs)

            # Check tolerances
            if mz_tol_da is not None:
                in_mz = dmz_abs_da <= mz_tol_da
            else:
                in_mz = dmz_ppm <= (mz_tol_ppm if mz_tol_ppm is not None else 10.0)

            in_rt = drt_s <= rt_tol_seconds
            in_ccs = dccs_A <= ccs_tol_A

            if require_all_within_tolerance:
                if not (in_mz and in_rt and in_ccs):
                    continue
            else:
                if not (in_mz or in_rt or in_ccs):
                    continue

            matches.append(AnnotationMatch(
                molecule=mol,
                ion=ion,
                annotations=mol.annotations,
                deltas={
                    "mz_da": dmz_da,
                    "mz_ppm": dmz_ppm,
                    "rt_s": drt_s,
                    "ccs_A": dccs_A
                }
            ))

        if require_unique_match:
            assert len(matches) == 1, \
                "Multiple matches found, even though one is required. Consider lowering the tolerances."
            return matches[0]

        # Sorting
        if sort_by == "ppm":
            matches.sort(key=lambda m: abs(m.deltas["mz_ppm"]))
        elif sort_by == "rt":
            matches.sort(key=lambda m: abs(m.deltas["rt_s"]))
        elif sort_by == "ccs":
            matches.sort(key=lambda m: abs(m.deltas["ccs_A"]))
        else:
            # combined: simple weighted sum, gives priority to ppm then RT then CCS
            def combined_key(m: AnnotationMatch) -> float:
                # Normalize ppm and rt, ccs loosely
                return abs(m.deltas["mz_ppm"]) + (m.deltas["rt_s"] / max(rt_tol_seconds, 1e-9)) + (m.deltas["ccs_A"] / max(ccs_tol_A, 1e-9))
            matches.sort(key=combined_key)

        return matches


if __name__ == "__main__":
    file = r"\\hlabstorage.dmz.marum.de\scratch\Yannick\Guaymas new method height recursive\MetaboScape\annotations.mca"
    mgr = McaImportManager(file)

    # multiple annotations?
    multi_ann = []
    for mol_ann in mgr.molecule_annotations:
        if len([ann for ann in mol_ann.annotations if ann.source.tool_name == 'SmartFormula']) > 1:
            # should both be spectral library annotations
            multi_ann.append(mol_ann)

    # results = mgr.find_annotations(mz=1319.34919, rt_seconds=62.15 * 60, ccs=405, mz_tol_ppm=5, rt_tol_seconds=0.02 * 60, ccs_tol_A=.2)
    results = mgr.find_annotations(mz=330.24252, rt_seconds=1.98 * 60, ccs=178.2, mz_tol_ppm=5, rt_tol_seconds=0.02 * 60, ccs_tol_A=.2)
    for r in results:
        print(f"Match: Ion m/z={r.ion.mz:.6f}, rt={r.ion.rt_seconds:.2f}s, ccs={r.ion.ccs:.2f} Å, deltas={r.deltas}")
        for a in r.annotations:
            print(f"  Annotation: {a.name} ({a.formula}), source={a.source.tool_name if a.source else None}, scores={[s.name for s in a.scores]}")
