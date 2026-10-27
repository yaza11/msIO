"""
Manual annotations of measured features: the ones stored in the measured database and the changes collected in the
GUI (annotations to add, annotations to discard / remove), which are only written to the database on request.
"""

from msIO.feature_managers.db import FeatureManagerDB
from msIO.gui.annotation.common import ANNOTATION_KIND_LABELS, AnnotationKey


class AnnotationStore:
    """
    State of the manual annotations (no GUI).

    saved:     annotations in the measured database {key: info}, info holds kind, name, formula, adduct, smiles,
               ms2_score, reference (library or parent feature id), library_name (file name of the library) and the
               candidate_ids of the database rows
    pending:   annotations to add {key: info}
    discarded: pending annotations not to add and saved annotations to remove
    """

    def __init__(self):
        self.saved: dict[AnnotationKey, dict] = {}
        self.pending: dict[AnnotationKey, dict] = {}
        self.discarded: set[AnnotationKey] = set()

    @staticmethod
    def library_key(f_id_meas: int, library_name: str, f_id_lib: int) -> AnnotationKey:
        """Feature ids of different libraries overlap, so the reference is "library file name:library feature id"."""
        return "library", int(f_id_meas), f"{library_name}:{int(f_id_lib)}"

    def reset(self, meas_db: FeatureManagerDB | None):
        """Forget all changes and read the annotations of a (new) measured database."""
        self.pending.clear()
        self.discarded.clear()
        self.load_saved(meas_db)

    def load_saved(self, meas_db: FeatureManagerDB | None):
        """Read the manual annotations (all kinds) stored in the measured database."""
        self.saved = {}
        if meas_db is None:
            return
        for candidate in meas_db.get_manual_annotations():
            kind: str = FeatureManagerDB.get_manual_annotation_kind(candidate)
            reference: int | None = FeatureManagerDB.get_referenced_feature_id(candidate)
            library_name: str = FeatureManagerDB.get_referenced_library(candidate) or ""
            if kind == "library":
                key: AnnotationKey = self.library_key(candidate.feature_id, library_name, reference)
            elif kind == "manual":
                key = (kind, candidate.feature_id, f"{candidate.name_sirius}|{candidate.adduct_sirius}")
            else:
                key = (kind, candidate.feature_id, str(reference))
            if key in self.saved:  # stored more than once
                self.saved[key]["candidate_ids"].append(candidate.id)
                continue
            self.saved[key] = dict(
                kind=kind,
                name=candidate.name_sirius,
                formula=candidate.formula_sirius,
                adduct=candidate.adduct_sirius,
                smiles=candidate.smiles,
                ms2_score=candidate.confidence_score,
                reference=reference,
                library_name=library_name if kind == "library" else None,
                candidate_ids=[candidate.id],
            )

    def state(self, key: AnnotationKey) -> str:
        """State of an annotation: "in database", "to remove", "to add" or "" (none)."""
        if key in self.saved:
            return "to remove" if key in self.discarded else "in database"
        if key in self.pending and key not in self.discarded:
            return "to add"
        return ""

    def unsaved_changes(self) -> tuple[list[AnnotationKey], list[AnnotationKey]]:
        """Annotations to add to and to remove from the measured database."""
        to_add = [key for key in self.pending if key not in self.discarded]
        to_remove = [key for key in self.saved if key in self.discarded]
        return to_add, to_remove

    def add(self, key: AnnotationKey, info: dict) -> bool:
        """
        Add an annotation to the pending ones. Returns False if it is already in the database (it is kept, if it was
        marked for removal).
        """
        self.discarded.discard(key)
        if key in self.saved:
            return False
        self.pending[key] = info
        return True

    def toggle(self, keys: list[AnnotationKey], keep: bool = None):
        """Mark annotations to keep or discard (keep=None toggles each of them)."""
        for key in keys:
            discard = (key not in self.discarded) if keep is None else not keep
            if discard:
                self.discarded.add(key)
            else:
                self.discarded.discard(key)

    def kept(self) -> dict[AnnotationKey, dict]:
        """Annotations as they are after writing the changes: saved and pending ones that are not discarded."""
        return {key: info for collection in (self.saved, self.pending) for key, info in collection.items()
                if key not in self.discarded}

    def annotated_feature_ids(self) -> set[int]:
        """Measured features with an annotation that is (or will be) in the database."""
        return {key[1] for collection in (self.saved, self.pending) for key in collection
                if key not in self.discarded}

    def feature_annotation(self, f_id: int) -> dict | None:
        """Annotation of a measured feature (pending or saved): manual > library match > in-source fragment."""
        annotations: dict[str, dict] = {}
        for collection in (self.saved, self.pending):
            for key, info in collection.items():
                if key[1] == f_id and key not in self.discarded:
                    annotations.setdefault(key[0], info)
        for kind in ("manual", "library", "in_source_fragment"):
            if kind in annotations:
                return annotations[kind]
        return None

    def _write_annotation(self, meas_db: FeatureManagerDB, key: AnnotationKey):
        """Add a pending annotation to the measured database."""
        kind, f_id_meas, _ = key
        info: dict = self.pending[key]
        if kind == "library":
            meas_db.add_annotations_from_library(
                library=info['library'], f_id_meas=f_id_meas, f_id_lib=info['reference'], ms2_score=info['ms2_score']
            )
        elif kind == "manual":
            meas_db.add_manual_annotation(
                f_id_meas, name=info['name'], formula=info['formula'], adduct=info['adduct'], smiles=info['smiles']
            )
        elif kind == "in_source_fragment":
            meas_db.add_in_source_fragment_annotation(f_id_meas, info['reference'], note=info['note'])
        else:
            raise ValueError(f'unknown annotation kind {kind}')

    def write(self, meas_db: FeatureManagerDB) -> tuple[int, int, list[str]]:
        """
        Add the new and remove the discarded annotations in the measured database. Failed additions stay pending.
        Returns the number of added and removed annotations and the error messages.
        """
        to_add, to_remove = self.unsaved_changes()
        errors: list[str] = []
        failed_additions: set[AnnotationKey] = set()
        n_removed: int = 0
        for key in to_add:
            try:
                self._write_annotation(meas_db, key)
            except Exception as error:
                failed_additions.add(key)
                errors.append(f"adding {ANNOTATION_KIND_LABELS[key[0]]} to feature {key[1]}: {error}")
        for key in to_remove:
            try:
                meas_db.remove_annotations_by_id(self.saved[key]['candidate_ids'])
                n_removed += 1
            except Exception as error:
                errors.append(f"removing {ANNOTATION_KIND_LABELS[key[0]]} from feature {key[1]}: {error}")

        # everything else is now reflected by the database
        self.pending = {key: value for key, value in self.pending.items() if key in failed_additions}
        self.discarded.clear()
        self.load_saved(meas_db)
        return len(to_add) - len(failed_additions), n_removed, errors
