from dataclasses import dataclass

import numpy as np
import pandas as pd

from msIO.features.base import FeatureBaseClass


@dataclass
class FormulaCandidate:
    feature_id: int
    formula_sirius: str
    formula_rank: int = None
    adduct_sirius: str = None
    zodiac_score: float = None
    sirius_score: float = None
    tree_score: float = None
    isotope_score: float = None
    num_explained_peaks: int = None
    explained_intensity: float = None
    median_mass_error_fragments_ppm: float = None
    mass_error_precursor_ppm: float = None
    lipid_class: str = None
    rt_seconds: float = None
    sirius_compound_folder: str = None


@dataclass
class CompoundCandidate:
    confidence_rank: int = None
    structure_per_id_rank: int = None
    formula_rank: int = None
    num_adducts: int = None
    num_predicted_fingerprints: int = None
    confidence_score: float = None
    finger_id_score: float = None
    zodiac_score: float = None
    sirius_score: float = None
    formula_sirius: str = None
    adduct_sirius: str = None
    inchi: str = None
    name_sirius: str = None
    smiles: str = None
    xlogp: float = None
    rt_seconds: float = None
    feature_id: int = None
    sirius_compound_folder: str = None


@dataclass
class CompoundGroup:
    sirius_compound_folder: str = None
    formula_sirius: str = None
    adduct_sirius: str = None
    npc_pathway_name: str = None
    npc_pathway_probability: float = None
    npc_superclass_name: str = None
    npc_superclass_probability: float = None
    npc_class_name: str = None
    npc_class_probability: float = None
    cf_most_specific_name: str = None
    cf_most_specific_probability: float = None
    cf_level5_name: str = None
    cf_level5_probability: float = None
    cf_subclass_name: str = None
    cf_subclass_probability: float = None
    cf_class_name: str = None
    cf_class_probability: float = None
    cf_superclass_name: str = None
    cf_superclass_probability: float = None
    cf_path: str = None
    feature_id: int = None


@dataclass
class FeatureSirius(FeatureBaseClass):
    feature_id: int
    formula_candidates: dict[int, FormulaCandidate] = None  # key: rank
    compound_candidates: dict[str, CompoundCandidate] = None  # key: formula
    compound_groups: dict[str, CompoundGroup] = None  # key: formula

    use_zodiac_scoring_for_best: bool = True
    highest_scoring_candidate_rank: int = None

    @classmethod
    def from_tables(cls, feature_id: int, tables: dict[str, pd.DataFrame]):
        def select_rows_by_feature(df: pd.DataFrame) -> pd.DataFrame:
            mask = df.feature_id == feature_id
            return df.loc[mask, :]

        formula_candidates = {}
        sub_df = select_rows_by_feature(tables['formula_identifications'])
        for _, row in sub_df.iterrows():
            f = FormulaCandidate(**row.to_dict())
            formula_candidates[row.formula_rank] = f

        compound_candidates = {}
        sub_df = select_rows_by_feature(tables['compound_identifications'])
        for _, row in sub_df.iterrows():
            f = CompoundCandidate(**row.to_dict())
            compound_candidates[row.formula_sirius] = f
        # TODO: fill nan values, where possible

        compound_groups = {}
        sub_df = select_rows_by_feature(tables['canopus_formula_summary'])
        for _, row in sub_df.iterrows():
            f = CompoundGroup(**row.to_dict())
            compound_groups[row.formula_sirius] = f

        return cls(feature_id=feature_id,
                   formula_candidates=formula_candidates,
                   compound_candidates=compound_candidates,
                   compound_groups=compound_groups)

    def __post_init__(self):
        """flatten attributes by taking properties from highest ranked formula"""
        add_attributes = {}
        # add attributes from highest scoring formula
        scores: list[float] = [c.zodiac_score if self.use_zodiac_scoring_for_best else c.sirius_score
                               for r, c in self.formula_candidates.items()]
        idx = np.argmax(scores)
        rank: int = list(self.formula_candidates.keys())[idx]
        self.highest_scoring_candidate_rank = rank
        candidate = self.formula_candidates[rank]
        formula: str = candidate.formula_sirius

        for obj_dict in [self.compound_groups, self.compound_candidates]:
            if formula in obj_dict:
                add_attributes |= obj_dict[formula].__dict__
        add_attributes |= self.formula_candidates[rank].__dict__

        self.__dict__ |= add_attributes

