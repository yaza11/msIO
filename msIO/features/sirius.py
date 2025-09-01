from dataclasses import dataclass

from typing import Optional, Self
from sqlalchemy import ForeignKey, String, Float, Integer, Boolean
from sqlalchemy.orm import Mapped, mapped_column, relationship

import numpy as np
import pandas as pd

from msIO.features.base import FeatureBaseClass, SqlBaseClass, CONVERTABLE_TYPES


class FormulaCandidate(SqlBaseClass, FeatureBaseClass):
    __tablename__ = "formula_candidate"

    id: Mapped[int] = mapped_column(primary_key=True)
    formula_sirius: Mapped[str] = mapped_column(String, nullable=True)
    formula_rank: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    adduct_sirius: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    zodiac_score: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    sirius_score: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    tree_score: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    isotope_score: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    num_explained_peaks: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    explained_intensity: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    median_mass_error_fragments_ppm: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    mass_error_precursor_ppm: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    lipid_class: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    rt_seconds: Mapped[Optional[float]] = mapped_column(Integer, nullable=True)
    sirius_compound_folder: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    feature_id: Mapped[int] = mapped_column(ForeignKey("feature_sirius.feature_id"))
    feature: Mapped["FeatureSirius"] = relationship(back_populates="formula_candidates")


class CompoundCandidate(SqlBaseClass, FeatureBaseClass):
    __tablename__ = "compound_candidate"

    id: Mapped[int] = mapped_column(primary_key=True)
    confidence_rank: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    structure_per_id_rank: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    formula_rank: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    num_adducts: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    num_predicted_fingerprints: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    confidence_score: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    finger_id_score: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    zodiac_score: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    sirius_score: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    formula_sirius: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    adduct_sirius: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    inchi: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    name_sirius: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    smiles: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    xlogp: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    rt_seconds: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    sirius_compound_folder: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    feature_id: Mapped[int] = mapped_column(ForeignKey("feature_sirius.feature_id"))
    feature: Mapped["FeatureSirius"] = relationship(back_populates="compound_candidates")


class CompoundGroup(SqlBaseClass, FeatureBaseClass):
    __tablename__ = "compound_group"

    id: Mapped[int] = mapped_column(primary_key=True)
    sirius_compound_folder: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    formula_sirius: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    adduct_sirius: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    npc_pathway_name: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    npc_pathway_probability: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    npc_superclass_name: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    npc_superclass_probability: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    npc_class_name: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    npc_class_probability: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    cf_most_specific_name: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    cf_most_specific_probability: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    cf_level5_name: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    cf_level5_probability: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    cf_subclass_name: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    cf_subclass_probability: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    cf_class_name: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    cf_class_probability: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    cf_superclass_name: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    cf_superclass_probability: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    cf_path: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    feature_id: Mapped[int] = mapped_column(ForeignKey("feature_sirius.feature_id"))
    feature: Mapped["FeatureSirius"] = relationship(back_populates="compound_groups")


class FeatureSirius(SqlBaseClass, FeatureBaseClass):
    __tablename__ = "feature_sirius"

    feature_id: Mapped[int] = mapped_column(primary_key=True)

    formula_candidates: Mapped[list[FormulaCandidate]] = relationship(
        back_populates="feature",
        cascade="all, delete-orphan"
    )
    compound_candidates: Mapped[list[CompoundCandidate]] = relationship(
        back_populates="feature",
        cascade="all, delete-orphan"
    )
    compound_groups: Mapped[list[CompoundGroup]] = relationship(
        back_populates="feature",
        cascade="all, delete-orphan"
    )

    use_zodiac_scoring_for_best: Mapped[bool] = mapped_column(Boolean, default=True)
    highest_scoring_formula: Mapped[Optional[String]] = mapped_column(String, nullable=True)

    combined_feature_id: Mapped[Optional[int]] = mapped_column(ForeignKey('features.id'))
    combined_feature: Mapped[Optional["FeatureCombined"]] = relationship(back_populates='sirius')

    def __init__(self, use_zodiac_scoring_for_best: bool = True, **kwargs):
        super().__init__(**kwargs)
        # TODO: find better method to set default value
        self.use_zodiac_scoring_for_best = use_zodiac_scoring_for_best

    def formula_candidates_by_rank(self) -> dict[int, FormulaCandidate]:
        return {fc.formula_rank: fc for fc in self.formula_candidates if fc.formula_rank is not None}

    def compound_candidates_by_formula(self) -> dict[str, CompoundCandidate]:
        return {cc.formula_sirius: cc for cc in self.compound_candidates if cc.formula_sirius}

    def compound_groups_by_formula(self) -> dict[str, CompoundGroup]:
        return {cg.formula_sirius: cg for cg in self.compound_groups if cg.formula_sirius}

    @classmethod
    def from_tables(cls, feature_id: int, tables: dict[str, pd.DataFrame]) -> Self:
        def select_rows_by_feature(df: pd.DataFrame) -> pd.DataFrame:
            mask = df.feature_id == feature_id
            return df.loc[mask, :]

        formula_candidates = []
        sub_df = select_rows_by_feature(tables['formula_identifications'])
        for _, row in sub_df.iterrows():
            f = FormulaCandidate(**row.to_dict())
            formula_candidates.append(f)

        compound_candidates = []
        sub_df = select_rows_by_feature(tables['compound_identifications'])
        for _, row in sub_df.iterrows():
            f = CompoundCandidate(**row.to_dict())
            compound_candidates.append(f)
        # TODO: fill nan values, where possible

        compound_groups = []
        sub_df = select_rows_by_feature(tables['canopus_formula_summary'])
        for _, row in sub_df.iterrows():
            f = CompoundGroup(**row.to_dict())
            compound_groups.append(f)

        return cls(feature_id=feature_id,
                   formula_candidates=formula_candidates,
                   compound_candidates=compound_candidates,
                   compound_groups=compound_groups)

    def _set_highest_scoring_formula(self):
        scores: list[float] = [c.zodiac_score if self.use_zodiac_scoring_for_best else c.sirius_score
                               for c in self.formula_candidates]
        idx = int(np.nanargmax(scores))
        best_candidate = self.formula_candidates[idx]
        self.highest_scoring_formula = best_candidate.formula_sirius

    def get_highest_scoring(self, prefer_npc: bool = True) -> Self:
        """flatten attributes by taking properties from highest ranked formula,
        returns new instance (shallow copy with added attributes)"""
        def add_from_highest_scoring(obj) -> None:
            new.__dict__ |= obj.__dict__

        if self.highest_scoring_formula is None:
            self._set_highest_scoring_formula()

        new = self.__class__()
        new.__dict__ |= self.__dict__

        # add attributes of objects in reverse order of chain in inference as later attributes override previous ones
        # find best for compound groups based on probability
        candidates_compound_groups = [c for c in self.compound_groups
                                      if c.formula_sirius == self.highest_scoring_formula]
        if len(candidates_compound_groups) > 0:
            # if there are multiple candidates, use the one with the highest probability
            prob_attrs = [a for a in CompoundGroup.__annotations__.keys()
                          if a.endswith('probability') and not (a.startswith('npc') ^ prefer_npc)]
            candidates_probabilities = [max([c.__getattribute__(a) for a in prob_attrs])
                                        for c in candidates_compound_groups]
            obj_idx = candidates_probabilities.index(max(candidates_probabilities))
            add_from_highest_scoring(candidates_compound_groups[obj_idx])

        candidates_compounds = [c for c in self.compound_candidates
                                if c.formula_sirius == self.highest_scoring_formula]
        # find best candidate based on confidence score or sirius/zodiac score if all confidence scores are nan
        if len(candidates_compounds) > 0:
            if any([c.confidence_score > 0 for c in candidates_compounds]):
                scores = [c.confidence_score for c in candidates_compounds]
            elif self.use_zodiac_scoring_for_best:
                scores = [c.zodiac_score for c in candidates_compounds]
            else:
                scores = [c.sirius_score for c in candidates_compounds]
            idx = np.nanargmax(scores)
            add_from_highest_scoring(candidates_compounds[idx])

        candidates_formula = [c for c in self.formula_candidates
                              if c.formula_sirius == self.highest_scoring_formula]
        # there should always be at least 1
        scores = [c.zodiac_score
                  if self.use_zodiac_scoring_for_best
                  else c.sirius_score
                  for c in candidates_formula]
        add_from_highest_scoring(candidates_formula[np.nanargmax(scores)])

        return new


if __name__ == '__main__':
    pass
    f = FeatureSirius()
