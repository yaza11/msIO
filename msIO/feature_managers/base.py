from typing import Iterable

from tqdm import tqdm

from msIO.features.base import FeatureBaseClass


class FeatureManager:
    _features: dict[int, FeatureBaseClass] = None

    @property
    def feature_ids(self):
        raise NotImplementedError()

    def _inner_missing_feature(self, f_id) -> None:
        raise NotImplementedError()

    def _fetch_missing_features(self, feature_ids: Iterable[int]):
        if self._features is None:
            self._features = {}

        # exclude non-existent ids and those already processed
        missing_features: set[int] = set(feature_ids).difference(set(self._features.keys())) & set(self.feature_ids)

        for f_id in tqdm(missing_features,
                         desc=f'getting values from {self.__class__.__name__}',
                         total=(n_missing := len(missing_features)),
                         disable=n_missing <= 10):
            self._inner_missing_feature(f_id)

    def clear_cache(self):
        self._features = None

    def get_feature(self, feature_id: int) -> FeatureBaseClass:
        feature_id = int(feature_id)  # need to convert e.g. numpy int to native int
        assert feature_id in self.feature_ids, f'found no feature of {self.__class__.__name__} with {feature_id=}'
        self._fetch_missing_features([feature_id])
        return self._features[feature_id]

    def get_features(self, feature_ids: Iterable[int] = None) -> dict[int, FeatureBaseClass]:
        if feature_ids is None:
            feature_ids = self.feature_ids

        self._fetch_missing_features(feature_ids)
        return {k: v for k, v in self._features.items() if k in feature_ids}
