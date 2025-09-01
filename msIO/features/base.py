from sqlalchemy.orm import DeclarativeBase


CONVERTABLE_TYPES = {int, float, bool, str}


def get_py_dtypes_for_obj(obj: object) -> dict[str, type]:
    dtypes: dict[str, type] = {}
    for attr_name, mappable in obj.__annotations__.items():
        # pull out the dtype from Mapped
        t = mappable.__args__[0]
        # pull out arg from optional
        if hasattr(t, '_name') and (t._name == 'Optional'):
            t = t.__args__[0]

        dtypes[attr_name] = t
    return dtypes


class SqlBaseClass(DeclarativeBase):
    pass


class FeatureBaseClass:
    """Some universal functionality for Feature objects."""

    @classmethod
    def py_types(cls) -> dict[str, type]:
        return get_py_dtypes_for_obj(cls)

    @classmethod
    def _convert_type(cls, attr: str, val):
        """Use annotations to get desired type"""
        f = cls.py_types()[attr]
        if f not in CONVERTABLE_TYPES:
            return val
        if f is int:
            try:
                v = f(val)
            except ValueError:
                v = None
            return v
        return f(val)

    def __init__(self, **kwargs):
        super().__init__()

        # types_from_annotations: dict[str, type] = self.py_types()
        #
        # def do_nothing_conversion(val):
        #     return val
        #
        # def convert_func(attr, val):
        #     f = types_from_annotations.get(attr, do_nothing_conversion)
        #     print(f'attempting to convert {attr} with {val=} to {f}')
        #     return f(val)

        # use type annotations to convert input kwargs to right types
        kwargs_converted = {k: self._convert_type(k, v) for k, v in kwargs.items()}

        for k, v in kwargs_converted.items():
            setattr(self, k, v)  # ensures ORM descriptors are used
        if getattr(self, 'rt_minutes', None) is not None and getattr(self, 'rt_seconds', None) is None:
            self.rt_seconds = self.rt_minutes * 60
        elif getattr(self, 'rt_seconds', None) is not None and getattr(self, 'rt_minutes', None) is None:
            self.rt_minutes = self.rt_seconds / 60

        # define hook manually such that we don't need to rely on dataclass
        self.__post_init__()

    def __post_init__(self):
        pass


if __name__ == '__main__':
    # test type conversion and other population
    fbc1 = FeatureBaseClass(rt_seconds=60)
    fbc2 = FeatureBaseClass(rt_seconds="60")

    assert fbc1.rt_minutes == fbc2.rt_minutes == 1
