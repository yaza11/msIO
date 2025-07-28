class FeatureBaseClass:
    rt_seconds: float = None
    rt_minutes: float = None

    @classmethod
    def _convert_type(cls, attr: str, val):
        """Use annotations to get desired type"""
        return cls.__annotations__[attr](val)

    def __post_init__(self):
        if not self.rt_seconds and self.rt_minutes:
            self.rt_seconds = self.rt_minutes * 60

        if not self.rt_minutes and self.rt_seconds:
            self.rt_minutes = self.rt_seconds / 60
