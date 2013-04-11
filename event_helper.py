import re
import obspy.core.event
from obspy.core.event import WaveformStreamID

ID_ROOT = "smi:de.erdbeben-in-bayern"


def camelcase2lower(name):
    """
    Convert CamelCase to lower_case_with_underscores.
    """
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


class CommonEventHelper():
    """
    Some common helper methods for Event type classes.
    """
    def newID(self):
        """
        Set new resource_id.
        """
        name = camelcase2lower(self.__class__.__name__)
        self.resource_id = ResourceIdentifier(name)


class ResourceIdentifier(obspy.core.event.ResourceIdentifier):
    def __init__(self, class_name):
        id_head = "/".join((ID_ROOT, class_name))
        super(ResourceIdentifier, self).__init__(prefix=id_head)


class FocalMechanism(obspy.core.event.FocalMechanism, CommonEventHelper):
    def __init__(self):
        super(FocalMechanism, self).__init__()
        self.newID()


class StationMagnitudeContribution(
        obspy.core.event.StationMagnitudeContribution, CommonEventHelper):
    def __init__(self):
        super(StationMagnitudeContribution, self).__init__()
        self.newID()


class StationMagnitude(obspy.core.event.StationMagnitude, CommonEventHelper):
    def __init__(self):
        super(StationMagnitude, self).__init__()
        self.newID()


class Magnitude(obspy.core.event.Magnitude, CommonEventHelper):
    def __init__(self):
        super(Magnitude, self).__init__()
        self.newID()


class Catalog(obspy.core.event.Catalog, CommonEventHelper):
    def __init__(self):
        super(Catalog, self).__init__()
        self.newID()


class Event(obspy.core.event.Event, CommonEventHelper):
    def __init__(self):
        super(Event, self).__init__()
        self.newID()


class Origin(obspy.core.event.Origin, CommonEventHelper):
    def __init__(self):
        super(Origin, self).__init__()
        self.newID()


class Pick(obspy.core.event.Pick, CommonEventHelper):
    def __init__(self, seed_string=None, phase_hint=None):
        super(Pick, self).__init__()
        if seed_string:
            self.waveform_id = WaveformStreamID(seed_string=seed_string)
        if phase_hint:
            self.phase_hint = phase_hint
        self.newID()

    def __setattr__(self, name, value):
        """
        Set new resource_id on any attribute change (other than setting a new
        resource_id).

        XXX TODO if we do all attribute changes in setter methods here, we can
        probably take care of this in the setter methods and avoid this
        override?!
        """
        if name != "resource_id":
            self.newID()
        return super(Pick, self).__setattr__(name, value)


class Arrival(obspy.core.event.Arrival, CommonEventHelper):
    def __init__(self, origin=None, pick=None):
        super(Arrival, self).__init__()
        self.newID()
        if origin:
            origin.arrivals.append(self)
        if pick:
            self.pick_id = pick.resource_id
