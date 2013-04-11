from obspy.core.event import Catalog, Event, Origin, Pick, Arrival, \
    Magnitude, StationMagnitude, StationMagnitudeContribution, \
    FocalMechanism, ResourceIdentifier, WaveformStreamID


ID_ROOT = "smi:de.erdbeben-in-bayern"


class ResourceIdentifier(ResourceIdentifier):
    def __init__(self, class_name):
        id_head = "/".join((ID_ROOT, class_name))
        super(ResourceIdentifier, self).__init__(prefix=id_head)


class FocalMechanism(FocalMechanism):
    def __init__(self):
        super(FocalMechanism, self).__init__()
        self.resource_id = ResourceIdentifier("focal_mechanism")


class StationMagnitudeContribution(StationMagnitudeContribution):
    def __init__(self):
        super(StationMagnitudeContribution, self).__init__()
        self.resource_id = ResourceIdentifier(
            "station_magnitude_contribution")


class StationMagnitude(StationMagnitude):
    def __init__(self):
        super(StationMagnitude, self).__init__()
        self.resource_id = ResourceIdentifier("station_magnitude")


class Magnitude(Magnitude):
    def __init__(self):
        super(Magnitude, self).__init__()
        self.resource_id = ResourceIdentifier("magnitude")


class Catalog(Catalog):
    def __init__(self):
        super(Catalog, self).__init__()
        self.resource_id = ResourceIdentifier("catalog")


class Event(Event):
    def __init__(self):
        super(Event, self).__init__()
        self.resource_id = ResourceIdentifier("event")


class Origin(Origin):
    def __init__(self):
        super(Origin, self).__init__()
        self.resource_id = ResourceIdentifier("origin")


class Pick(Pick):
    def __init__(self, seed_string=None, phase_hint=None):
        super(Pick, self).__init__()
        if seed_string:
            self.waveform_id = WaveformStreamID(seed_string=seed_string)
        if phase_hint:
            self.phase_hint = phase_hint
        self.resource_id = ResourceIdentifier("pick")

    def newID(self):
        self.resource_id = ResourceIdentifier("pick")


class Arrival(Arrival):
    def __init__(self, origin=None, pick=None):
        super(Arrival, self).__init__()
        self.resource_id = ResourceIdentifier("arrival")
        if origin:
            origin.arrivals.append(self)
        if pick:
            self.pick_id = pick.resource_id
