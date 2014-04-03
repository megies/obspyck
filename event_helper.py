import re
import obspy.core.event
from obspy.core.event import WaveformStreamID, ResourceIdentifier, TimeWindow

ID_ROOT = "smi:de.erdbeben-in-bayern"
# these classes get subclassed and need to be patched in readEvents
CLASSES_TO_PATCH = ['FocalMechanism',
    'StationMagnitudeContribution', 'StationMagnitude', 'Magnitude', 'Catalog',
    'Event', 'Origin', 'Pick', 'Arrival', 'Amplitude']


def camelcase2lower(name):
    """
    Convert CamelCase to lower_case_with_underscores.
    """
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


def newResourceIdentifier(class_name):
    id_head = "/".join((ID_ROOT, class_name))
    return ResourceIdentifier(prefix=id_head)


class CommonEventHelper():
    """
    Some common helper methods for Event type classes.
    """
    def newID(self):
        """
        Set new resource_id.
        """
        class_name = camelcase2lower(self.__class__.__name__)
        self.resource_id = newResourceIdentifier(class_name)


class FocalMechanism(obspy.core.event.FocalMechanism, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(FocalMechanism, self).__init__()
        self.newID()


class StationMagnitudeContribution(
        obspy.core.event.StationMagnitudeContribution, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(StationMagnitudeContribution, self).__init__()
        self.newID()


class StationMagnitude(obspy.core.event.StationMagnitude, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(StationMagnitude, self).__init__()
        self.newID()


class Magnitude(obspy.core.event.Magnitude, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(Magnitude, self).__init__()
        self.newID()


class Catalog(obspy.core.event.Catalog, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(Catalog, self).__init__()
        self.newID()


class Event(obspy.core.event.Event, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(Event, self).__init__()
        self.newID()


class Origin(obspy.core.event.Origin, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(Origin, self).__init__()
        self.newID()


class Pick(obspy.core.event.Pick, CommonEventHelper):
    def __init__(self, seed_string=None, phase_hint=None, *args, **kwargs):
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

    def setTime(self, time):
        self.time = time

    def setErrorTime(self, time):
        """
        Set upper or lower uncertainty, do nothing if pick has no time.
        """
        if not self.time:
            return
        delta = abs(time - self.time)
        # Determine if left or right Error Pick
        if time < self.time:
            self.time_errors.lower_uncertainty = delta
        elif time > self.time:
            self.time_errors.upper_uncertainty = delta
        # changing a subproperty, need to manually set new resource_id
        self.newID()


class Arrival(obspy.core.event.Arrival, CommonEventHelper):
    def __init__(self, origin=None, pick=None, *args, **kwargs):
        super(Arrival, self).__init__()
        self.newID()
        if origin:
            origin.arrivals.append(self)
        if pick:
            self.pick_id = pick.resource_id


class Amplitude(obspy.core.event.Amplitude, CommonEventHelper):
    def __init__(self, seed_string=None, *args, **kwargs):
        super(Amplitude, self).__init__()
        if seed_string:
            self.waveform_id = WaveformStreamID(seed_string=seed_string)
        self.newID()
        self.low = None
        self.high = None
        self.low_time = None
        self.high_time = None
        self.time_window = TimeWindow()

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
        return super(Amplitude, self).__setattr__(name, value)

    def setLow(self, time, value):
        self.low = value
        self.low_time = time
        self.updateValue()
        self.updateTimeWindow()

    def setHigh(self, time, value):
        self.high = value
        self.high_time = time
        self.updateValue()
        self.updateTimeWindow()

    def updateValue(self):
        if self.low and self.high:
            self.generic_amplitude = self.high - self.low
        else:
            self.generic_amplitude = None

    def updateTimeWindow(self):
        tw = self.time_window
        if self.low_time and self.high_time:
            tw.reference = self.low_time
            diff = self.high_time - self.low_time
            absdiff = abs(diff)
            if diff >= 0:
                tw.begin = 0.0
                tw.end = absdiff
            else:
                tw.begin = absdiff
                tw.end = 0.0
        else:
            tw.clear()


local = locals()


def readQuakeML(*args, **kwargs):
    """
    Patched readEvents function from obspy that creates instances of our
    subclassed event classes instead of the original obspy classes.
    """
    bkp = {}
    # replace original event classes with subclasses
    for classname in CLASSES_TO_PATCH:
        bkp[classname] = obspy.core.event.__dict__[classname]
        obspy.core.quakeml.__dict__[classname] = local[classname]
    from obspy.core.quakeml import readQuakeML
    ret = obspy.core.quakeml.readQuakeML(*args, **kwargs)
    # reset original event classes
    for classname, class_ in bkp.iteritems():
        obspy.core.quakeml.__dict__[classname] = bkp[classname]
    return ret
