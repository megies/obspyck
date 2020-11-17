# -*- coding: utf-8 -*-
import re
import warnings
from copy import deepcopy

import obspy.core.event
import obspy.io
import obspy.io.quakeml
import obspy.io.quakeml.core
from obspy import UTCDateTime
from obspy.core.event import WaveformStreamID, ResourceIdentifier, \
    TimeWindow, CreationInfo, Comment

from .util import VERSION_INFO

ID_ROOT = "smi:de.erdbeben-in-bayern"
AGENCY_ID = "Erdbebendienst Bayern"
AGENCY_URI = "%s/agency" % ID_ROOT
# these classes get subclassed and need to be patched in readEvents
CLASSES_TO_PATCH = [
    'FocalMechanism', 'StationMagnitudeContribution', 'StationMagnitude',
    'Magnitude', 'Catalog', 'Event', 'Origin', 'Pick', 'Arrival', 'Amplitude']


# we're setting some attributes for internal purposes and want to ignore those
# warnings:
warnings.filterwarnings(
    action='ignore', module=r'obspy\.core\.util\.attribdict',
    category=UserWarning,
    message=r'Setting attribute .* which is not a default attribute .*')


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

    def __set_creation_info(self):
        self.creation_info = CreationInfo()
        self.creation_info.agency_id = AGENCY_ID
        self.creation_info.agency_uri = AGENCY_URI
        self.creation_info.creation_time = UTCDateTime()
        self.creation_info.version = VERSION_INFO


class FocalMechanism(obspy.core.event.FocalMechanism, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(FocalMechanism, self).__init__()
        self.newID()
        self._CommonEventHelper__set_creation_info()


class StationMagnitudeContribution(
        obspy.core.event.StationMagnitudeContribution, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(StationMagnitudeContribution, self).__init__()
        self.newID()


class StationMagnitude(obspy.core.event.StationMagnitude, CommonEventHelper):
    # we use _amplitude_ids to store amplitude IDS of used Amplitudes, since
    # QuakeML only has room to store a single Amplitude for a StationMagnitude
    # unfortunately
    do_not_warn_on = ['extra', '_amplitude_ids']

    def __init__(self, *args, **kwargs):
        super(StationMagnitude, self).__init__()
        self.used = True
        self.newID()
        self._CommonEventHelper__set_creation_info()


class Magnitude(obspy.core.event.Magnitude, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(Magnitude, self).__init__()
        self.newID()
        self._CommonEventHelper__set_creation_info()


class Catalog(obspy.core.event.Catalog, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(Catalog, self).__init__()
        self.newID()
        self._CommonEventHelper__set_creation_info()


class Event(obspy.core.event.Event, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(Event, self).__init__()
        self.newID()
        self._CommonEventHelper__set_creation_info()

    def set_creation_info_username(self, username):
        if not self.creation_info:
            self._CommonEventHelper__set_creation_info()
        self.creation_info.author = username


class Origin(obspy.core.event.Origin, CommonEventHelper):
    def __init__(self, *args, **kwargs):
        super(Origin, self).__init__()
        self.newID()
        self._CommonEventHelper__set_creation_info()


class Pick(obspy.core.event.Pick, CommonEventHelper):
    def __init__(self, seed_string=None, phase_hint=None, *args, **kwargs):
        super(Pick, self).__init__()
        if seed_string:
            self.waveform_id = WaveformStreamID(seed_string=seed_string)
        if phase_hint:
            self.phase_hint = phase_hint
        self.newID()
        self._CommonEventHelper__set_creation_info()

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
        self._CommonEventHelper__set_creation_info()


class Amplitude(obspy.core.event.Amplitude, CommonEventHelper):
    # we use _station_magnitude_id to store ID of derived StationMagnitude,
    # since QuakeML does not link back and this makes it easier to find the
    # correct one
    do_not_warn_on = ['extra', '_station_magnitude_id']

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
        self._CommonEventHelper__set_creation_info()

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
        self.update()

    def setHigh(self, time, value):
        self.high = value
        self.high_time = time
        self.update()

    def setFromTimeWindow(self, tr):
        """
        Set all values internally from time window.
        :type tr: :class:`~obspy.core.trace.Trace`
        """
        tr = tr.copy()
        t_min = self.time_window.reference
        begin = self.time_window.begin
        end = self.time_window.end
        if begin == 0 and end > 0:
            t_max = t_min + end
            tr = tr.trim(t_min, t_max, nearest_sample=True)
            self.low = tr.data[0]
            self.high = tr.data[-1]
        elif end == 0 and begin > 0:
            t_max = t_min - begin
            tr = tr.trim(t_max, t_min, nearest_sample=True)
            self.low = tr.data[-1]
            self.high = tr.data[0]
        else:
            raise NotImplementedError()
        self.updateValue()
        self.low_time = t_min
        self.high_time = t_max
        self.updatePeriod()
        # XXX TODO "tw_bkp"-sanity check should be removed
        tw_bkp = deepcopy(self.time_window)
        self.updateTimeWindow()
        self.set_general_info()
        if tw_bkp != self.time_window:
            raise NotImplementedError()

    def update(self):
        self.updateValue()
        self.updateTimeWindow()
        self.updatePeriod()

    def updateValue(self):
        if self.low and self.high:
            self.generic_amplitude = self.high - self.low
        else:
            self.generic_amplitude = None

    def updatePeriod(self):
        if not (self.low_time and self.high_time):
            return
        period = 2.0 * abs(self.low_time - self.high_time)
        if 2.0 * (self.time_window.begin + self.time_window.end) != period:
            msg = "inconsistency in amplitude time handling!!!!"
            raise Exception(msg)
        self.period = period

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

    def get_p2p(self):
        if self.low is None or self.high is None:
            return None
        return self.high - self.low

    def get_timedelta(self):
        tw = self.time_window
        if self.low_time is None or self.high_time is None:
            return None
        return abs(self.low_time - self.high_time)

    def set_general_info(self):
        self.method_id = "/".join(
            [ID_ROOT, "amplitude_method", "obspyck", "2"])
        self.unit = "dimensionless"
        self.comments = [Comment(text="peak-to-peak amplitude in raw counts")]


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
        obspy.io.quakeml.core.__dict__[classname] = local[classname]
    from obspy.io.quakeml.core import _read_quakeml
    ret = obspy.io.quakeml.core._read_quakeml(*args, **kwargs)
    # reset original event classes
    for classname, class_ in bkp.items():
        obspy.io.quakeml.core.__dict__[classname] = bkp[classname]
    return ret


def merge_events_in_catalog(catalog):
    """
    Merge information contained in all events (origins, picks, magnitudes, ...)
    in the catalog into the first event in the catalog.
    WARNING: Works in place!
    """
    if len(catalog) < 2:
        return
    event = Event()
    for e in catalog:
        event.picks += e.picks
        event.origins += e.origins
        event.magnitudes += e.magnitudes
        event.amplitudes += e.amplitudes
        event.focal_mechanisms += e.focal_mechanisms
        event.station_magnitudes += e.station_magnitudes
        event.comments += e.comments
    catalog.events = [e]
