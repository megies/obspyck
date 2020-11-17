# -*- coding: utf-8 -*-
# -------------------------------------------------------------------
# Filename: util.py
#  Purpose: Helper functions for ObsPyck
#   Author: Tobias Megies, Lion Krischer
#    Email: megies@geophysik.uni-muenchen.de
#  License: GPLv2
#
# Copyright (C) 2010 Tobias Megies, Lion Krischer
# -------------------------------------------------------------------
import copy
import glob
import io
import math
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import warnings

from PyQt5 import QtWidgets
import numpy as np
import matplotlib as mpl
from matplotlib.colors import ColorConverter
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as QFigureCanvas
from matplotlib.widgets import MultiCursor as MplMultiCursor

import obspy
from obspy import Trace, Inventory
import obspy.clients.arclink
from obspy import UTCDateTime, read_inventory, read, Stream
from obspy.clients.arclink import Client as ArcLinkClient
from obspy.clients.fdsn import Client as FDSNClient
from obspy.clients.filesystem.sds import Client as SDSClient
from obspy.clients.seedlink import Client as SeedlinkClient
from obspy.clients.seishub import Client as SeisHubClient
from obspy.geodetics.base import gps2dist_azimuth
from obspy.io.xseed import Parser

from . import __version__
from .rotate_to_zne import (
    _rotate_specific_channels_to_zne, get_orientation_from_parser,
    get_orientation)

obspy.clients.arclink.client.MAX_REQUESTS = 200

mpl.rc('figure.subplot', left=0.05, right=0.98, bottom=0.10, top=0.92,
       hspace=0.28)
mpl.rcParams['font.size'] = 10


try:
    # obspy >= 1.1.0
    from obspy.core.util import MATPLOTLIB_VERSION
except ImportError:
    # obspy < 1.1.0
    from obspy.core.util import get_matplotlib_version
    MATPLOTLIB_VERSION = get_matplotlib_version()


# restrict to 64 characters to fit in allowed QuakeML maximum version field
# length
VERSION_INFO = 'ObsPyck {}, ObsPy {}'.format(__version__, obspy.__version__)[:64]


COMMANDLINE_OPTIONS = (
        # XXX wasn't working as expected
        #(("--debug"), {'dest': "debug", 'action': "store_true",
        #        'default': False,
        #        'help': "Switch on Ipython debugging in case of exception"}),
        (("-c", "--config"), {
            'dest': "config_file", "type": "str",
            'help': "Location of obspyck config file to use. If not "
                    "specified, obspyck looks up '~/.obspyckrc' (an example "
                    "configuration will be created if the file does not "
                    "exist)."}),
        (("-s", "--station-combination"), {
            'dest': "station_combinations", "action": "append", "default": [],
            'help': "Station combination to fetch data for (which must be "
                    "defined in config file). Can be specified multiple times "
                    "to combine station combinations."}),
        (("-i", "--id"), {
            'dest': "seed_ids", "action": "append", "default": [],
            'help': "Additional SEED ID to fetch. The ID should be specified "
                    "as a full SEED ID with a '?' for the component code "
                    "(e.g. 'GR.FUR..HH?'). The server to use for fetching the "
                    "specified SEED ID is looked up in the config file "
                    "section '[seed_id_lookup]' same as for SEED IDs "
                    "specified in station combinations. Can be specified "
                    "multiple times."}),
        (("-t", "--time"), {
            'dest': "time", "default": None,
            'help': "Starttime of seismogram to retrieve. It takes a "
                    "string which UTCDateTime can convert. E.g. "
                    "'2010-01-10T05:00:00'"}),
        (("-d", "--duration"), {
            'type': "float", 'dest': "duration", "default": None,
            'help': "Duration of seismogram in seconds"}),
        (("-o", "--starttime-offset"), {
            'type': "float", 'dest': "starttime_offset",
            'default': None, 'help': "Offset to add to specified starttime "
            "in seconds. Thus a time from an automatic picker can be used "
            "with a specified offset for the starttime. E.g. to request a "
            "waveform starting 30 seconds earlier than the specified time "
            "use '-o -30'."}),
        # (("--noevents",), {
        #     'action': "store_true",
        #     'dest': "noevents", 'default': False,
        #     'help': "Deactivate fetching event data using FDSNWS and "
        #             "plotting theoretical arrivals."}),
        (("--event",), {
            'dest': "event", 'type': "str",
            'default': '',
            'help': "Filename (or path) of event file to load."}),
        )
PROGRAMS = {
        'nlloc': {'filenames': {'exe': "NLLoc", 'phases': "nlloc.obs",
                                'summary': "nlloc.hyp",
                                'scatter': "nlloc.scat"}},
        'hyp_2000': {'filenames': {'exe': "hyp2000",'control': "bay2000.inp",
                                   'phases': "hyp2000.pha",
                                   'stations': "stations.dat",
                                   'summary': "hypo.prt"}},
        'focmec': {'filenames': {'exe': "rfocmec", 'phases': "focmec.dat",
                                 'stdout': "focmec.stdout",
                                 'summary': "focmec.out"}}}
COMPONENT_COLORS = {'Z': "k", 'N': "b", 'E': "r"}
WIDGET_NAMES = ("qToolButton_clearAll", "qToolButton_clearOrigMag",
        "qToolButton_clearFocMec", "qToolButton_doHyp2000",
        "qToolButton_doNlloc", "qComboBox_nllocModel",
        "qToolButton_doFocMec", "qToolButton_showMap",
        "qToolButton_showFocMec", "qToolButton_nextFocMec",
        "qToolButton_showWadati", "qToolButton_getNextEvent",
        "qToolButton_updateEventList", "qToolButton_sendNewEvent",
        "qToolButton_replaceEvent",
        "qToolButton_deleteEvent", "qCheckBox_public",
        "qComboBox_eventType",
        "qToolButton_sort_abc", "qToolButton_sort_distance",
        "qToolButton_previousStream", "qLabel_streamNumber",
        "qComboBox_streamName", "qToolButton_nextStream",
        "qToolButton_overview", "qComboBox_phaseType", "qToolButton_rotateLQT",
        "qToolButton_rotateZRT", "qToolButton_filter", "qToolButton_trigger",
        "qToolButton_physical_units", "qDoubleSpinBox_waterlevel",
        "qToolButton_arpicker", "qComboBox_filterType", "qCheckBox_zerophase",
        "qLabel_highpass", "qDoubleSpinBox_highpass", "qLabel_lowpass",
        "qDoubleSpinBox_lowpass",
        "qDoubleSpinBox_corners", "qLabel_corners", "qCheckBox_50Hz",
        "qTextEdit_qml", "qPushButton_qml_update",
        "qLabel_sta", "qDoubleSpinBox_sta",
        "qLabel_lta", "qDoubleSpinBox_lta", "qToolButton_spectrogram",
        "qCheckBox_spectrogramLog", "qLabel_wlen", "qDoubleSpinBox_wlen",
        "qLabel_perlap", "qDoubleSpinBox_perlap", "qPlainTextEdit_stdout",
        "qPlainTextEdit_stderr")
MAG_MARKER = {'marker': (8, 2, 0), 'edgewidth': 1.8, 'size': 20}
AXVLINEWIDTH = 1.5
# dictionary for key-bindings.
# XXX Qt:
#KEYS = {'setPick': "Key_A", 'setPickError': "Key_S", 'delPick': "Key_Q",
#        'setMagMin': "Key_A", 'setMagMax': "Key_S", 'delMagMinMax': "Key_Q",
#        'switchPhase': "Key_Control",
#        'prevStream': "Key_Y", 'nextStream': "Key_X", 'switchWheelZoomAxis': "Key_Shift",
#        'setWeight': {'Key_0': 0, 'Key_1': 1, 'Key_2': 2, 'Key_3': 3},
#        'setPol': {'Key_U': "up", 'Key_D': "down", 'Key_Plus': "poorup", 'Key_Minus': "poordown"},
#        'setOnset': {'Key_I': "impulsive", 'Key_E': "emergent"}}

ROTATE_LQT_COMP_MAP = {"Z": "L", "N": "Q", "E": "T"}
ROTATE_ZRT_COMP_MAP = {"Z": "Z", "N": "R", "E": "T"}
S_POL_MAP_ZRT = {'R': {'up': "forward", 'down': "backward",
                       'poorup': "forward", 'poordown': "backward"},
                 'T': {'up': "right", 'down': "left",
                       'poorup': "right", 'poordown': "left"}}
S_POL_PHASE_TYPE = {'R': "SV", 'T': "SH"}
POLARITY_2_FOCMEC = {'Z': {'positive': "U", 'negative': "D"},
                     'R': {'positive': "F", 'negative': "B"},
                     'T': {'positive': "R", 'negative': "L"}}

# only strings involved, so shallow copy is fine
POLARITY_CHARS = {'positive': "+", 'negative': "-", 'undecidable': "?",
                  None: "_"}
ONSET_CHARS = {'impulsive': "I", 'emergent': "E", 'questionable': "?",
               None: "_"}
LOGLEVELS = {'normal': "CRITICAL", 'verbose': "INFO", 'debug': "DEBUG",
             'quiet': 100}

ONE_SIGMA = 68.3
TWO_SIGMA = 95.4

NOT_REIMPLEMENTED_MSG = ("Feature was not reimplemented after major "
                         "change to QuakeML.")

class QMplCanvas(QFigureCanvas):
    """
    Class to represent the FigureCanvas widget.
    """
    def __init__(self, parent=None):
        # Standard Matplotlib code to generate the plot
        self.fig = Figure()
        # initialize the canvas where the Figure renders into
        QFigureCanvas.__init__(self, self.fig)
        self.setParent(parent)

def matplotlib_color_to_rgb(color):
    """
    Converts matplotlib colors to rgb.
    """
    rgb = ColorConverter().to_rgb(color)
    return [int(_i*255) for _i in rgb]

def check_keybinding_conflicts(keys):
    """
    check for conflicting keybindings. 
    we have to check twice, because keys for setting picks and magnitudes
    are allowed to interfere...
    """
    for ignored_key_list in [['setMagMin', 'setMagMax', 'delMagMinMax'],
                             ['setPick', 'setPickError', 'delPick']]:
        tmp_keys = copy.deepcopy(keys)
        tmp_keys2 = {}
        for ignored_key in ignored_key_list:
            tmp_keys.pop(ignored_key)
        while tmp_keys:
            key, item = tmp_keys.popitem()
            if isinstance(item, dict):
                while item:
                    k, v = item.popitem()
                    tmp_keys2["_".join([key, str(v)])] = k
            else:
                tmp_keys2[key] = item
        if len(set(tmp_keys2.keys())) != len(set(tmp_keys2.values())):
            err = "Interfering keybindings. Please check variable KEYS"
            raise Exception(err)


def resolve_complex_station_combination(config, section_title, time):
    """
    Resolve complex time-dependent station combination from config file
    """
    seed_ids = set()
    for key, value in config.items(section_title):
        # no value: always use that stream
        if not value:
            seed_ids.add(key)
        # two comma separated time stamps: only use stream if requested time is
        # in between given timestamps
        elif ',' in value:
            start, end = list(map(UTCDateTime, value.split(',')))
            if end < start:
                msg = ("Invalid configuration in section [{}], stream '{}'. "
                       "End time larger than start time for used time "
                       "span.").format(section_title, key)
                raise ValueError(msg)
            if start <= time <= end:
                seed_ids.add(key)
        # timestamp preceded by smaller-than sign: only use stream if requested
        # time is smaller than specified timestamp
        elif value[0] == '<':
            end = UTCDateTime(value[1:])
            if time <= end:
                seed_ids.add(key)
        # timestamp preceded by larger-than sign: only use stream if requested
        # time is after the specified timestamp
        elif value[0] == '>':
            start = UTCDateTime(value[1:])
            if time >= start:
                seed_ids.add(key)
        else:
            msg = ("Invalid configuration in section [{}], stream '{}'. "
                   "Value can be empty, two comma-separated time stamps or "
                   "one timestamp preceded by either '>' or '<'.").format(
                        section_title, key)
            raise ValueError(msg)
    return seed_ids


def _get_metadata(tr, inv):
    """
    Extract metadata for given Trace
    """
    try:
        coordinates = inv.get_coordinates(
            tr.id, tr.stats.starttime)
        orientation = get_orientation(
            inv, tr.id, tr.stats.starttime)
        response = inv.get_response(tr.id, tr.stats.starttime)
    except Exception as e:
        if str(e).startswith('No matching '):
            return None
        raise
    return response, coordinates, orientation


def _attach_metadata(st, inventories):
    """
    Attach response, coordinates and orientation to all traces in stream.
    Raise an exception if it fails. Show a warning if multiple matching
    metadata are found for any trace.
    """
    if isinstance(inventories, Inventory):
        inventories = [inventories]
    if isinstance(st, Trace):
        st = [st]
    for tr in st:
        metadata = []
        for inv in inventories:
            metadata_ = _get_metadata(tr, inv)
            if metadata_ is None:
                continue
            metadata.append(metadata_)
        if not metadata:
            msg = 'Failed to get response for {}!'.format(tr.id)
            raise Exception(msg)
        elif len(metadata) > 1:
            msg = ('Found multiple matching metadata entries for {}, using '
                   'first.').format(tr.id)
            warnings.warn(msg)
        response, coordinates, orientation = metadata[0]
        tr.stats.coordinates = coordinates
        tr.stats.orientation = orientation
        tr.stats.response = response


def fetch_waveforms_with_metadata(options, args, config):
    """
    Sets up obspy clients and fetches waveforms and metadata according to
    command line options.
    Now also fetches data via arclink if --arclink-ids is used.
    Args are tried to read as local waveform files or metadata files.

    XXX Notes: XXX
     - there is a problem in the arclink client with duplicate traces in
       fetched streams. therefore at the moment it might be necessary to use
       "-m overwrite" option.

    :returns: (dictionary with clients,
               list(:class:`obspy.core.stream.Stream`s),
               list of Inventory)
    """
    if not options.station_combinations and not options.seed_ids and not args:
        msg = ('No data to use specified. At least use one of option "-s", '
               '"-i" or specify local waveform/metadata files to load.')
        raise Exception(msg)

    no_metadata = config.getboolean("base", "no_metadata")

    time_ = UTCDateTime(options.time)

    if options.starttime_offset is None:
        t1 = time_ + config.getfloat("base", "starttime_offset")
    else:
        t1 = time_ + options.starttime_offset
    if options.duration is None:
        t2 = t1 + config.getfloat("base", "duration")
    else:
        t2 = t1 + options.duration

    seed_ids_to_fetch = set()
    for station_combination_key in options.station_combinations:
        seed_ids = config.get("station_combinations", station_combination_key)
        # starts with an exclamation mark: lookup seed ids from that section
        if seed_ids[0] == '!':
            seed_ids = resolve_complex_station_combination(
                config, seed_ids[1:], t1)
        # otherwise a list of comma-separated stream labels
        else:
            seed_ids = seed_ids.split(",")
        for seed_id in seed_ids:
            seed_ids_to_fetch.add(seed_id)
    for seed_id in options.seed_ids:
        seed_ids_to_fetch.add(seed_id)

    seed_id_lookup = {}
    seed_id_lookup_keys = config.options("seed_id_lookup")
    for seed_id in seed_ids_to_fetch:
        netstaloc = seed_id.rsplit(".", 1)[0]
        netsta = seed_id.rsplit(".", 2)[0]
        net = seed_id.rsplit(".", 3)[0]
        # look up by exact SEED ID:
        if seed_id in seed_id_lookup_keys:
            seed_id_lookup[seed_id] = config.get("seed_id_lookup", seed_id)
        # look up by SEED ID down to location code
        elif netstaloc in seed_id_lookup_keys:
            seed_id_lookup[seed_id] = config.get("seed_id_lookup", netstaloc)
        # look up by SEED ID down to station code
        elif netsta in seed_id_lookup_keys:
            seed_id_lookup[seed_id] = config.get("seed_id_lookup", netsta)
        # look up by network code
        elif net in seed_id_lookup_keys:
            seed_id_lookup[seed_id] = config.get("seed_id_lookup", net)

    clients = {}

    streams = []
    sta_fetched = set()
    # Local files:
    all_inventories = []
    inventories = []
    if args:
        print("=" * 80)
        print("Reading local files:")
        print("-" * 80)
        stream_tmp = Stream()
        for file in args:
            # try to read as metadata
            try:
                inv = read_inventory(file)
            except:
                pass
            else:
                print(f"{file}: Metadata")
                inventories.append(inv)
                continue
            # try to read as waveforms
            try:
                st = read(file, starttime=t1, endtime=t2,
                          verify_chksum=not config.getboolean("base", "ignore_gse2_chksum_error"))
            except TypeError:
                print(f"File {file} not recognized as dataless or waveform file. Skipped.")
                continue
            msg = "%s: Waveforms" % file
            if not st:
                msg += " (not matching requested time window)"
            print(msg)
            stream_tmp += st
        if not inventories and not no_metadata:
            msg = ("No station metadata for waveforms from local files. "
                   "(Set the following config option to start obspyck "
                   "regardless of missing station metadata: [base] "
                   "no_metadata = true)")
            raise Exception(msg)
        if not no_metadata:
            _attach_metadata(stream_tmp, inventories)
            for tr in stream_tmp:
                if tr.stats._format == 'GSE2':
                    apply_gse2_calib(tr)
        ids = set([(tr.stats.network, tr.stats.station, tr.stats.location) for tr in stream_tmp])
        for net, sta, loc in ids:
            stream_tmp_ = stream_tmp.select(
                network=net, station=sta, location=loc)
            # check whether to attempt rotation
            if config.has_section("rotate_channels"):
                net_sta_loc = ".".join((net, sta, loc))
                if net_sta_loc in config.options("rotate_channels"):
                    rotate_channels(stream_tmp_, net, sta, loc, config)
            streams.append(stream_tmp_)
    all_inventories += inventories

    print("=" * 80)
    print("Fetching waveforms and metadata from servers:")
    print("-" * 80)
    for seed_id, server in sorted(seed_id_lookup.items()):
        server_type = config.get(server, "type")
        if server_type not in ("seishub", "fdsn", "jane", "arclink",
                               "seedlink", "sds"):
            msg = ("Unknown server type '{}' in server definition section "
                   "'{}' in config file.").format(server_type, server)
            raise NotImplementedError(msg)
        client = connect_to_server(server, config, clients)
        net, sta, loc, cha = seed_id.split(".")
        net_sta_loc = "%s.%s.%s" % (net, sta, loc)
        # make sure we dont fetch a single station of
        # one network twice (could happen with wildcards)
        if any([char in net_sta_loc for char in '?*[]']):
            msg = ("Wildcards in SEED IDs to fetch are only allowed in "
                   "channel part: {}").format(seed_id)
            raise NotImplementedError(msg)
        if net_sta_loc in sta_fetched:
            print(f"{seed_id.ljust(15)} ({server_type}: {server}) skipped! "
                  "(Was already retrieved)")
            continue
        try:
            sys.stdout.write("\r%s (%s: %s) ..." % (
                seed_id.ljust(15), server_type, server))
            sys.stdout.flush()
            # SeisHub
            if server_type == "seishub":
                st = client.waveform.get_waveforms(
                    net, sta, loc, cha, t1, t2, apply_filter=True)
                if not no_metadata:
                    data = client.station.get_list(
                        network=net, station=sta, datetime=t1)
                    if len(data) == 0:
                        msg = "No station metadata on server."
                        raise Exception(msg)
                    inventories = []
                    for d in data:
                        bio = io.BytesIO(
                            client.station.get_resource(d['resource_name']))
                        bio.seek(0)
                        inv = read_inventory(bio, format='XSEED')
                        inventories.append(inv)
                    # look in all metadata available, but prefer metadata from
                    # same source (by putting it in front in list)
                    all_inventories = inventories + all_inventories
                    _attach_metadata(st, all_inventories)
            # ArcLink
            elif server_type == "arclink":
                st = client.get_waveforms(
                    network=net, station=sta, location=loc, channel=cha,
                    starttime=t1, endtime=t2)
                if not no_metadata:
                    inventories = []
                    for net_, sta_, loc_, cha_ in set([
                            tuple(tr.id.split(".")) for tr in st]):
                        bio = io.BytesIO()
                        client.save_response(bio, net_, sta_, loc_, cha_,
                                             t1-10, t2+10)
                        bio.seek(0)
                        inventories.append(
                            read_inventory(bio, format='SEED'))
                    # look in all metadata available, but prefer metadata from
                    # same source (by putting it in front in list)
                    all_inventories = inventories + all_inventories
                    _attach_metadata(st, all_inventories)
            # FDSN (or JANE)
            elif server_type in ("fdsn", "jane"):
                st = client.get_waveforms(
                    network=net, station=sta, location=loc, channel=cha,
                    starttime=t1, endtime=t2)
                if not no_metadata:
                    inventory = client.get_stations(
                        network=net, station=sta, location=loc,
                        level="response")
                    # look in all metadata available, but prefer metadata from
                    # same source (by putting it in front in list)
                    all_inventories = [inventory] + all_inventories
                    _attach_metadata(st, all_inventories)
            # Seedlink
            elif server_type == "seedlink":
                # XXX I think the wild card checks for net/sta/loc can be
                # XXX removed, are already done earlier..
                for wc in "*?":
                    for code in (net, sta, loc):
                        if wc in code:
                            msg = ("Wildcards are not allowed for network/"
                                   "station/location when fetching data via "
                                   "seedlink ('{}')")
                            raise ValueError(msg.format(code))
                if '*' in cha:
                    msg = ("Wildcard '*' not allowed for channel code "
                           "when fetching data via seedlink ('{}')")
                    raise ValueError(msg.format(cha))
                st = client.get_waveform(
                    network=net, station=sta, location=loc, channel=cha,
                    starttime=t1, endtime=t2)
                if not st:
                    msg = "Server returned no data."
                    raise Exception(msg)
                if not no_metadata:
                    meta_server = config.get(server, "metadata_server")
                    meta_server_type = config.get(meta_server, "type")
                    meta_client = connect_to_server(meta_server,
                                                    config, clients)
                    if meta_server_type in ('fdsn', 'jane'):
                        inventory = meta_client.get_stations(
                            network=net, station=sta, location=loc,
                            level="response")
                        # look in all metadata available, but prefer metadata from
                        # same source (by putting it in front in list)
                        all_inventories = [inventory] + all_inventories
                        _attach_metadata(st, all_inventories)
                    else:
                        raise NotImplementedError()
            # SDS
            elif server_type == "sds":
                # XXX I think the wild card checks for net/sta/loc can be
                # XXX removed, are already done earlier..
                for wc in "*?":
                    for code in (net, sta, loc):
                        if wc in code:
                            msg = ("Wildcards are not allowed for network/"
                                   "station/location when fetching data via "
                                   "SDS client ('{}')")
                            raise ValueError(msg.format(code))
                if '*' in cha:
                    msg = ("Wildcard '*' not allowed for channel code "
                           "when fetching data via SDS client ('{}')")
                    raise ValueError(msg.format(cha))
                st = client.get_waveforms(
                    network=net, station=sta, location=loc, channel=cha,
                    starttime=t1, endtime=t2)
                if not st:
                    msg = "Server returned no data."
                    raise Exception(msg)
                if not no_metadata:
                    meta_server = config.get(server, "metadata_server")
                    meta_server_type = config.get(meta_server, "type")
                    meta_client = connect_to_server(meta_server,
                                                    config, clients)
                    if meta_server_type in ('fdsn', 'jane'):
                        inventory = meta_client.get_stations(
                            network=net, station=sta, location=loc,
                            level="response")
                        # look in all metadata available, but prefer metadata from
                        # same source (by putting it in front in list)
                        all_inventories = [inventory] + all_inventories
                        _attach_metadata(st, all_inventories)
                    else:
                        raise NotImplementedError()
            sta_fetched.add(net_sta_loc)
            sys.stdout.write("\r%s (%s: %s) fetched.\n" % (
                seed_id.ljust(15), server_type, server))
            sys.stdout.flush()
        except Exception as e:
            sys.stdout.write(
                "\r%s (%s: %s) skipped! (Exception: %s)\n" % (
                    seed_id.ljust(15), server_type, server, e))
            sys.stdout.flush()
            continue
        # check whether to attempt rotation
        if not no_metadata:
            if config.has_section("rotate_channels"):
                if net_sta_loc in config.options("rotate_channels"):
                    rotate_channels(st, net, sta, loc, config)
        # SeisHub
        if server_type == "seishub":
            for tr in st:
                if tr.stats._format == 'GSE2':
                    apply_gse2_calib(tr)
                tr.stats['_format'] = "SeisHub"
        # ArcLink
        elif server_type == "arclink":
            for tr in st:
                tr.stats['_format'] = "ArcLink"
        # FDSN (or JANE)
        elif server_type in ("fdsn", "jane"):
            for tr in st:
                tr.stats['_format'] = "FDSN"
        # seedlink
        elif server_type == "seedlink":
            for tr in st:
                tr.stats['_format'] = "SeedLink"
        # SDS
        elif server_type == "sds":
            for tr in st:
                tr.stats['_format'] = "SDS"
        streams.append(st)
    print("=" * 80)
    return (clients, streams, all_inventories)


def rotate_channels(st, net, sta, loc, config):
    net_sta_loc = ".".join((net, sta, loc))
    channels = config.get("rotate_channels", net_sta_loc).split(",")
    tr = st.select(id=".".join((net_sta_loc, channels[0])))[0]
    parser = tr.stats.get("parser")
    coordinates = tr.stats.get("coordinates")
    response = tr.stats.get("response")
    st = _rotate_specific_channels_to_zne(
        st, net, sta, loc, channels)
    for tr in st:
        if parser is not None:
            tr.stats.parser = copy.deepcopy(parser)
        if coordinates is not None:
            tr.stats.coordinates = copy.deepcopy(coordinates)
        if response is not None:
            tr.stats.response = copy.deepcopy(response)


def connect_to_server(server_name, config, clients):
    """
    Return existing client for given server name or connect to server, add it
    to the clients dictionary and return client.

    :type server_name: str
    :param server_name: Server name as specified by configuration file
    :type config: :class:`ConfigParser.ConfigParser`
    :param config: ObsPyck configuration.
    :type clients: dict
    :param clients: Dictionary with already connected servers, keys are server
        names as specified in configuration, values are corresponding client
        instances.
    """
    if server_name in clients:
        return clients[server_name]

    server_type = config.get(server_name, "type")
    # # doesnt work on obspy <1.0, so set it above on module level:
    # ArcLinkClient.max_status_requests = 2000

    client_classes = {
        "arclink": ArcLinkClient,
        "fdsn": FDSNClient,
        "jane": FDSNClient,
        "seishub": SeisHubClient,
        "seedlink": SeedlinkClient,
        "sds": SDSClient,
        }

    if server_type not in client_classes.keys():
        msg = ("Unknown server type '{}' in server definition section "
               "'{}' in config file.").format(server_type, server_name)
        raise NotImplementedError(msg)

    config_keys = {
        "arclink": (
            "host", "port", "user", "password", "institution", "timeout",
            "dcid_key_file", "debug", "command_delay", "status_delay"),
        "fdsn": (
            "base_url", "user", "password", "user_agent", "debug", "timeout"),
        "jane": (
            "base_url", "user", "password", "user_agent", "debug", "timeout"),
        "seishub": (
            "base_url", "user", "password", "timeout", "debug", "retries"),
        "seedlink": (
            "server", "port", "timeout", "debug"),
        "sds": (
            "sds_root", "sds_type", "format", "fileborder_seconds",
            "fileborder_samples"),
        }
    config_getters = {
        "port": config.getint, "timeout": config.getfloat,
        "debug": config.getboolean, "command_delay": config.getfloat,
        "status_delay": config.getfloat, "retries": config.getint,
        "fileborder_seconds": config.getfloat,
        "fileborder_samples": config.getint,
        }

    kwargs = {}
    for key in config_keys[server_type]:
        value = config_getters.get(key, config.get)(server_name, key) or None
        if value is not None:
            kwargs[key] = value

    client = client_classes[server_type](**kwargs)
    clients[server_name] = client
    return client


def merge_check_and_cleanup_streams(streams, options, config):
    """
    Cleanup given list of streams so that they conform with what ObsPyck
    expects.

    Conditions:
    - either one Z or three ZNE traces
    - no two streams for any station (of same network)
    - no streams with traces of different stations

    :returns: (warn_msg, merge_msg, list(:class:`obspy.core.stream.Stream`s))
    """
    # we need to go through streams/dicts backwards in order not to get
    # problems because of the pop() statement
    warn_msg = ""
    merge_msg = ""
    # Merge on every stream if this option is passed on command line:
    for st in streams:
        st.merge(method=-1)
    merge_type = config.get("base", "merge")
    if merge_type:
        if merge_type == "safe":
            for st in streams:
                st.merge(method=0)
        elif merge_type == "overwrite":
            for st in streams:
                if st.get_gaps() and max([gap[-1] for gap in st.get_gaps()]) < 5:
                    msg = 'Interpolated over gap(s) with less than 5 ' + \
                          'samples for station: %s.%s'
                    msg = msg % (st[0].stats.network, st[0].stats.station)
                    warn_msg += msg + "\n"
                    st.merge(method=1, fill_value="interpolate")
                else:
                    st.merge(method=1)
        else:
            err = "Unrecognized config value for merging traces. Try " + \
                  "'safe' or 'overwrite' (or leave empty)."
            raise Exception(err)

    # Sort streams again, if there was a merge this could be necessary 
    for st in streams:
        st.traces = sorted(st.traces, key=_seed_id_keyfunction)
    sta_list = set()
    # XXX we need the list() because otherwise the iterator gets garbled if
    # XXX removing streams inside the for loop!!
    for st in list(streams):
        # check for streams with mixed stations/networks and remove them
        if len(st) != len(st.select(network=st[0].stats.network,
                                    station=st[0].stats.station)):
            msg = "Warning: Stream with a mix of stations/networks. " + \
                  "Discarding stream:\n%s" % str(st)
            print(msg)
            warn_msg += msg + "\n"
            streams.remove(st)
            continue
        net_sta = "%s.%s" % (st[0].stats.network.strip(),
                             st[0].stats.station.strip())
        # Here we make sure that a station/network combination is not
        # present with two streams.
        if net_sta in sta_list:
            msg = ("Warning: Station/Network combination '%s' "
                   "already in stream list. Discarding stream.") % net_sta
            print(msg)
            warn_msg += msg + "\n"
            streams.remove(st)
            continue
        sta_list.add(net_sta)
    # demean traces if not explicitly deactivated on command line
    if config.getboolean("base", "zero_mean"):
        for st in streams:
            try:
                st.detrend('simple')
                st.detrend('constant')
            except NotImplementedError as e:
                if "Trace with masked values found." in e.message:
                    msg = 'Detrending/demeaning not possible for station ' + \
                          '(masked Traces): %s' % net_sta
                    warn_msg += msg + "\n"
                else:
                    raise
    return (warn_msg, merge_msg, streams)


def cleanup_streams_without_metadata(streams):
    """
    Function to remove streams that do not provide the necessary metadata.

    :returns: list of streams
    """
    # we need to go through streams/dicts backwards in order not to get
    # problems because of the pop() statement
    for i in range(len(streams))[::-1]:
        ok = True
        for tr in streams[i]:
            for key in ('orientation', 'coordinates'):
                if key not in tr.stats:
                    ok = False
            # for dataless we have parser, for stationxml/fdsnws we have response
            if "parser" not in tr.stats and "response" not in tr.stats:
                ok = False
            if not ok:
                print(f'Error: Missing metadata for "{tr.id}". Discarding '
                      'stream.')
                streams.pop(i)
                break
    return streams


def setup_external_programs(options, config):
    """
    Sets up temdir, copies program files, fills in PROGRAMS dict, sets up
    system calls for programs.
    Depends on command line options, returns temporary directory.

    :param options: Command line options of ObsPyck
    :type options: options as returned by :meth:`optparse.OptionParser.parse_args`
    :returns: String representation of temporary directory with program files.
    """
    pluginpath = config.get("base", "pluginpath") or None
    if pluginpath is None:
        pluginpath = os.path.dirname(os.path.abspath(__file__))
    tmp_dir = tempfile.mkdtemp(prefix="obspyck-")
    if not os.path.isdir(pluginpath):
        msg = ("'pluginpath' in config file is no valid directory: '%s'\n"
               "Localization methods/functions are deactivated.") % pluginpath
        warnings.warn(msg)
        return tmp_dir
    # set binary names to use depending on architecture and platform...
    env = os.environ
    architecture = platform.architecture()[0]
    system = platform.system()
    global SHELL
    if system == "Windows":
        SHELL = True
    else:
        SHELL = False
    # link velocity models ################################################
    if os.path.exists(os.path.join(pluginpath, "VELOCITY_MODELS")):
        os.symlink(os.path.join(pluginpath, "VELOCITY_MODELS"),
                   os.path.join(tmp_dir, "VELOCITY_MODELS"))
    # Setup external programs #############################################
    for prog_basename, prog_dict in PROGRAMS.items():
        prog_srcpath = os.path.join(pluginpath, prog_basename)
        prog_tmpdir = os.path.join(tmp_dir, prog_basename)
        prog_dict['dir'] = prog_tmpdir
        shutil.copytree(prog_srcpath, prog_tmpdir, symlinks=True)
        prog_dict['files'] = {}
        for key, filename in prog_dict['filenames'].items():
            prog_dict['files'][key] = os.path.join(prog_tmpdir, filename)
        prog_dict['files']['exe'] = "__".join(\
                [prog_dict['filenames']['exe'], system, architecture])
        # setup clean environment
        prog_dict['env'] = {}
        prog_dict['env']['PATH'] = prog_dict['dir'] + os.pathsep + env['PATH']
        if 'SystemRoot' in env:
            prog_dict['env']['SystemRoot'] = env['SystemRoot']
    # Hyp2000 #############################################################
    prog_dict = PROGRAMS['hyp_2000']
    prog_dict['env']['HYP2000_DATA'] = prog_dict['dir'] + os.sep
    def tmp(prog_dict):
        files = prog_dict['files']
        for file in [files['phases'], files['stations'], files['summary']]:
            if os.path.isfile(file):
                os.remove(file)
        return
    prog_dict['PreCall'] = tmp
    def tmp(prog_dict):
        sub = subprocess.Popen(prog_dict['files']['exe'], shell=SHELL,
                cwd=prog_dict['dir'], env=prog_dict['env'],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        input = open(prog_dict['files']['control'], "rb").read()
        (msg, err) = sub.communicate(input)
        if system == "Darwin":
            returncode = sub.returncode
        else:
            returncode = sub.wait()
        return (msg, err, returncode)
    prog_dict['Call'] = tmp
    # NLLoc ###############################################################
    prog_dict = PROGRAMS['nlloc']
    def tmp(prog_dict):
        filepattern = os.path.join(prog_dict['dir'], "nlloc*")
        print(filepattern)
        for file in glob.glob(filepattern):
            os.remove(file)
        return
    prog_dict['PreCall'] = tmp
    def tmp(prog_dict, controlfilename):
        sub = subprocess.Popen([prog_dict['files']['exe'], controlfilename],
                cwd=prog_dict['dir'], env=prog_dict['env'], shell=SHELL,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if system == "Darwin":
            returncode = sub.returncode
        else:
            returncode = sub.wait()
        msg = sub.stdout.read().decode('UTF-8')
        err = sub.stderr.read().decode('UTF-8')
        for pattern, key in [("nlloc.*.*.*.loc.scat", 'scatter'),
                             ("nlloc.*.*.*.loc.hyp", 'summary')]:
            pattern = os.path.join(prog_dict['dir'], pattern)
            newname = os.path.join(prog_dict['dir'], prog_dict['files'][key])
            for file in glob.glob(pattern):
                os.rename(file, newname)
        return (msg, err, returncode)
    prog_dict['Call'] = tmp
    # focmec ##############################################################
    prog_dict = PROGRAMS['focmec']
    def tmp(prog_dict):
        sub = subprocess.Popen(prog_dict['files']['exe'], shell=SHELL,
                cwd=prog_dict['dir'], env=prog_dict['env'],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if system == "Darwin":
            returncode = sub.returncode
        else:
            returncode = sub.wait()
        msg = sub.stdout.read().decode('UTF-8')
        err = sub.stderr.read().decode('UTF-8')
        return (msg, err, returncode)
    prog_dict['Call'] = tmp
    #######################################################################
    return tmp_dir

#Monkey patch (need to remember the ids of the mpl_connect-statements to remove them later)
#See source: http://matplotlib.sourcearchive.com/documentation/0.98.1/widgets_8py-source.html
class MultiCursor(MplMultiCursor):
    def __init__(self, canvas, axes, useblit=True, **lineprops):
        if hasattr(self, "id1"):
            self.canvas.mpl_disconnect(self.id1)
            self.canvas.mpl_disconnect(self.id2)
        self.canvas = canvas
        self.axes = axes
        xmin, xmax = axes[-1].get_xlim()
        xmid = 0.5*(xmin+xmax)
        self.hlines = []
        self.vlines = []
        self.horizOn = False
        self.vertOn = True
        self.lines = [ax.axvline(xmid, visible=False, **lineprops) for ax in axes]
        self.visible = True
        self.useblit = useblit
        self.background = None
        self.needclear = False
        self.id1 = self.canvas.mpl_connect('motion_notify_event', self.onmove)
        self.id2 = self.canvas.mpl_connect('draw_event', self.clear)

    @property
    def lines(self):
        if MATPLOTLIB_VERSION < [1, 3, 0]:
            return self.__dict__["lines"]
        else:
            return self.vlines

    @lines.setter
    def lines(self, value):
        if MATPLOTLIB_VERSION < [1, 3, 0]:
            self.__dict__["lines"] = value
        else:
            self.vlines = value


def gk2lonlat(x, y, m_to_km=True):
    """
    This function converts X/Y Gauss-Krueger coordinates (zone 4, central
    meridian 12 deg) to Longitude/Latitude in WGS84 reference ellipsoid.
    We do this using pyproj (python bindings for proj4) which can be installed
    using 'easy_install pyproj' from pypi.python.org.
    Input can be single coordinates or coordinate lists/arrays.
    
    Useful Links:
    http://pyproj.googlecode.com/svn/trunk/README.html
    http://trac.osgeo.org/proj/
    http://www.epsg-registry.org/
    """
    import pyproj

    # convert to meters first
    if m_to_km:
        x = x * 1000.
        y = y * 1000.
    # pyproj has deprecated the old 'init="epsg:4326"' syntax
    try:
        transformer = pyproj.Transformer.from_crs(31468, 4326, always_xy=True)
        lon, lat = transformer.transform(x, y)
    except Exception:
        proj_wgs84 = pyproj.Proj(init="epsg:4326")
        proj_gk4 = pyproj.Proj(init="epsg:31468")
        lon, lat = pyproj.transform(proj_gk4, proj_wgs84, x, y)
    return (lon, lat)


def readNLLocScatter(scat_filename, textviewStdErrImproved):
    """
    This function reads location and values of pdf scatter samples from the
    specified NLLoc *.scat binary file (type "<f4", 4 header values, then 4
    floats per sample: x, y, z, pdf value) and converts X/Y Gauss-Krueger
    coordinates (zone 4, central meridian 12 deg) to Longitude/Latitude in
    WGS84 reference ellipsoid.
    Messages on stderr are written to specified GUI textview.
    Returns an array of xy pairs.
    """
    # read data, omit the first 4 values (header information) and reshape
    data = np.fromfile(scat_filename, dtype="<f4").astype("float")[4:]
    data = data.reshape((-1, 4)).swapaxes(0, 1)
    data[0], data[1] = gk2lonlat(data[0], data[1])
    return data.T


def errorEllipsoid2CartesianErrors(azimuth1, dip1, len1, azimuth2, dip2, len2,
                                   len3):
    """
    This method converts the location error of NLLoc given as the 3D error
    ellipsoid (two azimuths, two dips and three axis lengths) to a cartesian
    representation.
    We calculate the cartesian representation of each of the ellipsoids three
    eigenvectors and use the maximum of these vectors components on every axis.
    """
    z = len1 * np.sin(np.radians(dip1))
    xy = len1 * np.cos(np.radians(dip1))
    x = xy * np.sin(np.radians(azimuth1))
    y = xy * np.cos(np.radians(azimuth1))
    v1 = np.array([x, y, z])

    z = len2 * np.sin(np.radians(dip2))
    xy = len2 * np.cos(np.radians(dip2))
    x = xy * np.sin(np.radians(azimuth2))
    y = xy * np.cos(np.radians(azimuth2))
    v2 = np.array([x, y, z])

    v3 = np.cross(v1, v2)
    v3 /= np.sqrt(np.dot(v3, v3))
    v3 *= len3

    v1 = np.abs(v1)
    v2 = np.abs(v2)
    v3 = np.abs(v3)

    error_x = max([v1[0], v2[0], v3[0]])
    error_y = max([v1[1], v2[1], v3[1]])
    error_z = max([v1[2], v2[2], v3[2]])
    
    return (error_x, error_y, error_z)

def formatXTicklabels(x, *pos):
    """
    Make a nice formatting for y axis ticklabels: minutes:seconds.microsec
    """
    # x is of type numpy.float64, the string representation of that float
    # strips of all tailing zeros
    # pos returns the position of x on the axis while zooming, None otherwise

    if x == 0:
        return '0s'

    info = []

    if x < 0:
        info.append('-')
    x = abs(x)

    minutes = int(x / 60.)
    seconds = x % 60

    if minutes > 0:
        info.append('%imin' % minutes)

    if seconds != 0:
        info.append("%ss" % seconds)

    return ' '.join(info)


def coords2azbazinc(stream, origin):
    """
    Returns azimuth, backazimuth and incidence angle from station coordinates
    given in first trace of stream and from event location specified in origin
    dictionary.
    """
    sta_coords = stream[0].stats.coordinates
    dist, bazim, azim = gps2dist_azimuth(sta_coords.latitude,
            sta_coords.longitude, float(origin.latitude),
            float(origin.longitude))
    elev_diff = sta_coords.elevation - float(origin.depth)
    inci = math.atan2(dist, elev_diff) * 180.0 / math.pi
    return azim, bazim, inci

class SplitWriter():
    """
    Implements a write method that writes a given message on all children
    """
    def __init__(self, *objects):
        """
        Remember provided objects as children.
        """
        self.children = objects

    def write(self, msg):
        """
        Sends msg to all childrens write method.
        """
        for obj in self.children:
            if isinstance(obj, QtWidgets.QPlainTextEdit):
                if msg == '\n':
                    return
                if msg.endswith('\n'):
                    msg_ = msg[:-1]
                else:
                    msg_ = msg
                obj.appendPlainText(msg_)
            else:
                obj.write(msg)

    def flush(self, *args, **kwargs):
        """
        Delegates `flush()` calls to children if applicable
        """
        for obj in self.children:
            try:
                obj.flush(*args, **kwargs)
            except AttributeError:
                pass


def getArrivalForPick(arrivals, pick):
    """
    searches given arrivals for an arrival that references the given
    pick and returns it (empty Arrival object otherwise).
    """
    for a in arrivals:
        if a.pick_id == pick.resource_id:
            return a
    return None


def getPickForArrival(picks, arrival):
    """
    searches list of picks for a pick that matches the arrivals pick_id
    and returns it (empty Pick object otherwise).
    """
    pick = None
    for p in picks:
        if arrival.pick_id == p.resource_id:
            pick = p
            break
    return pick


def get_event_info(starttime, endtime, streams):
    events = []
    arrivals = {}
    try:
        msg = ('getTravelTimes not available in obspy.taup in newer obspy '
               'versions anymore, so this functionality is currently not '
               'implemented.')
        raise NotImplementedError(msg)
        # client = FDSNClient("NERIES")
        # events = client.get_events(starttime=starttime - 20 * 60,
        #                            endtime=endtime)
        # for ev in events[::-1]:
        #     has_arrivals = False
        #     origin = ev.origins[0]
        #     origin_time = origin.time
        #     lon1 = origin.longitude
        #     lat1 = origin.latitude
        #     depth = abs(origin.depth / 1e3)
        #     for st in streams:
        #         sta = st[0].stats.station
        #         lon2 = st[0].stats.coordinates['longitude']
        #         lat2 = st[0].stats.coordinates['latitude']
        #         dist = locations2degrees(lat1, lon1, lat2, lon2)
        #         tts = getTravelTimes(dist, depth)
        #         list_ = arrivals.setdefault(sta, [])
        #         for tt in tts:
        #             tt['time'] = origin_time + tt['time']
        #             if starttime < tt['time'] < endtime:
        #                 has_arrivals = True
        #                 list_.append(tt)
        #     if not has_arrivals:
        #         events[:] = events[:-1]
    except Exception as e:
        msg = ("Problem while fetching events or determining theoretical "
               "phases: %s: %s" % (e.__class__.__name__, str(e)))
        return None, None, msg
    return events, arrivals, None


def apply_gse2_calib(tr):
    """
    Applies GSE2 specific calibration to overall sensitivity.
    Not valid for accelerometer data!
    """
    try:
        calibration = tr.stats.calib * ((2.0 * np.pi / tr.stats.gse2.calper) ** 1) * 1e-9
        # tr.stats.paz.sensitivity = tr.stats.paz.sensitivity / calibration
        # XXX this is not implemented anymore after switching to use parser in
        # Trace.simulate() with evalresp
        raise NotImplementedError()
    except Exception as e:
        msg = ("Warning: Failed to apply GSE2 calibration factor to overall "
               "sensitivity (%s, %s). Continuing anyway.")
        msg = msg % (e.__class__.__name__, str(e))
        print(msg)


def map_rotated_channel_code(channel, rotation):
    """
    Modifies a channel code according to given rotation (e.g. EHN -> EHR)

    :type channel: str
    :type rotation: str
    """
    if rotation in ("LQT", "ZRT"):
        while len(channel) < 3:
            msg = ("Channel code ('%s') does not have three characters. "
                   "Filling with leading spaces.")
            warnings.warn(msg % channel)
            channel = " " + channel
        if rotation == "LQT":
            mapping = ROTATE_LQT_COMP_MAP
        elif rotation == "ZRT":
            mapping = ROTATE_ZRT_COMP_MAP
        channel = channel[0:2] + mapping[channel[2]]
    elif rotation is None:
        pass
    return channel


# copied from obspy master pre 1.1.0
def _seed_id_keyfunction(tr):
    """
    Keyfunction to use in sorting two Traces by SEED ID

    Assumes that the last (or only) "."-separated part is a channel code.
    Assumes the last character is a the component code and sorts it
    "Z"-"N"-"E"-others_lexical.
    """
    x = tr.id
    # for comparison we build a list of 5 SEED code pieces:
    # [network, station, location, band+instrument, component]
    # with partial codes (i.e. not 4 fields after splitting at dots),
    # we go with the following assumptions (these seem a bit random, but that's
    # what can be encountered in string representations of the Inventory object
    # hierarchy):
    #  - no dot means network code only (e.g. "IU")
    #  - one dot means network.station code only (e.g. "IU.ANMO")
    #  - two dots means station.location.channel code only (e.g. "ANMO.10.BHZ")
    #  - three dots: full SEED ID (e.g. "IU.ANMO.10.BHZ")
    #  - more dots: sort after any of the previous, plain lexical sort
    # if no "." in the string: assume it's a network code

    # split to get rid of the description that that is added to networks and
    # stations which might also contain dots.
    number_of_dots = x.strip().split()[0].count(".")

    x = x.upper()
    if number_of_dots == 0:
        x = [x] + [""] * 4
    elif number_of_dots == 1:
        x = x.split(".") + [""] * 3
    elif number_of_dots in (2, 3):
        x = x.split(".")
        if number_of_dots == 2:
            x = [""] + x
        # split channel code into band+instrument code and component code
        x = x[:-1] + [x[-1][:-1], x[-1] and x[-1][-1] or '']
        # special comparison for component code, convert "ZNE" to integers
        # which compare less than any character
        comp = "ZNE".find(x[-1])
        # last item is component code, either the original 1-char string, or an
        # int from 0-2 if any of "ZNE". Python3 does not allow comparison of
        # int and string anymore (Python 2 always compares ints smaller than
        # any string), so we need to work around this by making this last item
        # a tuple with first item False for ints and True for strings.
        if comp >= 0:
            x[-1] = (False, comp)
        else:
            x[-1] = (True, x[-1])
    # all other cases, just convert the upper case string to a single item
    # list - it will compare greater than any of the split lists.
    else:
        x = [x, ]

    return x


def set_matplotlib_defaults(config):
    """
    Sets matplotlib rc customizations specified in config file, section
    ``[matplotlibrc]``.

    :type config: :class:`ConfigParser.SafeConfigParser`
    """
    if config.has_section('matplotlibrc'):
        for key, value in config.items('matplotlibrc'):
            mpl.rcParams[key] = value


def _save_input_data(streams, inventories, directory):
    """
    :type streams: list of Stream
    :type inventories: list of Inventory
    :type directory: str
    """
    if not os.path.isdir(directory):
        raise OSError('Not a directory: ' + str(directory))
    if isinstance(streams, (list, tuple)):
        pass
    elif isinstance(streams, Stream):
        streams = [streams]
    else:
        raise TypeError
    # write raw waveforms
    st = Stream()
    for st_ in streams:
        st += st_
    with warnings.catch_warnings(record=True):
        st.write(os.path.join(directory, 'waveforms.mseed'), format='MSEED')
    # write inventories
    inv = Inventory(networks=[], source='')
    for inv_ in inventories:
        inv += inv_
    inv.write(os.path.join(directory, 'inventory.xml'), format='STATIONXML')
