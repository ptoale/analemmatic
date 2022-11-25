#!/usr/bin/env python
"""
solartime: module to help with handling local solar time

Solar time is based on the position of the sun in the local sky. There are two
kinds of solar time:
    1) Apparent Solar Time: this is the time told by a sundial
    2) Mean Solar Time: this is the local time told by a watch (kinda)

Apparent solar noon is defined by the Sun's transit of the local meridian. For a
location in the northern hemisphere, the Sun will be due south when this occurs.
This will typically not correspond to mean solar noon. One hour of apparent 
solar time is the time it takes for the Sun to move 15 degrees through the sky.
One apparent solar day is therefore the time it takes the Sun to travel 360
degrees and return to the local meridian. This amount of time is not uniform 
throughout the year due to the eccentricity of the Earth's orbit and to the
tilt of its axis. The equation of time describes the difference between these 
two time scales.

Mean solar time is therefore defined in terms of a mean Sun, that travels along
the celestial equator at a constant rate equal to the mean rate of the true Sun. 
Even the length of a mean solar day slowly changes (over 100 years). The current 
length of the mean solar day is 86400.002 seconds. The use of time zones means that
mean solar time is constant over regions approximately 15 degrees of longitude
wide. Mean solar noon closely corresponds to apparent solar noon only on the 
standard meridian of the zone. In the western hemisphere, this is the easterly
side of the zone. 

Daylight savings time (DST) in the US currently runs from 2:00am local time on the
second Sunday of March to 2:00am local time on the first Sunday of November, covering
34 weeks of the year. 

The Sun is on the local meridian at noon on a certain day of the year. What UTC time
does this correspond to?
    1) Apparent solar time is noon, which means AST=0. Start by adding 12
    2) Correct for location within the time zone:
        dt = (4 min/deg)*(longitude - LSTM)
       where LSTM can be found from: LSTM = (15 deg/hour)*(delta_GMT)
    3) Make the equation of time correction. This depends on the day of the year
       and must be computed.
    4) Adjust for daylight savings time, if necessary. This requires adding an
       hour during the summer-half of the year. This is now equal to local
       standard time (LST)
    5) Convert to UTC using the time zone offset

Routines to convert between local standard time and UTC:
    1) get_tzinfo: interface to pytz tzinfo objects
    2) make_aware: add tzinfo objects to naive datetime
    3) local_to_utc: convert from local standard time to utc
    4) utc_to_local: convert from utc to local standard time
    
Routines to convert between local standard time, mean solar time, and (apparent) solar time:
    1) local_standard_time_meridian: find the meridian of the local time zone
    2) longitude_correction: correct for location within time zone
    3) eq_of_time: the equation of time
    4) dst: account for daylight savings time
    5) local_to_solar: convert from local standard time to apparent solar time
    6) solar_to_local: convert from apparent solar time to local standard time
    7) local_to_mean_solar: convert from local standard time to mean solar time
    8) mean_solar_to_local: convert from mean solar time to local standard time
    9) hour_angle: convert from local solar time to hour angle
    
Routines to convert between mean solar time and utc:
    1) mean_solar_to_utc: convert from mean solar time to utc
    2) utc_to_mean_solar: convert from utc to mean solar time

"""
import math
from datetime import datetime, timedelta
import pytz


def get_tzinfo(zone):
    """
    Get a tzinfo object for a given time zone.
    
    Args:
        zone (str): the time zone official name
        
    Returns:
        A tzinfo object, or None

    Try getting a tzinfo object for US Central time
    >>> tz = get_tzinfo('US/Central')
    >>> print(tz.zone)
    US/Central
        
    What happens if we mess it up?
    >>> tz = get_tzinfo('US/Centrol')
    ERROR: Unknown time zone: US/Centrol
    >>> if tz: print(tz.zone)
            
    """

    "Create the tzinfo object, or return None"
    tz = None
    try:
        tz = pytz.timezone(zone)
    except pytz.UnknownTimeZoneError:
        print('ERROR: Unknown time zone: {}'.format(zone))
    return tz


def make_aware(dt, zone):
    """
    Make a datetime object aware of a timezone. If it's a naive datetime,
    it will be localized to the timezone. If it already has a tzinfo, it will
    be translated to the timezone.
    
    Args:
        dt (datetime): the naive datetime object
        zone (str): the local timezone official name
        
    Returns:
        The input datetime object, but now with timezone awareness. If the
        zone is not a correct timezone name, None is returned.
    
    Tell it its in CDT    
    >>> dt_CDT = make_aware(datetime(2016, 4, 15, 13, 43), 'US/Central')
    >>> print(dt_CDT)
    2016-04-15 13:43:00-05:00

    Wait, I meant EDT...
    >>> dt_EDT = make_aware(dt_CDT, 'US/Eastern')
    >>> print(dt_EDT)
    2016-04-15 14:43:00-04:00

    Or a madeup timezone
    >>> dt_CDT = make_aware(datetime(2016, 4, 15, 13, 43), 'US/Somewhere')
    ERROR: Unknown time zone: US/Somewhere
    >>> print(dt_CDT)
    None
    
    """

    "Create the tz object"
    tz = get_tzinfo(zone)
    if not tz:
        return None

    "If the datetime object already has a tzinfo..."
    if dt.tzinfo:
        return dt.astimezone(tz)
        
    "Otherwise, try to localize it"
    try:
        return tz.localize(dt, is_dst=None)
    except pytz.NonExistentTimeError:
        print('ERROR: There is no UTC time for {} in {}'.format(dt, tz.zone))
        return None
    except pytz.AmbiguousTimeError:
        print('ERROR: There are multiple UTC times for {} in {}'.format(dt, tz.zone))
        return None


def local_to_utc(dt, zone):
    """
    Convert a datetime object from local standard time to UTC.
    
    Args:
        dt (datetime): the datetime object, assumed to be naive
        zone (str): the local timezone official name
        
    Returns:
        The datetime object adjusted to UTC, 
        or None if there is no corresponding UTC (usually related to DST)
    
    There should be a 6-hour offset between CST and UTC
    >>> xmas_2016_local = datetime(2016, 12, 25, 0, 0, 0, 0)
    >>> xmas_2016_utc = local_to_utc(xmas_2016_local, 'US/Central')
    >>> print(xmas_2016_utc)
    2016-12-25 06:00:00+00:00
    
    And a 5-hour offset between CDT and UTC
    >>> j4_2016_local = datetime(2016, 7, 4, 0, 0, 0, 0)
    >>> j4_2016_utc = local_to_utc(j4_2016_local, 'US/Central')
    >>> print(j4_2016_utc)
    2016-07-04 05:00:00+00:00
    
    There are local times that don't correspond to any UTC time    
    >>> dst_2016_local = datetime(2016, 3, 13, 2, 30, 0, 0)
    >>> dst_2016_utc = local_to_utc(dst_2016_local, 'US/Central')
    ERROR: There is no UTC time for 2016-03-13 02:30:00 in US/Central
    >>> print(dst_2016_utc)
    None

    And there can be multiple local times that correspond to the same UTC time
    >>> dst_2016_local = datetime(2016, 11, 6, 1, 30, 0, 0)
    >>> dst_2016_utc = local_to_utc(dst_2016_local, 'US/Central')
    ERROR: There are multiple UTC times for 2016-11-06 01:30:00 in US/Central
    >>> print(dst_2016_utc)
    None

    """

    "Make the dt timezone aware, if possible"
    local_dt = make_aware(dt, zone)
    if not local_dt:
        return None
    
    "Now convert to UTC"
    utc_dt = local_dt.astimezone(pytz.utc)
    return utc_dt


def utc_to_local(dt, zone):
    """
    Convert a datetime object from UTC to local standard time.
    
    Args:
        dt (datetime): the datetime object, assumed to be naive
        zone (str): the local timezone official name
        
    Returns:
        The datetime object adjusted to local time, or None 
    
    There should be a 6-hour offset between CST and UTC
    >>> xmas_2016_utc = datetime(2016, 12, 25, 0, 0, 0, 0)
    >>> xmas_2016_local = utc_to_local(xmas_2016_utc, 'US/Central')
    >>> print(xmas_2016_local)
    2016-12-24 18:00:00-06:00
    
    And a 5-hour offset between CDT and UTC
    >>> j4_2016_utc = datetime(2016, 7, 4, 0, 0, 0, 0)
    >>> j4_2016_local = utc_to_local(j4_2016_utc, 'US/Central')
    >>> print(j4_2016_local)
    2016-07-03 19:00:00-05:00

    """

    "Create the tz object"
    tz = get_tzinfo(zone)
    if not tz:
        return None

    "Make the dt timezone aware, if possible"
    utc_dt = make_aware(dt, 'UTC')
    if not utc_dt:
        return None
    
    "Now convert to UTC"
    local_dt = utc_dt.astimezone(tz)
    return local_dt


def local_standard_time_meridian(zone):
    """
    Calculate the local standard time meridian for a zone. 
    Be careful with DST!
    
    Args:
        zone (str): the local timezone official name
        
    Returns:
        The longitude, in degrees, of the standard meridian for this zone.
    
    >>> print(local_standard_time_meridian('US/Central'))
    -90.0
    >>> print(local_standard_time_meridian('US/Pacific'))
    -120.0
    
    """
    
    "Get the timezone"
    tz = get_tzinfo(zone)
    if not tz:
        return None
    
    "Pick a time that is not effected by DST and localize it"
    dt = datetime(2015, 12, 25, 0, 0, 0)
    local = make_aware(dt, zone)
    if not local:
        return None
    
    "Find the longitude from the hour offset from utc"
    lon = local.utcoffset().total_seconds()*15./3600.
    
    return lon


def dst(dt, zone):
    """
    Returns the number of minutes of any DST offset.
    
    Args:
        dt (datetime): the time to check
        zone (str): the local timezone official name
        
    Returns:
        The number of minutes to adjust the datetime by due to DST

    No DST in Febuary    
    >>> print(dst(datetime(2016, 2, 1, 12, 0, 0), 'US/Central'))
    0.0

    DST in June
    >>> print(dst(datetime(2016, 6, 1, 12, 0, 0), 'US/Central'))
    60.0

    No DST in December
    >>> print(dst(datetime(2016, 12, 1, 12, 0, 0), 'US/Central'))
    0.0
    
    """
    
    "Get the timezone"
    tz = get_tzinfo(zone)
    if not tz:
        return None

    "Localize the datetime"
    local = make_aware(dt, zone)
    if not local:
        return None
    
    "Return the DST time offset, in minutes"
    offset = local.dst().total_seconds()/60.
    return offset


def equation_of_time(dt):
    """
    Simple implementation of the equation of time.
    
    Args:
        dt (datetime): the datetime object used to evaluate the equation of time
        
    Returns:
        The number of minutes (as a float) between the length of the standard
        day and the length of the solar day, on the given date
    
    >>> print(equation_of_time(datetime(2016, 1, 1, 12, 0, 0)))
    -3.7286600067084477
    >>> print(equation_of_time(datetime(2016, 2, 13, 12, 0, 0)))
    -14.599483104066074
    >>> print(equation_of_time(datetime(2016, 4, 15, 12, 0, 0)))
    0.011148309685002133
    >>> print(equation_of_time(datetime(2016, 10, 30, 12, 0, 0)))
    16.45318866131135
    
    """
    
    "Get the day of the year"
    d = dt.timetuple().tm_yday

    "The argument to the trig functions..."
    arg = math.radians(360.*(d-81)/365.24)

    "And the EoT itself"
    return 9.87*math.sin(2*arg) - 7.53*math.cos(arg) - 1.5*math.sin(arg)


def longitude_correction(zone, lon):
    """
    Calculate the difference between mean solar time and local standard time
    based on location within time zone (longitude).

    Args:
        zone (str): the local timezone official name
        lon (float): the longitude of the location, in decimal degrees
        
    Returns:
        The number of minutes beteen mean solar time and local standard time,
        due to location within time zone.

    Longitude correction is a constant throughout the year    
    >>> print(longitude_correction('US/Central', -87.6086066062))
    9.565573575200006
    >>> print(longitude_correction('US/Central', -75.))
    60.0
    >>> print(longitude_correction('US/Central', -90.))
    0.0
    >>> print(longitude_correction('US/Central', -105.))
    -60.0

    """

    "Get the LSTM for the local time zone"
    lstm = local_standard_time_meridian(zone)
    if not lstm:
        return None
    
    "Calculate the time difference, in minutes, between the LSTM and the location."
    return 4.*(lon - lstm)


def local_to_mean_solar(dt, zone, lon):
    """
    Convert a datetime object from local standard time to mean solar time.

    Args:
        dt (datetime): the local standard time
        zone (str): the local timezone official name
        lon (float): the longitude of the location, in decimal degrees

    Returns:
        The mean solar time

    >>> print(local_to_mean_solar(datetime(2016, 1, 1, 12, 0, 0), 'US/Central', -87.6086066062))
    2016-01-01 12:09:33.934415
    """

    "Calculate the time difference, in minutes, between the LSTM and the location."
    delta1 = longitude_correction(zone, lon)

    "Adjust for DST"
    delta2 = -dst(dt, zone)

    "Convert the local time"
    corr = delta1 + delta2
    dt_solar = dt + timedelta(minutes=corr)

    return dt_solar


def local_to_solar(dt, zone, lon):
    """
    Convert a datetime object from local standard time to solar time.

    Args:
        dt (datetime): the local standard time
        zone (str): the local timezone official name
        lon (float): the longitude of the location, in decimal degrees

    Returns:
        The solar time. Or solar time offset by 12 hours...

    >>> print(local_to_solar(datetime(2016, 1, 1, 12, 0, 0), 'US/Central', -87.6086066062))
    2016-01-01 12:05:50.214814
    >>> print(local_to_solar(datetime(2016, 2, 13, 12, 0, 0), 'US/Central', -87.6086066062))
    2016-02-13 11:54:57.965428
    >>> print(local_to_solar(datetime(2016, 4, 15, 12, 0, 0), 'US/Central', -87.6086066062))
    2016-04-15 11:09:34.603313
    >>> print(local_to_solar(datetime(2016, 10, 30, 12, 0, 0), 'US/Central', -87.6086066062))
    2016-10-30 11:26:01.125734

    """

    "Calculate the time difference, in minutes, between the LSTM and the location."
    delta1 = longitude_correction(zone, lon)

    "Next get the EoT correction"
    delta2 = equation_of_time(dt)

    "Adjust for DST"
    delta3 = -dst(dt, zone)

    "Convert the local time"
    corr = delta1 + delta2 + delta3
    dt_solar = dt + timedelta(minutes=corr)

    return dt_solar


def solar_to_local(dt, zone, lon):
    """
    Convert from solar time to local standard time.

    Args:
        dt (datetime): the local standard time
        zone (str): the local timezone official name
        lon (float): the longitude of the location, in decimal degrees

    Returns:
        The local standard time.

    >>> print(solar_to_local(datetime(2016, 1, 1, 12, 0, 0), 'US/Central', -87.6086066062))
    2016-01-01 11:54:09.785186
    >>> print(solar_to_local(datetime(2016, 2, 13, 12, 0, 0), 'US/Central', -87.6086066062))
    2016-02-13 12:05:02.034572
    >>> print(solar_to_local(datetime(2016, 4, 15, 12, 0, 0), 'US/Central', -87.6086066062))
    2016-04-15 12:50:25.396687
    >>> print(solar_to_local(datetime(2016, 10, 30, 12, 0, 0), 'US/Central', -87.6086066062))
    2016-10-30 12:33:58.874266

    """

    "Get the LSTM for the local time zone"
    lstm = local_standard_time_meridian(zone)
    if not lstm:
        return None

    "Calculate the time difference, in minutes, between the LSTM and the location."
    delta1 = 4. * (lon - lstm)

    "Next get the EoT correction"
    delta2 = equation_of_time(dt)

    "Adjust for DST"
    delta3 = -dst(dt, zone)

    "Convert the local time"
    corr = delta1 + delta2 + delta3
    dt_local = dt - timedelta(minutes=corr)

    return dt_local


def mean_solar_to_local(dt, zone, lon):
    """
    Convert from mean solar time to local standard time.

    Args:
        dt (datetime): the local standard time
        zone (str): the local timezone official name
        lon (float): the longitude of the location, in decimal degrees

    Returns:
        The local standard time.

    >>> print(solar_to_local(datetime(2016, 1, 1, 12, 0, 0), 'US/Central', -87.6086066062))
    2016-01-01 11:54:09.785186

    """

    "Get the LSTM for the local time zone"
    lstm = local_standard_time_meridian(zone)
    if not lstm:
        return None

    "Calculate the time difference, in minutes, between the LSTM and the location."
    delta1 = 4. * (lon - lstm)

    "Adjust for DST"
    delta2 = -dst(dt, zone)

    "Convert the local time"
    corr = delta1 + delta2
    dt_local = dt - timedelta(minutes=corr)

    return dt_local


def utc_to_mean_solar(dt, zone, lon):
    """
    Convert from UTC to mean solar.

    Args:
        dt (datetime): the utc time
        zone (str): the local timezone official name
        lon (float): the longitude of the location, in decimal degrees
        
    Returns:
        The mean solar time corresponding to this utc.

    >>> print(utc_to_mean_solar(datetime(2016, 7, 4, 0, 0, 0, 0), 'US/Central', -87.6086066062))
    2016-07-03 18:09:33.934415-05:00
    
    """
    
    "Convert from utc to local standard time"
    local = utc_to_local(dt, zone)
    if not local:
        return None
        
    "Convert from local standard time to solar"
    solar = local_to_mean_solar(local, zone, lon)
    
    return solar
    

def mean_solar_to_utc(dt, zone, lon):
    """
    Convert from mean solar time to UTC.

    Args:
        dt (datetime): the mean solar time
        zone (str): the local timezone official name
        lon (float): the longitude of the location, in decimal degrees
        
    Returns:
        The utc time corresponding to this mean solar time.

    >>> print(mean_solar_to_utc(datetime(2016, 7, 4, 12, 0, 0, 0), 'US/Central', -87.6086066062))
    2016-07-04 17:50:26.065585+00:00
    
    """
    
    "Convert from mean solar to local standard time"
    local = mean_solar_to_local(dt, zone, lon)
    if not local:
        return None
        
    "Convert from local standard time to utc"
    utc = local_to_utc(local, zone)
    
    return utc
    

def hour_angle(solar_dt):
    """
    Calculate the hour angle of the (mean) sun from the local solar time.
    
    Args:
        solar_dt (datetime): the local solar datetime
        
    Returns:
        The hour angle at this solar time, in decimal degrees.
        
    >>> print(hour_angle(datetime(2016, 5, 18, 8, 30, 0, 0)))
    -52.5
    >>> print(hour_angle(datetime(2016, 5, 18, 12, 0, 0, 0)))
    0.0
    >>> print(hour_angle(datetime(2016, 5, 18, 18, 45, 0, 0)))
    101.25
    
    """
    hour = solar_dt.hour + solar_dt.minute/60. + solar_dt.second/3600. + solar_dt.microsecond/3600000000.
    return 15.*(hour - 12)


if __name__ == '__main__':
    import doctest
    doctest.testmod()
