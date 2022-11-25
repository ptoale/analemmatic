#!/usr/bin/env python
"""
sun: module to calculate position of the Sun

"""
import numpy as np
from astropy.time import Time

J2000 = 2451545
M0 = 357.5291
M1 = 0.98560028
C1 = 1.9148
C2 = 0.0200
C3 = 0.0003
C4 = 0
C5 = 0
C6 = 0
J0 = 0.0009
J1 = 0.0053
J2 = -0.0068
J3 = 1
PI = 102.9372
EPS = 23.4393

def jd_from_iso(t_iso):
    """
    Convert an ISO time string into a Julian date.

    :param t_iso: time in ISO format
    :type t_iso: str
    :return: Julian date
    :rtype: float

    >>> jd_from_iso('2022-11-25 13:09:30')
    2459909.048263889
    """
    t = Time(t_iso)
    return t.jd


def iso_from_jd(jd):
    """
    Convert a Julian date into an ISO time string.

    :param jd: Julian date
    :type jd: float
    :return: time in ISO format
    :rtype: str

    >>> jd = jd_from_iso('2022-11-25 13:09:30')
    >>> iso_from_jd(jd)
    '2022-11-25 13:09:30.000'
    """
    t = Time(jd, format='jd')
    return t.iso


def julian_day(jd):
    """
    Calculate the Julian day from a Julian date.

    :param jd: Julian date
    :type jd: float
    :return: Julian day number
    :rtype: int

    >>> julian_day(jd_from_iso('2022-11-25 13:09:30'))
    8365
    """
    n = int(np.ceil(jd - J2000 + 0.0008))
    return n


def mean_solar_time(jd, lon):
    """
    Get the mean solar time.

    :param jd: Julian date
    :type jd: float
    :param lon: longitude of site, in degrees
    :type lon: float
    :return: mean solar time
    :rtype: float

    >>> mean_solar_time(jd_from_iso('2022-11-25 13:09:30'), -87.61)
    8365.243361111112
    """
    n = julian_day(jd)
    return n - lon/360.


def solar_mean_anomaly(jd, lon):
    """
    Get the mean anomaly of the Sun.

    :param jd: Julian date
    :type jd: float
    :param lon: longitude of site, in degrees
    :type lon: float
    :return: mean anomaly, in degrees
    :rtype: float

    >>> solar_mean_anomaly(jd_from_iso('2022-11-25 13:09:30'), -87.61)
    322.3152989792525
    """
    js = mean_solar_time(jd, lon)
    m = M0 + M1*js
    return m % 360


def equation_of_center(jd, lon):
    """
    Equation of the center.

    :param jd: Julian date
    :type jd: float
    :param lon: longitude of site, in degrees
    :type lon: float
    :return: equation of center, in degrees
    :rtype: float

    >>> equation_of_center(jd_from_iso('2022-11-25 13:09:30'), -87.61)
    -1.190174922925307
    """
    m = np.radians(solar_mean_anomaly(jd, lon))
    return C1*np.sin(m) + C2*np.sin(2*m) + C3*np.sin(3*m) + C4*np.sin(4*m) + C5*np.sin(5*m) + C6*np.sin(6*m)


def solar_true_anomaly(jd, lon):
    """
    Get the true anomaly of the Sun.

    :param jd: Julian date
    :type jd: float
    :param lon: longitude of site, in degrees
    :type lon: float
    :return: true anomaly, in degrees
    :rtype: float

    >>> solar_true_anomaly(jd_from_iso('2022-11-25 13:09:30'), -87.61)
    321.1251240563272
    """
    m = solar_mean_anomaly(jd, lon)
    c = equation_of_center(jd, lon)
    return (m + c) % 360


def ecliptic_longitude(jd, lon):
    """
    Ecliptic longitude of Sun.

    :param jd: Julian date
    :type jd: float
    :param lon: longitude of site, in degrees
    :type lon: float
    :return: ecliptic longitude, in degrees
    :rtype: float

    >>> ecliptic_longitude(jd_from_iso('2022-11-25 13:09:30'), -87.61)
    244.06232405632727
    """
    m = solar_mean_anomaly(jd, lon)
    c = equation_of_center(jd, lon)
    return (m + c + PI + 180) % 360


def solar_transit(jd, lon):
    """
    Julian date of solar transit.

    :param jd: Julian date
    :type jd: float
    :param lon: longitude of site, in degrees
    :type lon: float
    :return: time of solar transit
    :rtype: float

    >>> iso_from_jd(solar_transit(jd_from_iso('2022-11-25 13:09:30'), -87.61))
    '2022-11-26 17:38:04.282'
    """
    js = mean_solar_time(jd, lon)
    m = solar_mean_anomaly(jd, lon)
    lam = ecliptic_longitude(jd, lon)
    return J2000 + js + J1*np.sin(np.radians(m)) + J2*np.sin(np.radians(2*lam))


def solar_declination(jd, lon):
    """
    Declination of Sun.

    :param jd: Julian date
    :type jd: float
    :param lon: longitude of site, in degrees
    :type lon: float
    :return: declination, in degrees
    :rtype: float

    >>> solar_declination(jd_from_iso('2022-11-25 13:09:30'), -87.61)
    -20.959584798935246
    """
    lam = ecliptic_longitude(jd, lon)
    sd = np.sin(np.radians(lam))*np.sin(np.radians(EPS))
    return np.degrees(np.arcsin(sd))


def solar_hour_angle(jd, lat, lon, elev):
    """
    Hour angle of Sun.

    :param jd: Julian date
    :type jd: float
    :param lat: latitude of site, in degrees
    :type lat: float
    :param lon: longitude of site, in degrees
    :type lon: float
    :param elev: elevation of site, in degrees
    :type elev: float
    :return: hour angle, in degrees
    :rtype: float

    >>> solar_hour_angle(jd_from_iso('2022-11-25 13:09:30'), 33.30, -87.61, 85.)
    76.94251996554361
    """
    dec = solar_declination(jd, lon)
    ref = -0.83 - np.sqrt(elev)*2.076/60
    cw = ((np.sin(np.radians(ref)) - np.sin(np.radians(dec))*np.sin(np.radians(lat)))
          / (np.cos(np.radians(dec))*np.cos(np.radians(lat))))
    return np.degrees(np.arccos(cw))


def sun_rise(jd, lat, lon, elev):
    """
    Julian date of Sun rise.

    :param jd: Julian date
    :type jd: float
    :param lat: latitude of site, in degrees
    :type lat: float
    :param lon: longitude of site, in degrees
    :type lon: float
    :param elev: elevation of site, in degrees
    :type elev: float
    :return: time of sun rise
    :rtype: float

    >>> iso_from_jd(sun_rise(jd_from_iso('2022-11-25 13:09:30'), 33.30, -87.61, 85.))
    '2022-11-26 12:30:18.077'
    """
    st = solar_transit(jd, lon)
    ha = solar_hour_angle(jd, lat, lon, elev)
    return st - ha/360


def sun_set(jd, lat, lon, elev):
    """
    Julian date of Sun set.

    :param jd: Julian date
    :type jd: float
    :param lat: latitude of site, in degrees
    :type lat: float
    :param lon: longitude of site, in degrees
    :type lon: float
    :param elev: elevation of site, in degrees
    :type elev: float
    :return: time of sun set
    :rtype: float

    >>> iso_from_jd(sun_set(jd_from_iso('2022-11-25 13:09:30'), 33.30, -87.61, 85.))
    '2022-11-26 22:45:50.487'
    """
    st = solar_transit(jd, lon)
    ha = solar_hour_angle(jd, lat, lon, elev)
    return st + ha/360


if __name__ == '__main__':
    import doctest
    doctest.testmod()
