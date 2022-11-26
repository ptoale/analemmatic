#!/usr/bin/env python
"""
Sundial module.

Need to differentiate between states:
- latitude -> size/dimensions of dial, gnomon offsets
- longitude -> location of hour markers
- timezone -> local standard time, daylight saving time
- time -> shadow length and direction

"""
from datetime import datetime, timedelta
import numpy as np
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import EarthLocation, get_body, HADec, AltAz
import matplotlib.pyplot as plt
from solartime import mean_solar_to_utc, longitude_correction, local_to_utc


class Sundial(object):
    """
    Sundial class.

    :param width: the width of the dial, in meters
    :type width: float
    :param latitude: the latitude of the dial, in degrees
    :type latitude: float
    :param longitude: the longitude of the dial, in degrees
    :type longitude: float
    :param timezone: the timezone of the dial, as in 'US/Central'
    :type timezone: str
    :param elevation: the elevation of the dial, in meters (optional)
    :type elevation: float
    :param epoch: the year to use for any calculations (optional)
    :type epoch: int

    """

    def __init__(self, width, latitude, longitude, timezone, elevation=0, epoch=2023):
        self.width = width
        self.latitude = latitude
        self.longitude = longitude
        self.timezone = timezone
        self.elevation = elevation
        self.epoch = epoch

        self.location = EarthLocation(lat=self.latitude*u.deg, lon=self.longitude*u.deg, height=self.elevation*u.m)
        self.semi_major = self.width/2.
        self.semi_minor = self.semi_major*np.sin(np.radians(self.latitude))
        self.eclipse = None
        self.gnomon_offsets = None
        self.hour_offset = None
        self.height = None
        self.ax = None

    def add(self, day, hour):
        if self.ax is None:
            fig, self.ax = plt.subplots(1, 1, figsize=(8, 8))
            fig.suptitle('Analemmatic Sundial', fontsize=14, fontweight='bold')
            self.ax.set_title('(${:6.1f}^\circ$,${:6.1f}^\circ$,{:6.1f}m) {}'.format(self.latitude, self.longitude,
                                                                                     self.elevation, self.timezone))
            ax_max = 1.1*self.semi_major
            self.ax.axis([-ax_max, ax_max, -ax_max, ax_max])
            self.ax.plot([-ax_max, ax_max], [0, 0], 'k')
            self.ax.plot([0, 0], [-ax_max, ax_max], 'k')

        self.draw_eclipse()
        self.draw_gnomon_offsets()
        self.draw_time_tics()
        self.draw_shadow(day, hour)

    def draw(self, day=None, hour=None):
        fig, self.ax = plt.subplots(1, 1, figsize=(8, 8))
        fig.suptitle('Analemmatic Sundial', fontsize=14, fontweight='bold')
        self.ax.set_title('(${:6.1f}^\circ$,${:6.1f}^\circ$,{:6.1f}m) {}'.format(self.latitude, self.longitude,
                                                                                 self.elevation, self.timezone))
        ax_max = 1.1*self.semi_major
        self.ax.axis([-ax_max, ax_max, -ax_max, ax_max])
        self.ax.plot([-ax_max, ax_max], [0, 0], 'k')
        self.ax.plot([0, 0], [-ax_max, ax_max], 'k')

        self.draw_eclipse()
        self.draw_gnomon_offsets()
        self.draw_time_tics()
        if day:
            self.draw_shadow(day, hour)

    def show(self):
        plt.show()

    def draw_eclipse(self):
        if self.eclipse is None:
            x = np.linspace(-self.semi_major, self.semi_major, 1000)
            y = self.semi_minor*np.sqrt(1 - (x/self.semi_major)**2)
            self.eclipse = (x, y)
        self.ax.plot(self.eclipse[0], self.eclipse[1], 'r')
        self.ax.plot(self.eclipse[0], -self.eclipse[1], 'r')

    def draw_gnomon_offsets(self):
        if self.gnomon_offsets is None:
            self.gnomon_offsets = []
            date1 = datetime.fromisoformat("%s-01-01 12:00:00.000" % str(self.epoch))
            date2 = datetime.fromisoformat("%s-01-01 12:00:00.000" % str(self.epoch + 1))
            n_days = int((date2 - date1).total_seconds() / (24 * 60 * 60))

            "loop over days"
            for d in range(n_days):
                dt = date1 + timedelta(days=d)
                utc = Time(mean_solar_to_utc(dt, self.timezone, self.longitude), format='datetime')

                ha_dec = HADec(location=self.location, obstime=utc,
                               pressure=1013*u.hPa, temperature=18*u.Celsius, relative_humidity=0.5,
                               obswl=550*u.nm)
                sun = get_body('Sun', time=utc, location=self.location).transform_to(ha_dec)
                sun_dec = np.radians(sun.dec)

                goff = self.semi_major * np.cos(np.radians(self.latitude)) * np.tan(sun_dec)
                self.gnomon_offsets.append(goff)

        goff_min = min(self.gnomon_offsets)
        goff_max = max(self.gnomon_offsets)

        w = 0.1*self.semi_major
        self.ax.plot([-w, w], [goff_min, goff_min], 'r')
        self.ax.plot([-w, w], [goff_max, goff_max], 'r')
        self.ax.plot([0, 0], [goff_min, goff_max], 'r')

    def draw_time_tics(self):
        if self.hour_offset is None:
            self.hour_offset = longitude_correction(self.timezone, self.longitude) / 60.

        w = 3 * 0.015 * self.semi_major
        for h in range(24):
            ha = 15*(h - 12 + self.hour_offset)
            x0 = self.semi_major * np.sin(np.radians(ha))
            y0 = self.semi_minor * np.cos(np.radians(ha))

            ang2 = np.arctan2(x0, y0)
            x1 = x0 + w * np.sin(ang2)
            y1 = y0 + w * np.cos(ang2)

            self.ax.plot([x0, x1], [y0, y1], '-r')

    def draw_shadow(self, day, hour):
        if self.height is None:
            date1 = datetime.fromisoformat("%s-01-01 12:00:00.000" % str(self.epoch))
            date2 = datetime.fromisoformat("%s-01-01 12:00:00.000" % str(self.epoch + 1))
            n_days = int((date2 - date1).total_seconds() / (24 * 60 * 60))

            min_heights = []

            "loop over days"
            for d in range(n_days):
                dt = date1 + timedelta(days=d)
                utc = Time(mean_solar_to_utc(dt, self.timezone, self.longitude), format='datetime')

                az_alt = AltAz(location=self.location, obstime=utc,
                               pressure=1013*u.hPa, temperature=18*u.Celsius, relative_humidity=0.5,
                               obswl=550*u.nm)
                sun = get_body('Sun', time=utc, location=self.location).transform_to(az_alt)
                sun_alt = np.radians(sun.alt)

                loh = 1 / np.tan(sun_alt)
                dist = self.semi_minor - self.gnomon_offsets[d]
                min_heights.append(dist / loh)

            self.height = max(min_heights)

        x0 = 0
        y0 = self.gnomon_offsets[day]
        date1 = datetime.fromisoformat("%s-01-01 00:00:00.000" % str(self.epoch))
        dt = date1 + timedelta(days=day) + timedelta(hours=hour)
        utc = Time(local_to_utc(dt, self.timezone), format='datetime')
        print('local = %s' % dt.isoformat())
        print('utc   = %s' % utc.iso)
        az_alt = AltAz(location=self.location, obstime=utc,
                       pressure=1013*u.hPa, temperature=18*u.Celsius, relative_humidity=0.5, obswl=550*u.nm)
        sun = get_body('Sun', time=utc, location=self.location).transform_to(az_alt)
        print(sun.az, sun.alt)
        if sun.alt > 5*u.deg:
            length = self.height/np.tan(np.radians(sun.alt))
            a = sun.az - 180*u.deg
            if a < 0:
                a += 360*u.deg
            x1 = x0 + length*np.sin(np.radians(a))
            y1 = y0 + length*np.cos(np.radians(a))
            print(x0, y0)
            print(x1, y1)
            self.ax.plot([x0, x1], [y0, y1], '-g')


if __name__ == '__main__':

    width = 1.0
    latitude = 33.30167
    longitude = -87.60750
    timezone = 'US/Central'
    elevation = 85
    year = 2023
    sundial = Sundial(width=width, latitude=latitude, longitude=longitude, timezone=timezone,
                      elevation=elevation, epoch=year)

    "These are the days that EoT=0"
    #sundial.add(106, 16.25)
    #sundial.add(165, 16.25)
    #sundial.add(245, 16.25)
    #sundial.add(360, 16.25)

    """
    for m in range(12):
        d = 30*m
        for h in range(14):
            sundial.add(d, 6+h)
    """

    """
    for m in range(24*60):
        sundial.add(360, m/60.)
    """
    sundial.draw()
    sundial.show()
