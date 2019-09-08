Sunlight applet:

Calculate the sunlight received on the surface of a planet at some
latitude as a fraction of the maximum possible exposure at that latitude.
The arguments to this function are the date, latitude, and tilt of the
planet's axis of rotation relative to the orbital plane.  The function is
invoked at regular angular intervals as the planet revolves around the sun.
It is assumed for simplicity that the planet moves at a constant rate in a
perfect circle about the sun, and that it has a 365 day year.

Terms:

Angle of latitude: angle north or south from the equatorial plane toward
     the axis of rotation.
Angle of tilt: initial angle of axis of rotation from perpendicular to
     the orbital plane (tilt of north pole away from sun).

User interface:

The user can select and dynamically modify the angles of latitude and tilt.
The planet is shown orbiting the sun from a perspective north of the orbital
plane, graphically depicting both the selected latitude and its sunlight
exposure.  The date and amount of sunlight exposure are also printed.

Calculation procedure:

Consider three geometric objects: (1) the sphere of the planet, (2) a plane
cutting the planet at the desired latitude, and (3) a plane which cuts the
planet in half along the line which demarcates sunlight from darkness.  This
latter plane will rotate on an axis through the center of the planet as the
planet revolves around the sun, thus (possibly) varying sunlight exposure at a
latitude throughout the year.

The intersection of the planes and the sphere may yield either no solution,
a single point, the entire circle of latitude, or a pair of points on the circle.
n the first three cases, the latitude either lies entirely exposed or unexposed
to sunlight.  In the last case, it is partially exposed and the two points can
be used to compute the fraction of the circle exposed to sunlight.

Original C++ program begun on Earth Day, April 22, 1994.
