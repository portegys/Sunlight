char *ManPage[] = 
{
"sunlight(1)                                                         sunlight(1)",
"										",
"Calculate the sunlight received on the surface of a planet at some		",
"latitude as a fraction of the maximum possible exposure at that latitude.	",
"The arguments to this function are the date, latitude, tilt of the		",
"planet's axis of rotation relative to the orbital plane.  The function is	",
"invoked at regular angular intervals as the planet revolves around the sun.	",
"It is assumed for simplicity that the planet moves at a constant rate in a	",
"perfect circle about the sun, and that it has a 365 day year.			",
"										",
"Options: -l <angle of latitude (north pole = 90 degrees)>			",
"         -t <angle of planetary tilt (of north pole away from sun in degrees)>	",
"         [-f <frequency of sunlight calculation in days (default = 1 day)>]	",
"     or									",
"         -i (interactive mode)							",
"         [-d <diameter of orbit display in pixels>]				",
"         [-X (X-Windows options to follow)]					",
"     or									",
"         -m (print manual page)						",
"										",
"Terms:										",
"										",
"Angle of latitude: angle north or south from the equatorial plane toward	",
"     the axis of rotation.							",
"Angle of tilt: initial angle of axis of rotation from perpendicular to		",
"     the orbital plane (tilt of north pole away from sun).			",
"										",
"User interface:								",
"										",
"In non-interactive mode, the annual sunlight exposure data is printed on	",
"standard output as a series of date/sunlight pairs, starting at the winter	",
"solstice of December 21, printed at the specified frequency.  Since no		",
"assumptions are made about the rotation rate of the planet, the sunlight is	",
"actually the proportion of the latitude exposed to sunlight, not the length	",
"of the day.  However, for earth it is approximately the length of a day.	",
"										",
"In interactive mode, the user can select and dynamically modify the angles	",
"of latitude and tilt.  The planet is shown orbiting the sun from a perspective	",
"north of the orbital plane, graphically depicting both the selected latitude and",
"its sunlight exposure.  The date and amount of sunlight exposure are also printed.",
"The display is done with the X Windows System.					",
"										",
"Calculation procedure:								",
"										",
"Consider three geometric objects: (1) the sphere of the planet, (2) a plane	",
"cutting the planet at the desired latitude, and (3) a plane which cuts the	",
"planet in half along the line which demarcates sunlight from darkness.  This	",
"latter plane will rotate on an axis through the center of the planet as the	",
"planet revolves around the sun, thus (possibly) varying sunlight exposure at a	",
"latitude throughout the year.							",
"										",
"The intersection of the planes and the sphere may yield either no solution,	",
"a single point, the entire circle of latitude, or a pair of points on the circle.",
"In the first three cases, the latitude either lies entirely exposed or unexposed",
"to sunlight.  In the last case, it is partially exposed and the two points can	",
"be used to compute the fraction of the circle exposed to sunlight.		",
"										",
"Begun on Earth Day, April 22, 1994.						",
(char *)0
};

#define XTFUNCPROTO	// for X function prototypes

#include <stdio.h>
#include <stdlib.h>
#include <sunmath.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <X11/StringDefs.h>
#include <X11/Intrinsic.h>
#include <Xol/OpenLook.h>
#include <Xol/ControlAre.h>
#include <Xol/Form.h>
#include <Xol/OblongButt.h>
#include <Xol/Slider.h>
#include <Xol/StaticText.h>
#include <Xol/Stub.h>

// sunlight calculator for day, planet latitude and tilt
class SunlightCalculator
{

public:

	typedef double Radian;	// radian angle
	typedef struct Pt	// point
	{
		double x,y,z;
	} Point;
	typedef	enum		// type of latitude intersection
	{
		NO_POINT,
		ONE_POINT,
		TWO_POINT,
		ALL_POINT
	} Intersect;

	// calculator inputs: day, planet latitude and tilt
	void calculate(int day, int latitude, int tilt);

	// calculated quantities
	int day;		// day
	Radian orbit_angle;	// planet orbit angle for day (radians)
	Radian latitude;	// latitude (radians)
	Radian tilt;		// tilt (radians)
	double light;		// fraction of latitude in sunlight
	Intersect type;		// type of intersection
	Point int1;		// first intersection point
	Point int2;		// second intersection point
	double xmax,xmin;	// latitude circle projections:
	double ymax,ymin;	// x,y,z maximum and minimum
	double zmax,zmin;
	char *date;		// date string

private:

	// 3-D Euclidean distance
	double dist(Point,Point);

	// construct date string for day
	char *day2date(int);
};
class SunlightCalculator SunLightCalc;	// sunlight calculator object

// sunlight display
class SunlightDisplay
{

public:

	// constructor
	SunlightDisplay()
	{
		this->SUNSCALE = .1;
		this->PLANETSCALE = .1;
	}

	// functions
	void initialize(Widget space,SunlightCalculator *);
	void drawSun();
	void drawPlanet();
	void erasePlanet();

private:

	// sunlight calculator
	SunlightCalculator *calculator;

	// dimensions
	Dimension orbitDiameter;	// diameter of orbit display in pixels
	int sunDiameter;		// sun diameter
	int sunRadius;			// sun radius
	int planetDiameter;		// planet diameter
	int planetRadius;		// planet radius
	int sunPlanetDistance;		// distance from sun to planet
	double SUNSCALE;		// scale of sun to orbit diameter
	double PLANETSCALE;		// scale of planet to orbit diameter

	// X display data
	Window window;			// window
	int screenDepth;		// screen depth
	Widget space;			// space widget
	GC spaceGC;			// graphics for space
	GC sunGC;			// graphics for sun
	GC lightSideGC;			// graphics for light side of planet
	GC darkSideGC;			// graphics for dark side of planet
	GC markerGC;			// graphics for latitude and pole
	GC indicatorGC;			// graphics for light indicator
	GC bwGC;			// black and white graphics

	// specialized drawing functions
	void drawPlanetBW();
	void drawPlanetColor();
};
class SunlightDisplay SunLightDisplay;	// sunlight display object

/* parameters */
int Day = 0;			// current day
int Latitude = 0;		// planet latitude (degrees)
int Tilt = 0;			// planet tilt (degrees)
Dimension Diameter = 700;	// diameter of orbit display in pixels
const int DELAYMAX = 100;	// maximum display delay (.1 secs)
int Delay = DELAYMAX;		// display delay
int Freq = 1;			// frequency of sunlight calculation/display (days)

// usage message
char *Usemsg = "Usage: %s\n\
\t-l <angle of latitude (north pole = 90 degrees)>\n\
\t-t <angle of planetary tilt (of north pole away from sun in degrees)>\n\
\t[-f <frequency of sunlight calculation in days (default = 1 day)>]\n\
or\n\t-i (interactive mode - X Windows graphics)\n\
\t[-d <diameter of orbit display in pixels>]\n\
\t[-X (X-Windows options to follow)]\n\
or\n\t-m (print manual page)\n";

int Argc;	// arguments
char **Argv;

void Interact();	// interactive function

int
main(int argc, char *argv[])
{
	register int i,d,f;
	Bool interactive;

	// get command line options
	Argc = argc;
	Argv = argv;
	Latitude = Tilt = -1;
	d = f = -1;
	interactive = False;
	while ((i = getopt(argc,argv,"d:f:il:mt:X")) != EOF)
	{
		switch(i)
		{
			case 'd':	if (d != -1)
					{
						fprintf(stderr,Usemsg,Argv[0]);
						exit(1);
					}
					if ((d = atoi(optarg)) <= 0)
					{
						fprintf(stderr,Usemsg,Argv[0]);
						fprintf(stderr,"diameter of orbit display must be greater than 0\n");
						exit(1);
					}
					Diameter = (Dimension)d;
					break;
			case 'f':	if (f != -1)
					{
						fprintf(stderr,Usemsg,Argv[0]);
						exit(1);
					}
					if ((f = atoi(optarg)) < 1 || f > 365)
					{
						fprintf(stderr,Usemsg,Argv[0]);
						fprintf(stderr,"frequency of sunlight calculation must be 1-365 days\n");
						exit(1);
					}
					Freq = f;
					break;
			case 'i':	if (interactive == True)
					{
						fprintf(stderr,Usemsg,Argv[0]);
						exit(1);
					}
					interactive = True;
					break;
			case 'l':	if (Latitude != -1)
					{
						fprintf(stderr,Usemsg,Argv[0]);
						exit(1);
					}
					if ((Latitude = atoi(optarg)) < 0 || Latitude > 90)
					{
						fprintf(stderr,Usemsg,Argv[0]);
						fprintf(stderr,"angle of latitude must be 0-90 degrees\n");
						exit(1);
					}
					break;
			case 'm':	for (i = 0; ManPage[i] != NULL; i++)
					{
						printf("%s\n",ManPage[i]);
					}
					exit(0);
			case 't':	if (Tilt != -1)
					{
						fprintf(stderr,Usemsg,Argv[0]);
						exit(1);
					}
					if ((Tilt = atoi(optarg)) < 0 || Tilt > 90)
					{
						fprintf(stderr,Usemsg,Argv[0]);
						fprintf(stderr,"angle of tilt must be 0-90 degrees\n");
						exit(1);
					}
					break;
			case 'X':	// X-Windows options follow
					if (interactive == False)
					{
						fprintf(stderr,Usemsg,Argv[0]);
						exit(1);
					}
					goto gotopt;
			default:	fprintf(stderr,Usemsg,Argv[0]);
					exit(1);
		}
	}

	// re-arrange argc and argv to pass to X
gotopt:	for (i = optind, argc = 1; i < Argc; i++, argc++)
	{
		argv[argc] = argv[i];
	}
	Argc = argc;
	if (interactive == False)
	{
		if (Latitude == -1 || Tilt == -1)
		{
			fprintf(stderr,Usemsg,Argv[0]);
			exit(1);
		}
		if (d != -1)
		{
			fprintf(stderr,Usemsg,Argv[0]);
			exit(1);
		}
	} else {
		if (Latitude != -1 || Tilt != -1)
		{
			fprintf(stderr,Usemsg,Argv[0]);
			exit(1);
		}
		if (f != -1)
		{
			fprintf(stderr,Usemsg,Argv[0]);
			exit(1);
		}
	}

	// interactive?
	if (interactive == True)
	{
		Interact();
		exit(0);
	}

	// non-interactive: print date and sunlight exposure for year,
	// starting at winter solstice.
	for (Day = 0; Day < 365; Day += Freq)
	{
		SunLightCalc.calculate(Day, Latitude, Tilt);
		printf("%s  %d%%\n",SunLightCalc.date,(int)(SunLightCalc.light*100.0));
	}

	exit(0);
}

/*
 * Calculate sunlight plane intersecting latitude circle.
 *
 * Assumptions:
 * a. The origin of the coordinate system is the planet center.
 * b. The radius of the planet is 1.
 * c. The y-axis connects the planet center with the sun (positive toward sun).
 * d. The z-axis is perpendicular to the orbital plane (positive toward north).
 *
 * Use these equations to solve the intersection:
 *
 * 1. Planet sphere:
 *
 *    x**2 + y**2 + z**2 = 1
 *
 * 2. Plane intersecting latitude:
 *
 *    (y * sin(tilt)) + (z * cos(tilt)) = dist_lat
 *
 *    dist_lat == distance to latitude plane from equatorial plane
 *    tilt == angle of planet tilt
 *
 *    This is the plane through the point:
 *
 *        (0,dist_lat*sin(tilt),dist_lat*cos(tilt))
 *
 *    and which is perpendicular to the line:
 *
 *        (0,0,0),(0,sin(tilt),cos(tilt))
 *
 * 3. Rotating sunlight plane:
 *
 *    x = (y / tan(orbit_angle)), when orbit_angle != 0,90,180,270 degrees
 *    y = 0, when orbit_angle = 0 or 180 degrees
 *    x = 0, when orbit_angle = 90 or 270 degrees
 */
void
SunlightCalculator::calculate(int day, int latitude, int tilt)
{
	double lat_dist,rad_dist;
	double a,b,c,d;

	// convert and store inputs
	this->day = day;
	this->orbit_angle = (Radian)(((double)day*2.0*M_PI)/365.0);
	this->latitude = (Radian)(((double)latitude*M_PI)/180.0);
	this->tilt = (Radian)(((double)tilt*M_PI)/180.0);
	this->date = day2date(day);

	// distance to latitude plane from equatorial plane
	if ((lat_dist = sin(this->latitude)) > 1.0) { lat_dist = 1.0; }

	// radius of planet at given latitude
	rad_dist = cos(this->latitude);

	// determine the projection maxima/minima
	this->xmax = rad_dist;
	this->xmin = -rad_dist;
	this->ymax = -(lat_dist*sin(this->tilt))+(rad_dist*cos(this->tilt));
	this->ymin = -(lat_dist*sin(this->tilt))-(rad_dist*cos(this->tilt));
	this->zmax = (lat_dist*cos(this->tilt))+(rad_dist*sin(this->tilt));
	this->zmin = (lat_dist*cos(this->tilt))-(rad_dist*sin(this->tilt));

	// case where sunlight plane is the xz plane
	// (tan(this->orbit_angle) is zero, y = 0)
	if (this->orbit_angle == 0.0 || this->orbit_angle == M_PI)
	{
		// case where tilt == 90 degrees
		if (this->tilt == M_PI_2)
		{
			if (lat_dist == 0.0)
			{
				this->type = ALL_POINT;
				this->light = .5;
				return;
			} else {
				this->type = NO_POINT;
			}
			if (this->orbit_angle == 0.0)
			{
				this->light = 0.0;	// facing away from sun
			} else {
				this->light = 1.0;	// facing toward sun
			}
			return;
		}

		// tilt < 90 degrees
		this->int1.y = 0.0;
		if ((this->int1.z = lat_dist/cos(this->tilt)) > 1.0) { this->int1.z = 1.0; }
		a = this->int1.z;
		a = 1.0-(a*a);
		if (a < 0.0)	// no sqrt, no solution
		{
			this->type = NO_POINT;
			if (this->orbit_angle == 0.0)
			{
				this->light = 0.0;
			} else {
				this->light = 1.0;
			}
			return;
		}
		if ((this->int1.x = sqrt(a)) > 0.0)
		{
			this->type = TWO_POINT;
			this->int2.x = -this->int1.x;
			this->int2.y = this->int1.y;
			this->int2.z = this->int1.z;
			if ((d = dist(int1,int2)/(2.0*rad_dist)) >= 1.0)
			{
				this->light = .5;
				return;
			}
			if (this->orbit_angle == 0.0)
			{
				this->light = asin(d)/M_PI;
			} else {
				this->light = (M_PI-asin(d))/M_PI;
			}
			return;
		}
		this->type = ONE_POINT;
		if (lat_dist >= 1.0)
		{
			// north pole
			this->light = .5;
			return;
		}
		if (this->orbit_angle == 0.0)
		{
			this->light = 0.0;
		} else {
			this->light = 1.0;
		}
		return;
	}

	// case where sunlight plane is the yz plane
	// (tan(this->orbit_angle) is undefined, x = 0)
	if (this->orbit_angle == M_PI_2 || this->orbit_angle == (M_PI+M_PI_2))
	{
		// case where tilt == 90 degrees
		if (this->tilt == M_PI_2)
		{
			this->int1.x = 0.0;
			this->int1.y = lat_dist;
			if ((this->int1.z = sqrt(1.0-(this->int1.y*this->int1.y))) > 0.0)
			{
				this->type = TWO_POINT;
				this->int2.x = this->int1.x;
				this->int2.y = this->int1.y;
				this->int2.z = -this->int1.z;
				this->light = .5;
				return;
			}
			this->type = ONE_POINT;
			this->light = .5;
			return;
		}

		// tilt < 90 degrees - use quadratic equation to solve intersection
		d = sin(this->tilt)/cos(this->tilt);
		d *= d;
		a = 1.0+d;
		b = -2.0*lat_dist*sin(this->tilt);
		d = cos(this->tilt);
		d *= d;
		b /= d;
		c = ((lat_dist*lat_dist)/d)-1.0;
		if ((d = (b*b)-(4.0*a*c)) < 0.0) { d = 0.0; }
		this->int1.x = 0.0;
		this->int1.y = (-b+sqrt(d))/(2.0*a);
		this->int1.z = (lat_dist-(this->int1.y*sin(this->tilt)))/cos(this->tilt);
		if (d > 0.0)
		{
			this->type = TWO_POINT;
			this->int2.x = 0.0;
			this->int2.y = (-b-sqrt(d))/(2.0*a);
			this->int2.z = (lat_dist-(this->int2.y*sin(this->tilt)))/cos(this->tilt);
			this->light = .5;
			return;
		}
		this->type = ONE_POINT;
		this->light = .5;
		return;
	}
	// end of cases where sunlight plane coincident with xyz planes

	// case where tilt == 90 degrees
	if (this->tilt == M_PI_2)
	{
		if ((this->int1.x = lat_dist/tan(this->orbit_angle)) > 1.0) { this->int1.x = 1.0; }
		this->int1.y = lat_dist;
		a = 1.0-(this->int1.x*this->int1.x)-(this->int1.y*this->int1.y);
		if (a < 0.0)
		{
			this->type = NO_POINT;
			if (this->orbit_angle > M_PI_2 && this->orbit_angle < (M_PI+M_PI_2))
			{
				this->light = 1.0;
			} else {
				this->light = 0.0;
			}
			return;
		}
		if ((this->int1.z = sqrt(a)) > 0.0)
		{
			this->type = TWO_POINT;
			this->int2.x = this->int1.x;
			this->int2.y = this->int1.y;
			this->int2.z = -this->int1.z;
			if ((d = dist(int1,int2)/(2.0*rad_dist)) >= 1.0)
			{
				this->light = .5;
				return;
			}
			if (this->orbit_angle > M_PI_2 && this->orbit_angle < (M_PI+M_PI_2))
			{
				this->light = (M_PI-asin(d))/M_PI;
			} else {
				this->light = asin(d)/M_PI;
			}
			return;
		}
		this->type = ONE_POINT;
		if (lat_dist >= 1.0)
		{
			this->light = .5;
			return;
		}
		if (this->orbit_angle > M_PI_2 && this->orbit_angle < (M_PI+M_PI_2))
		{
			this->light = 1.0;
		} else {
			this->light = 0.0;
		}
		return;
	}
	
	// "main" case - use quadratic equation to solve intersection
	a = tan(this->orbit_angle);
	a *= a;
	a = 1.0/a;
	a += 1.0;
	d = sin(this->tilt)/cos(this->tilt);
	d *= d;
	a += d;
	b = -2.0*lat_dist*sin(this->tilt);
	d = cos(this->tilt);
	d *= d;
	b /= d;
	c = ((lat_dist*lat_dist)/d)-1.0;
	d = (b*b)-(4.0*a*c);
	if (d < 0.0)	// no solution?
	{
		this->type = NO_POINT;
		if (this->orbit_angle > M_PI_2 && this->orbit_angle < (M_PI+M_PI_2))
		{
			this->light = 1.0;
		} else {
			this->light = 0.0;
		}
		return;
	}
	this->int1.y = (-b+sqrt(d))/(2.0*a);
	this->int1.x = this->int1.y/tan(this->orbit_angle);
	this->int1.z = (lat_dist-(this->int1.y*sin(this->tilt)))/cos(this->tilt);
	if (d > 0.0)
	{
		this->type = TWO_POINT;
		this->int2.y = (-b-sqrt(d))/(2.0*a);
		this->int2.x = this->int2.y/tan(this->orbit_angle);
		this->int2.z = (lat_dist-(this->int2.y*sin(this->tilt)))/cos(this->tilt);
		if ((d = dist(int1,int2)/(2.0*rad_dist)) >= 1.0)
		{
			this->light = .5;
			return;
		}
		if (this->orbit_angle > M_PI_2 && this->orbit_angle < (M_PI+M_PI_2))
		{
			this->light = (M_PI-asin(d))/M_PI;
		} else {
			this->light = asin(d)/M_PI;
		}
		return;
	}
	this->type = ONE_POINT;
	if (lat_dist >= 1.0)
	{
		this->light = .5;
		return;
	}
	if (this->orbit_angle > M_PI_2 && this->orbit_angle < (M_PI+M_PI_2))
	{
		this->light = 1.0;
	} else {
		this->light = 0.0;
	}
	return;
}

// 3-D Euclidean distance
double
SunlightCalculator::dist(Point int1, Point int2)
{
	double d,t;

	t = int1.x-int2.x;
	t *= t;
	d = t;
	t = int1.y-int2.y;
	t *= t;
	d += t;
	t = int1.z-int2.z;
	t *= t;
	d += t;
	return(sqrt(d));
}

// Construct date string for given day
// Day 0 is considered the winter solstice, December 21.
char *
SunlightCalculator::day2date(int day)
{
	long t;
	char *s;
	static char retbuf[20];
#ifdef NUMTIME
	struct tm *ts;
#endif
	t = (((day+354)%365)+1)*86400;
#ifdef NUMTIME
	ts = localtime(&t);
	sprintf(retbuf,"%02d%02d",ts->tm_mon+1,ts->tm_mday);
#else
	s = ctime(&t);
	s[10] = '\0';
	sprintf(retbuf,"%s",&s[4]);
#endif
	return(retbuf);
}

// widgets
Widget shell,space,container,controlBox,quitButton,
	latitudeSlider,latitudeLabel,latitudeText,
	tiltSlider,tiltLabel,tiltText,
	delaySlider,delayLabel,delayText,
	frequencySlider,frequencyLabel,frequencyText,
	dateSlider,dateLabel,dateText;

// functions
void exposeProc(Widget,XEvent *,Region);
void changeLatitudeCB(Widget,XtPointer,XtPointer);
void changeTiltCB(Widget,XtPointer,XtPointer);
void changeDelayCB(Widget,XtPointer,XtPointer);
void changeFreqCB(Widget,XtPointer,XtPointer);
void changeDateCB(Widget,XtPointer,XtPointer);
void updatePlanet();
void quit(Widget,XtPointer,XtPointer);
void tproc(XtPointer,XtIntervalId);

XtAppContext context;
XtIntervalId timer;

// interactive mode
void
Interact()
{
	Cardinal n;
	Arg args[20];
	char *s;

	// create widgets
	s = Argv[0];
	Argv[0] = "sunlight";
	shell = OlInitialize(Argv[0],NULL,NULL,0,(unsigned *)&Argc,Argv);
	Argv[0] = s;
	n = 0;
	container = XtCreateManagedWidget("container",formWidgetClass,shell,args,n);
	n = 0;
	XtSetArg(args[n],XtNborderWidth,0); n++;
	XtSetArg(args[n],XtNhPad,Diameter/16); n++;
	XtSetArg(args[n],XtNhSpace,Diameter/16); n++;
	XtSetArg(args[n],XtNlayoutType,OL_FIXEDCOLS); n++;
	XtSetArg(args[n],XtNmeasure,3); n++;
	controlBox = XtCreateManagedWidget("controlBox",
			controlAreaWidgetClass,container,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNstring,"Latitude:"); n++;
	latitudeLabel = XtCreateManagedWidget("latitude_label",staticTextWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNorientation,OL_HORIZONTAL); n++;
	XtSetArg(args[n],XtNsliderMax,90); n++;
	XtSetArg(args[n],XtNsliderMin,0); n++;
	XtSetArg(args[n],XtNsliderValue,0); n++;
	XtSetArg(args[n],XtNwidth,Diameter/2); n++;
	latitudeSlider = XtCreateManagedWidget("latitude",sliderWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNstring,"0 degrees"); n++;
	latitudeText = XtCreateManagedWidget("latitude_text",staticTextWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNstring,"Tilt:"); n++;
	tiltLabel = XtCreateManagedWidget("tilt_label",staticTextWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNorientation,OL_HORIZONTAL); n++;
	XtSetArg(args[n],XtNsliderMax,90); n++;
	XtSetArg(args[n],XtNsliderMin,0); n++;
	XtSetArg(args[n],XtNsliderValue,0); n++;
	XtSetArg(args[n],XtNwidth,Diameter/2); n++;
	tiltSlider = XtCreateManagedWidget("tilt",sliderWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNstring,"0 degrees"); n++;
	tiltText = XtCreateManagedWidget("tilt_text",staticTextWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNstring,"Delay:"); n++;
	delayLabel = XtCreateManagedWidget("delay_label",staticTextWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNorientation,OL_HORIZONTAL); n++;
	XtSetArg(args[n],XtNsliderMax,DELAYMAX); n++;
	XtSetArg(args[n],XtNsliderMin,0); n++;
	XtSetArg(args[n],XtNsliderValue,DELAYMAX); n++;
	XtSetArg(args[n],XtNwidth,Diameter/2); n++;
	delaySlider = XtCreateManagedWidget("delay",sliderWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNstring,"STOP"); n++;
	delayText = XtCreateManagedWidget("delay_text",staticTextWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNstring,"Frequency:"); n++;
	frequencyLabel = XtCreateManagedWidget("frequency_label",staticTextWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNorientation,OL_HORIZONTAL); n++;
	XtSetArg(args[n],XtNsliderMax,365); n++;
	XtSetArg(args[n],XtNsliderMin,1); n++;
	XtSetArg(args[n],XtNsliderValue,1); n++;
	XtSetArg(args[n],XtNwidth,Diameter/2); n++;
	frequencySlider = XtCreateManagedWidget("frequency",sliderWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNstring,"1 day"); n++;
	frequencyText = XtCreateManagedWidget("frequency_text",staticTextWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNstring,"Date:"); n++;
	dateLabel = XtCreateManagedWidget("date_label",staticTextWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNorientation,OL_HORIZONTAL); n++;
	XtSetArg(args[n],XtNsliderMax,364); n++;
	XtSetArg(args[n],XtNsliderMin,0); n++;
	XtSetArg(args[n],XtNsliderValue,0); n++;
	XtSetArg(args[n],XtNwidth,Diameter/2); n++;
	dateSlider = XtCreateManagedWidget("date",sliderWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNstring,"Dec 21"); n++;
	dateText = XtCreateManagedWidget("date_text",staticTextWidgetClass,
			controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNlabelJustify,OL_CENTER); n++;
	quitButton = XtCreateManagedWidget("Quit",
			oblongButtonGadgetClass,controlBox,args,n);
	n = 0;
	XtSetArg(args[n],XtNexpose,exposeProc); n++;
	XtSetArg(args[n],XtNborderWidth,2); n++;
	XtSetArg(args[n],XtNyRefName,"controlBox"); n++;
	XtSetArg(args[n],XtNyAddHeight,True); n++;
	XtSetArg(args[n],XtNwidth,Diameter); n++;
	XtSetArg(args[n],XtNheight,Diameter); n++;
	space = XtCreateManagedWidget("space",stubWidgetClass,
			container,args,n);

	// add call backs
	XtAddCallback(latitudeSlider,XtNsliderMoved,
		(XtCallbackProc)changeLatitudeCB,NULL);
	XtAddCallback(tiltSlider,XtNsliderMoved,
		(XtCallbackProc)changeTiltCB,NULL);
	XtAddCallback(delaySlider,XtNsliderMoved,
		(XtCallbackProc)changeDelayCB,NULL);
	XtAddCallback(frequencySlider,XtNsliderMoved,
		(XtCallbackProc)changeFreqCB,NULL);
	XtAddCallback(dateSlider,XtNsliderMoved,
		(XtCallbackProc)changeDateCB,NULL);
	XtAddCallback(quitButton,XtNselect,(XtCallbackProc)quit,NULL);

	context = XtWidgetToApplicationContext(shell);

	XtRealizeWidget(shell);

	// initial sunlight calculation
	SunLightCalc.calculate(Day, Latitude, Tilt);

	// initialize sunlight display
	SunLightDisplay.initialize(space,&SunLightCalc);

	XtAppMainLoop(context);
}

// refresh after window expose
void
exposeProc(Widget w,XEvent *xevent,Region region)
{
	SunLightDisplay.drawSun();
	SunLightDisplay.drawPlanet();
}

// call back for latitude slider
void
changeLatitudeCB(Widget w,XtPointer client_data,XtPointer call_data)
{
	char buf[20];
	Arg arg;

	// change planet latitude
	Latitude = *(int *)call_data;
	updatePlanet();

	// update latitude text
	sprintf(buf,"%d degrees",Latitude);
	XtSetArg(arg,XtNstring,buf);
	XtSetValues(latitudeText,&arg,1);
}

// call back for tilt slider
void
changeTiltCB(Widget w,XtPointer client_data,XtPointer call_data)
{
	char buf[20];
	Arg arg;

	// change planet tilt
	Tilt = *(int *)call_data;
	updatePlanet();

	// update tilt text
	sprintf(buf,"%d degrees",Tilt);
	XtSetArg(arg,XtNstring,buf);
	XtSetValues(tiltText,&arg,1);
}

// call back for delay slider
void
changeDelayCB(Widget w,XtPointer client_data,XtPointer call_data)
{
	char buf[20];
	Arg arg;

	// change timer
	if (Delay < DELAYMAX)
	{
		XtRemoveTimeOut(timer);		// max is stop setting
	}
	if ((Delay = *(int *)call_data) == 0)
	{
		// to prevent tight loop
		timer = XtAppAddTimeOut(context,1,(XtTimerCallbackProc)tproc,NULL);

	} else if (Delay < DELAYMAX)
	{
		timer = XtAppAddTimeOut(context,Delay*100,(XtTimerCallbackProc)tproc,NULL);
	}
	if (Delay < DELAYMAX)
	{
		sprintf(buf,"%d.%d secs",Delay/10,Delay%10);
	} else {
		sprintf(buf,"STOP");
	}
	XtSetArg(arg,XtNstring,buf);
	XtSetValues(delayText,&arg,1);
}

// call back for display frequency slider
void
changeFreqCB(Widget w,XtPointer client_data,XtPointer call_data)
{
	char buf[20];
	Arg arg;

	Freq = *(int *)call_data;
	sprintf(buf,"%d days",Freq);
	XtSetArg(arg,XtNstring,buf);
	XtSetValues(frequencyText,&arg,1);
}

// call back for date slider
void
changeDateCB(Widget w,XtPointer client_data,XtPointer call_data)
{
	char buf[20];
	Arg arg;

	// change date
	Day = *(int *)call_data;
	updatePlanet();

	// update date text
	strcpy(buf,SunLightCalc.date);
	XtSetArg(arg,XtNstring,buf);
	XtSetValues(dateText,&arg,1);
}

// update planet
void
updatePlanet()
{
	// disable timer
	if (Delay < DELAYMAX)
	{
		XtRemoveTimeOut(timer);
	}

	// erase current planet
	SunLightDisplay.erasePlanet();

	// re-calculate
	SunLightCalc.calculate(Day,Latitude,Tilt);

	// display updated planet
	SunLightDisplay.drawPlanet();

	// reset timer
	if (Delay == 0)
	{
		// to prevent tight loop
		timer = XtAppAddTimeOut(context,1,(XtTimerCallbackProc)tproc,NULL);

	} else if (Delay < DELAYMAX)
	{
		timer = XtAppAddTimeOut(context,Delay*100,(XtTimerCallbackProc)tproc,NULL);
	}
}

// call back for quit button
void
quit(Widget w,XtPointer client_data,XtPointer call_data)
{
	exit(0);
}

// process timer expiration
void
tproc(XtPointer client_data,XtIntervalId id)
{
	char buf[20];
	Arg arg;

	// erase old and display new planet position
	SunLightDisplay.erasePlanet();
	Day = (Day+Freq)%365;
	SunLightCalc.calculate(Day,Latitude,Tilt);
	SunLightDisplay.drawPlanet();

	// update date slider and text
	XtSetArg(arg,XtNsliderValue,Day);
	XtSetValues(dateSlider,&arg,1);
	strcpy(buf,SunLightCalc.date);
	XtSetArg(arg,XtNstring,buf);
	XtSetValues(dateText,&arg,1);
	
	// reset timer
	if (Delay == 0)
	{
		timer = XtAppAddTimeOut(context,1,(XtTimerCallbackProc)tproc,NULL);
	} else if (Delay < DELAYMAX)
	{
		timer = XtAppAddTimeOut(context,Delay*100,(XtTimerCallbackProc)tproc,NULL);
	}
}

// initialize sunlight display
void
SunlightDisplay::initialize(Widget space,SunlightCalculator *calculator)
{
	Arg arg;
	int screen;
	Pixel black,white,pixel;
	XGCValues gcValues;
	Colormap colorMap;
	XColor rgbColor,color;

	// store input information
	this->space = space;
	this->calculator = calculator;

	// get diameter of space
	XtSetArg(arg,XtNwidth,&(this->orbitDiameter));
	XtGetValues(space,&arg,1);
	
	// compute dimensions
	this->sunDiameter = (int)((double)(this->orbitDiameter)*(this->SUNSCALE));
	if ((this->sunRadius = this->sunDiameter/2) <= 0)
	{
		fprintf(stderr,"%s: sun cannot be displayed: diameter = %d\n",Argv[0],this->sunDiameter);
		exit(1);
	}
	this->planetDiameter = (int)((double)(this->orbitDiameter)*(this->PLANETSCALE));
	if ((this->planetRadius = this->planetDiameter/2) <= 0)
	{
		fprintf(stderr,"%s: planet cannot be displayed: diameter = %d\n",Argv[0],this->planetDiameter);
		exit(1);
	}
	this->sunPlanetDistance = (this->orbitDiameter/2)-(this->planetDiameter);

	// create graphics contexts for drawing
	this->window = XtWindow(this->space);
	screen = DefaultScreen(OlDefaultDisplay);
	black = BlackPixel(OlDefaultDisplay,screen);
	white = WhitePixel(OlDefaultDisplay,screen);
	this->screenDepth = DefaultDepth(OlDefaultDisplay,screen);
	if (this->screenDepth == 1)
	{
		// monochrome
		gcValues.background = white;
		gcValues.foreground = black;
		this->sunGC = XtGetGC(this->space,GCBackground|GCForeground,&gcValues);
		this->bwGC = XtGetGC(this->space,GCBackground|GCForeground,&gcValues);
		this->indicatorGC = XtGetGC(this->space,GCBackground|GCForeground,&gcValues);
		gcValues.function = GXinvert;
		this->markerGC = XtGetGC(this->space,GCBackground|GCForeground|GCFunction,&gcValues);
		gcValues.background = black;
		gcValues.foreground = white;
		this->spaceGC = XtGetGC(this->space,GCBackground|GCForeground,&gcValues);

	} else {

		// color
		colorMap = DefaultColormap(OlDefaultDisplay,screen);
		if (XAllocNamedColor(OlDefaultDisplay,colorMap,"Yellow",&rgbColor,&color) == 0)
		{
			pixel = white;
		} else {
			pixel = color.pixel;
		}
		gcValues.foreground = pixel;
		this->sunGC = XtGetGC(this->space,GCForeground,&gcValues);
		if (XAllocNamedColor(OlDefaultDisplay,colorMap,"GreenYellow",&rgbColor,&color) == 0)
		{
			pixel = white;
		} else {
			pixel = color.pixel;
		}
		gcValues.foreground = pixel;
		this->lightSideGC = XtGetGC(this->space,GCForeground,&gcValues);
		if (XAllocNamedColor(OlDefaultDisplay,colorMap,"Green",&rgbColor,&color) == 0)
		{
			pixel = black;
		} else {
			pixel = color.pixel;
		}
		gcValues.foreground = pixel;
		this->darkSideGC = XtGetGC(this->space,GCForeground,&gcValues);
		gcValues.foreground = white;
		this->indicatorGC = XtGetGC(this->space,GCForeground,&gcValues);
		gcValues.foreground = black;
		this->markerGC = XtGetGC(this->space,GCForeground,&gcValues);
		if (XAllocNamedColor(OlDefaultDisplay,colorMap,"Blue",&rgbColor,&color) == 0)
		{
			pixel = black;
		} else {
			pixel = color.pixel;
		}
		gcValues.foreground = pixel;
		this->spaceGC = XtGetGC(this->space,GCForeground,&gcValues);
	}
}

extern "C" void sincos(double,double *,double *);

// draw sun in space
void
SunlightDisplay::drawSun()
{
	int x1,y1,x2,y2,a,i;
	double s,c;

	// fill space
	XFillRectangle(OlDefaultDisplay,this->window,this->spaceGC,0,0,
		this->orbitDiameter,this->orbitDiameter);

	// draw sun and rays
	if (this->screenDepth == 1)
	{
		XDrawArc(OlDefaultDisplay,this->window,this->sunGC,
			(this->orbitDiameter/2)-(this->sunRadius),(this->orbitDiameter/2)-(this->sunRadius),
			this->sunDiameter,this->sunDiameter,0,360*64);
	} else {
		XFillArc(OlDefaultDisplay,this->window,this->sunGC,
			(this->orbitDiameter/2)-(this->sunRadius),(this->orbitDiameter/2)-(this->sunRadius),
			this->sunDiameter,this->sunDiameter,0,360*64);
	}
	for (a = i = 0; a < 360; a += 10, i++)
	{
		sincos(((double)a*M_PI)/180.0,&s,&c);
		x1 = (int)((double)(this->sunRadius)*c)+(this->orbitDiameter/2);
		y1 = (int)((double)(this->sunRadius)*s)+(this->orbitDiameter/2);
		if ((i%2) == 0)
		{
			x2 = (int)((double)(2*(this->sunRadius))*c)+(this->orbitDiameter/2);
			y2 = (int)((double)(2*(this->sunRadius))*s)+(this->orbitDiameter/2);
		} else {
			x2 = (int)((double)(3*(this->sunRadius))*c)+(this->orbitDiameter/2);
			y2 = (int)((double)(3*(this->sunRadius))*s)+(this->orbitDiameter/2);
		}
		XDrawLine(OlDefaultDisplay,this->window,this->sunGC,x1,y1,x2,y2);
	}
	XFlush(OlDefaultDisplay);
}

// draw planet
// uses values in calculator
void
SunlightDisplay::drawPlanet()
{
	if (this->screenDepth == 1)
	{
		this->drawPlanetBW();
	} else {
		this->drawPlanetColor();
	}
}

// draw planet in black and white
void
SunlightDisplay::drawPlanetBW()
{
	int x,y;		// planet circle corner
	int a;			// dark side angle
	int lx,ly;		// projected latitude ellipse corner
	int lw,lh;		// latitude width and height
	double s,c;		// sin, cos
	int od;			// orbit diameter
	int pd;			// planet diameter
	int pr;			// planet radius
	int s2p;		// distance from sun to planet
	char buf[10];

	// locate planet and draw
	od = this->orbitDiameter;
	pd = this->planetDiameter;
	pr = this->planetRadius;
	s2p = this->sunPlanetDistance;
	sincos(((double)(this->calculator->day)*2.0*M_PI)/365.0,&s,&c);
	x = (int)((double)s2p*c)+(od/2)-pr;
	y = (int)((double)s2p*s)+(od/2)-pr;
	XDrawArc(OlDefaultDisplay,this->window,this->bwGC,x,y,pd,pd,0,360*64);
	a = 270-((this->calculator->day*360)/365);
	if (a < 0) { a += 360; }
	XFillArc(OlDefaultDisplay,this->window,this->bwGC,x,y,pd,pd,a*64,180*64);

	// use latitude projections to draw ellipse
	if (this->calculator->latitude > 0.0 || this->calculator->tilt > 0.0)
	{
		lx = x+pr-(int)(this->calculator->ymax*(double)pr);
		ly = y+pr-(int)(this->calculator->xmax*(double)pr);
		lw = x+pr-(int)(this->calculator->ymin*(double)pr)-lx;
		lh = y+pr-(int)(this->calculator->xmin*(double)pr)-ly;
		XDrawArc(OlDefaultDisplay,this->window,this->markerGC,lx,ly,lw,lh,0,360*64);
	}

	// draw north pole
	XDrawLine(OlDefaultDisplay,this->window,this->markerGC,
		x+pr+(int)((double)pr*sin(this->calculator->tilt)),y+pr,
		x+pr+(int)((double)pr*1.5*sin(this->calculator->tilt)),y+pr);

	// draw amount of light
	sprintf(buf,"%d%%",(int)(this->calculator->light*100.0));
	XDrawString(OlDefaultDisplay,this->window,this->indicatorGC,x,y,buf,strlen(buf));

	XFlush(OlDefaultDisplay);
}

// draw planet in color
void
SunlightDisplay::drawPlanetColor()
{
	int x,y;		// planet circle corner
	int a;			// light/dark side angles
	int lx,ly;		// projected latitude ellipse corner
	int lw,lh;		// latitude width and height
	double s,c;		// sin, cos
	int od;			// orbit diameter
	int pd;			// planet diameter
	int pr;			// planet radius
	int s2p;		// distance from sun to planet
	char buf[10];

	// locate planet and draw
	od = this->orbitDiameter;
	pd = this->planetDiameter;
	pr = this->planetRadius;
	s2p = this->sunPlanetDistance;
	sincos(((double)(this->calculator->day)*2.0*M_PI)/365.0,&s,&c);
	x = (int)((double)s2p*c)+(od/2)-pr;
	y = (int)((double)s2p*s)+(od/2)-pr;
	XDrawArc(OlDefaultDisplay,this->window,this->darkSideGC,x,y,pd,pd,0,360*64);
	a = 270-((this->calculator->day*360)/365);
	if (a < 0) { a += 360; }
	XFillArc(OlDefaultDisplay,this->window,this->darkSideGC,x,y,pd,pd,a*64,180*64);
	a += 180;
	if (a > 360) { a -= 360; }
	XFillArc(OlDefaultDisplay,this->window,this->lightSideGC,x,y,pd,pd,a*64,180*64);

	// use latitude projections to draw ellipse
	if (this->calculator->latitude > 0.0 || this->calculator->tilt > 0.0)
	{
		lx = x+pr-(int)(this->calculator->ymax*(double)pr);
		ly = y+pr-(int)(this->calculator->xmax*(double)pr);
		lw = x+pr-(int)(this->calculator->ymin*(double)pr)-lx;
		lh = y+pr-(int)(this->calculator->xmin*(double)pr)-ly;
		XDrawArc(OlDefaultDisplay,this->window,this->markerGC,lx,ly,lw,lh,0,360*64);
	}

	// draw north pole
	XDrawLine(OlDefaultDisplay,this->window,this->markerGC,
		x+pr+(int)((double)pr*sin(this->calculator->tilt)),y+pr,
		x+pr+(int)((double)pr*1.5*sin(this->calculator->tilt)),y+pr);

	// draw amount of light
	sprintf(buf,"%d%%",(int)(this->calculator->light*100.0));
	XDrawString(OlDefaultDisplay,this->window,this->indicatorGC,x,y,buf,strlen(buf));

	XFlush(OlDefaultDisplay);
}

// erase planet
// uses values in calculator
void
SunlightDisplay::erasePlanet()
{
	int x,y;		// planet circle corner
	int lx,ly;		// projected latitude ellipse corner
	int lw,lh;		// latitude width and height
	int a;			// dark side angle
	double s,c;		// sin, cos
	int od;			// orbit diameter
	int pd;			// planet diameter
	int pr;			// planet radius
	int s2p;		// distance from sun to planet
	char buf[10];

	// locate planet and erase
	od = this->orbitDiameter;
	pd = this->planetDiameter;
	pr = this->planetRadius;
	s2p = this->sunPlanetDistance;
	sincos(((double)(this->calculator->day)*2.0*M_PI)/365.0,&s,&c);
	x = (int)((double)s2p*c)+(od/2)-pr;
	y = (int)((double)s2p*s)+(od/2)-pr;
	if (this->screenDepth == 1)
	{
		XDrawArc(OlDefaultDisplay,this->window,this->spaceGC,x,y,pd,pd,0,360*64);
		a = 270-((this->calculator->day*360)/365);
		if (a < 0) { a += 360; }
		XFillArc(OlDefaultDisplay,this->window,this->spaceGC,x,y,pd,pd,a*64,180*64);
	} else {
		XDrawArc(OlDefaultDisplay,this->window,this->spaceGC,x,y,pd,pd,0,360*64);
		XFillArc(OlDefaultDisplay,this->window,this->spaceGC,x,y,pd,pd,0,360*64);
	}

	// use latitude projections to erase ellipse
	if (this->calculator->latitude > 0.0 || this->calculator->tilt > 0.0)
	{
		lx = x+pr-(int)(this->calculator->ymax*(double)pr);
		ly = y+pr-(int)(this->calculator->xmax*(double)pr);
		lw = x+pr-(int)(this->calculator->ymin*(double)pr)-lx;
		lh = y+pr-(int)(this->calculator->xmin*(double)pr)-ly;
		XDrawArc(OlDefaultDisplay,this->window,this->spaceGC,lx,ly,lw,lh,0,360*64);
	}

	// erase north pole
	XDrawLine(OlDefaultDisplay,this->window,this->spaceGC,
		x+pr+(int)((double)pr*sin(this->calculator->tilt)),y+pr,
		x+pr+(int)((double)pr*1.5*sin(this->calculator->tilt)),y+pr);

	// erase amount of light
	sprintf(buf,"%d%%",(int)(this->calculator->light*100.0));
	XDrawString(OlDefaultDisplay,this->window,this->spaceGC,x,y,buf,strlen(buf));

	XFlush(OlDefaultDisplay);
}
