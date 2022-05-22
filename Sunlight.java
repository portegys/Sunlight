/*
 *
 * Sunlight applet:
 *
 * Calculate the sunlight received on the surface of a planet at some
 * latitude as a fraction of the maximum possible exposure at that latitude.
 * The arguments to this function are the date, latitude, and tilt of the
 * planet's axis of rotation relative to the orbital plane.  The function is
 * invoked at regular angular intervals as the planet revolves around the sun.
 * It is assumed for simplicity that the planet moves at a constant rate in a
 * perfect circle about the sun, and that it has a 365 day year.
 *
 * Terms:
 *
 * Angle of latitude: angle north or south from the equatorial plane toward
 *   the axis of rotation.
 * Angle of tilt: initial angle of axis of rotation from perpendicular to
 *   the orbital plane (tilt of north pole away from sun).
 *
 * User interface:
 *
 * The user can select and dynamically modify the angles of latitude and tilt.
 * The planet is shown orbiting the sun from a perspective north of the orbital
 * plane, graphically depicting both the selected latitude and its sunlight
 * exposure.  The date and amount of sunlight exposure are also printed.
 *
 * Calculation procedure:
 *
 * Consider three geometric objects: (1) the sphere of the planet, (2) a plane
 * cutting the planet at the desired latitude, and (3) a plane which cuts the
 * planet in half along the line which demarcates sunlight from darkness.  This
 * latter plane will rotate on an axis through the center of the planet as the
 * planet revolves around the sun, thus (possibly) varying sunlight exposure at a
 * latitude throughout the year.
 *
 * The intersection of the planes and the sphere may yield either no solution,
 * a single point, the entire circle of latitude, or a pair of points on the circle.
 * n the first three cases, the latitude either lies entirely exposed or unexposed
 * to sunlight.  In the last case, it is partially exposed and the two points can
 * be used to compute the fraction of the circle exposed to sunlight.
 *
 * Original C++ program begun on Earth Day, April 22, 1994.
 *
 */

import java.applet.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.image.*;
import java.text.*;
import java.util.*;
import javax.swing.*;

// Main applet.
public class Sunlight extends Applet
{
   // Applet information.
   public String getAppletInfo()
   {
      return("Sunlight - Version 1.0, May 2001 (Tom Portegys, portegys@lucent.com)");
   }

   // Parameters.
   private final int    DELAYMAX    = 100;      // maximum display delay (.01 secs)
   private final double SUNSCALE    = .1;       // scale of sun to orbit diameter
   private final double PLANETSCALE = .1;       // scale of planet to orbit diameter

   // Components.
   private Display    display;
   private Controls   controls;
   private Calculator calculator;

   // State.
   private int Day       = 0;           // current day
   private int Latitude  = 0;           // planet latitude (degrees)
   private int Tilt      = 0;           // planet tilt (degrees)
   private int Delay     = DELAYMAX;    // display delay
   private int Frequency = 1;           // frequency of sunlight calculation/display (days)

   // Image.
   private Image    image;
   private Graphics imageGraphics;
   private int      orbitDiameter;      // orbit diameter
   private int      sunDiameter;        // sun diameter
   private int      sunRadius;          // sun radius
   private int      planetDiameter;     // planet diameter
   private int      planetRadius;       // planet radius
   private int      sunPlanetDistance;  // distance from sun to planet

   // Initialize.
   public void init()
   {
      Dimension d;

      // Set image drawing dimensions.
      d                 = getSize();
      orbitDiameter     = Math.min(d.width, d.height);
      sunDiameter       = (int)((double)orbitDiameter * SUNSCALE);
      sunRadius         = sunDiameter / 2;
      planetDiameter    = (int)((double)orbitDiameter * PLANETSCALE);
      planetRadius      = planetDiameter / 2;
      sunPlanetDistance = (orbitDiameter / 2) - planetDiameter;

      // Get sunlight calculator.
      calculator = new Calculator();
      calculator.calculate(Day, Latitude, Tilt);

      // Draw initial space.
      image         = createImage(orbitDiameter, orbitDiameter);
      imageGraphics = image.getGraphics();
      drawSpace();

      // Add display and controls.
      display = new Display(new Dimension(orbitDiameter, orbitDiameter));
      add(display);
      controls = new Controls();
      add(controls);
   }


   // Start.
   public void start()
   {
      display.start();
   }


   // Stop.
   public void stop()
   {
      display.stop();
   }


   // Update planet: erase, calculate new position, and draw.
   synchronized void updatePlanet()
   {
      erasePlanet();
      calculator.calculate(Day, Latitude, Tilt);
      drawPlanet();
      display.display();
   }


   // Draw space.
   private void drawSpace()
   {
      int    x1, y1, x2, y2, a, i;
      double s, c;

      // Fill space.
      imageGraphics.setColor(Color.blue.brighter());
      imageGraphics.fillRect(0, 0, orbitDiameter, orbitDiameter);

      // Draw sun and rays.
      imageGraphics.setColor(Color.yellow);
      imageGraphics.fillOval((orbitDiameter / 2) - sunRadius, (orbitDiameter / 2) - sunRadius,
                             sunDiameter, sunDiameter);
      for (a = i = 0; a < 360; a += 10, i++)
      {
         s  = Math.sin(((double)a * Math.PI) / 180.0);
         c  = Math.cos(((double)a * Math.PI) / 180.0);
         x1 = (int)((double)(sunRadius) * c) + (orbitDiameter / 2);
         y1 = (int)((double)(sunRadius) * s) + (orbitDiameter / 2);
         if ((i % 2) == 0)
         {
            x2 = (int)((double)(2 * (sunRadius)) * c) + (orbitDiameter / 2);
            y2 = (int)((double)(2 * (sunRadius)) * s) + (orbitDiameter / 2);
         }
         else
         {
            x2 = (int)((double)(3 * (sunRadius)) * c) + (orbitDiameter / 2);
            y2 = (int)((double)(3 * (sunRadius)) * s) + (orbitDiameter / 2);
         }
         imageGraphics.drawLine(x1, y1, x2, y2);
      }

      // Draw planet.
      drawPlanet();
   }


   // Draw planet.
   // Uses values in calculator.
   private void drawPlanet()
   {
      int    x, y;                      // planet circle corner
      int    a;                         // light/dark side angles
      int    lx, ly;                    // projected latitude ellipse corner
      int    lw, lh;                    // latitude width and height
      double s, c;                      // sin, cos
      int    od;                        // orbit diameter
      int    pd;                        // planet diameter
      int    pr;                        // planet radius
      int    s2p;                       // distance from sun to planet

      // Locate planet and draw.
      od  = orbitDiameter;
      pd  = planetDiameter;
      pr  = planetRadius;
      s2p = sunPlanetDistance;
      s   = Math.sin(((double)(calculator.day) * 2.0 * Math.PI) / 365.0);
      c   = Math.cos(((double)(calculator.day) * 2.0 * Math.PI) / 365.0);
      x   = (int)((double)s2p * c) + (od / 2) - pr;
      y   = (int)((double)s2p * s) + (od / 2) - pr;
      imageGraphics.setColor(Color.green.darker());
      a = 270 - ((calculator.day * 360) / 365);
      if (a < 0) { a += 360; }
      imageGraphics.fillArc(x, y, pd, pd, a, 180);
      a += 180;
      if (a > 360) { a -= 360; }
      imageGraphics.setColor(Color.green.brighter());
      imageGraphics.fillArc(x, y, pd, pd, a, 180);

      // Use latitude projections to draw ellipse.
      imageGraphics.setColor(Color.black);
      if ((calculator.latitude > 0.0) || (calculator.tilt > 0.0))
      {
         lx = x + pr - (int)(calculator.ymax * (double)pr);
         ly = y + pr - (int)(calculator.xmax * (double)pr);
         lw = x + pr - (int)(calculator.ymin * (double)pr) - lx;
         lh = y + pr - (int)(calculator.xmin * (double)pr) - ly;
         imageGraphics.drawArc(lx, ly, lw, lh, 0, 360);
      }

      // Draw north pole.
      imageGraphics.drawLine(x + pr + (int)((double)pr * Math.sin(calculator.tilt)), y + pr,
                             x + pr + (int)((double)pr * 1.5 * Math.sin(calculator.tilt)), y + pr);

      // Draw amount of light.
      imageGraphics.setColor(Color.white);
      imageGraphics.drawString((int)(calculator.light * 100.0) + "%", x, y);
   }


   // Erase planet.
   // Uses values in calculator.
   private void erasePlanet()
   {
      int    x, y;                      // planet circle corner
      double s, c;                      // sin, cos
      int    od;                        // orbit diameter
      int    pd;                        // planet diameter
      int    pr;                        // planet radius
      int    s2p;                       // distance from sun to planet

      // Locate planet and erase.
      od  = orbitDiameter;
      pd  = planetDiameter;
      pr  = planetRadius;
      s2p = sunPlanetDistance;
      s   = Math.sin(((double)(calculator.day) * 2.0 * Math.PI) / 365.0);
      c   = Math.cos(((double)(calculator.day) * 2.0 * Math.PI) / 365.0);
      x   = (int)((double)s2p * c) + (od / 2) - pr;
      y   = (int)((double)s2p * s) + (od / 2) - pr;
      imageGraphics.setColor(Color.blue);
      imageGraphics.fillRect(x, y, pd + (int)((double)pr * 1.5) + 1, pd + 1);
      imageGraphics.drawString((int)(calculator.light * 100.0) + "%", x, y);
   }


   // Display.
   class Display extends Canvas implements Runnable
   {
      private Thread thread;
      Dimension      size;

      // Constructor.
      public Display(Dimension d)
      {
         size = d;
         setBounds(0, 0, size.width, size.height);
      }


      // Start.
      public void start()
      {
         if (thread == null)
         {
            thread = new Thread(this);
            thread.setPriority(Thread.MIN_PRIORITY);
            thread.start();
         }
      }


      // Stop.
      public synchronized void stop()
      {
         thread = null;
      }


      // Display space.
      void display()
      {
         Graphics g;

         if ((g = getGraphics()) != null)
         {
            paint(g);
         }
      }


      // Paint.
      public void paint(Graphics g)
      {
         g.drawImage(image, (size.width / 2) - (orbitDiameter / 2),
                     (size.height / 2) - (orbitDiameter / 2), this);
      }


      // Run.
      public void run()
      {
         Thread me;

         long timer, counter;

         if ((me = Thread.currentThread()) != thread) { return; }

         // Display space.
         display();

         // Action loop.
         counter = 0;
         while (thread == me)
         {
            // Delay.
            timer = Math.min((Delay * 10), 1000);
            try
            {
               Thread.sleep(timer);
            }
            catch (InterruptedException e) { break; }

            // Time to advance date?
            if (Delay < DELAYMAX) { counter += (timer / 10); }
            if (counter < Delay) { continue; }
            counter = 0;

            // Increment day.
            Day = (Day + Frequency) % 365;

            // Update controls for new date.
            controls.update();

            // Update planet.
            updatePlanet();
         }
         thread = null;
      }
   }      // End Display class.

   // Controls.
   class Controls extends Panel implements AdjustmentListener
   {
      private Panel     latitudePanel;
      private Label     latitudeLabel;
      private Label     latitudeState;
      private Scrollbar latitudeScrollbar;
      private Panel     tiltPanel;
      private Label     tiltLabel;
      private Label     tiltState;
      private Scrollbar tiltScrollbar;
      private Panel     delayPanel;
      private Label     delayLabel;
      private Label     delayState;
      private Scrollbar delayScrollbar;
      private Panel     frequencyPanel;
      private Label     frequencyLabel;
      private Label     frequencyState;
      private Scrollbar frequencyScrollbar;
      private Panel     datePanel;
      private Label     dateLabel;
      private Label     dateState;
      private Scrollbar dateScrollbar;
      private final int bubble = 5;

      // Constructor.
      public Controls()
      {
         setLayout(new GridLayout(5, 2));
         latitudePanel = new Panel();
         latitudePanel.setLayout(new BorderLayout());
         latitudeLabel = new Label("Latitude:", Label.LEFT);
         latitudePanel.add("West", latitudeLabel);
         latitudeState = new Label("0 degrees", Label.RIGHT);
         latitudePanel.add("East", latitudeState);
         add(latitudePanel);
         latitudeScrollbar = new Scrollbar(Scrollbar.HORIZONTAL, 0, bubble, 0, 90 + bubble);
         latitudeScrollbar.addAdjustmentListener(this);
         add(latitudeScrollbar);
         tiltPanel = new Panel();
         tiltPanel.setLayout(new BorderLayout());
         tiltLabel = new Label("Tilt:", Label.LEFT);
         tiltPanel.add("West", tiltLabel);
         tiltState = new Label("0 degrees", Label.RIGHT);
         tiltPanel.add("East", tiltState);
         add(tiltPanel);
         tiltScrollbar = new Scrollbar(Scrollbar.HORIZONTAL, 0, bubble, 0, 90 + bubble);
         tiltScrollbar.addAdjustmentListener(this);
         add(tiltScrollbar);
         delayPanel = new Panel();
         delayPanel.setLayout(new BorderLayout());
         delayLabel = new Label("Delay:", Label.LEFT);
         delayPanel.add("West", delayLabel);
         delayState = new Label("   STOP", Label.RIGHT);
         delayPanel.add("East", delayState);
         add(delayPanel);
         delayScrollbar = new Scrollbar(Scrollbar.HORIZONTAL, DELAYMAX, bubble, 0, DELAYMAX + bubble);
         delayScrollbar.addAdjustmentListener(this);
         add(delayScrollbar);
         frequencyPanel = new Panel();
         frequencyPanel.setLayout(new BorderLayout());
         frequencyLabel = new Label("Frequency:", Label.LEFT);
         frequencyPanel.add("West", frequencyLabel);
         frequencyState = new Label("   1 day", Label.RIGHT);
         frequencyPanel.add("East", frequencyState);
         add(frequencyPanel);
         frequencyScrollbar = new Scrollbar(Scrollbar.HORIZONTAL, 1, bubble, 1, 365 + bubble);
         frequencyScrollbar.addAdjustmentListener(this);
         add(frequencyScrollbar);
         datePanel = new Panel();
         datePanel.setLayout(new BorderLayout());
         dateLabel = new Label("Date:", Label.LEFT);
         datePanel.add("West", dateLabel);
         dateState = new Label("Dec 21", Label.RIGHT);
         datePanel.add("East", dateState);
         add(datePanel);
         dateScrollbar = new Scrollbar(Scrollbar.HORIZONTAL, 0, bubble, 0, 364 + bubble);
         dateScrollbar.addAdjustmentListener(this);
         add(dateScrollbar);
      }


      // Scrollbar listener.
      public void adjustmentValueChanged(AdjustmentEvent evt)
      {
         Object o;

         o = evt.getSource();
         if (o == latitudeScrollbar)
         {
            Latitude = latitudeScrollbar.getValue();
            updatePlanet();
            latitudeState.setText(Latitude + " degrees");
         }
         else if (o == tiltScrollbar)
         {
            Tilt = tiltScrollbar.getValue();
            updatePlanet();
            tiltState.setText(Tilt + " degrees");
         }
         else if (o == delayScrollbar)
         {
            Delay = delayScrollbar.getValue();
            if (Delay < DELAYMAX)
            {
               delayState.setText(Delay / 100 + "." + Delay % 100 + " secs");
            }
            else
            {
               delayState.setText("STOP");
            }
         }
         else if (o == frequencyScrollbar)
         {
            Frequency = frequencyScrollbar.getValue();
            frequencyState.setText(Frequency + " days");
         }
         else if (o == dateScrollbar)
         {
            Day = dateScrollbar.getValue();
            updatePlanet();
            dateState.setText(calculator.date);
         }
      }


      // Update date.
      void update()
      {
         dateScrollbar.setValue(Day);
         dateState.setText(calculator.date);
      }
   }      // End Controls class.

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

   // Sunlight calculator for day, planet latitude and tilt.
   class Calculator
   {
      // 3D Point.
      private class Point
      {
         double x, y, z;

         Point() { x = y = z = 0.0; }
         Point(double x, double y, double z)
         {
            this.x = x;
            this.y = y;
            this.z = z;
         }
      };

      // Intersection types.
      final int NO_POINT  = 0;
      final int ONE_POINT = 1;
      final int TWO_POINT = 2;
      final int ALL_POINT = 3;

      // Calculated quantities.
      int    day;                       // day
      double orbit_angle;               // planet orbit angle for day (radians)
      double latitude;                  // latitude (radians)
      double tilt;                      // tilt (radians)
      double light;                     // fraction of latitude in sunlight
      int    intType;                   // type of intersection
      Point  int1;                      // first intersection point
      Point  int2;                      // second intersection point
      double xmax, xmin;                // latitude circle projections:
      double ymax, ymin;                // x,y,z maximum and minimum
      double zmax, zmin;
      String date;                      // date

      // Calculate for day, planet latitude and tilt.
      void calculate(int day, int l, int t)
      {
         double lat_dist, rad_dist;
         double a, b, c, d;
         double pi  = Math.PI;
         double pi2 = Math.PI / 2.0;

         int1 = new Point();
         int2 = new Point();

         // Convert and store inputs.
         this.day         = day;
         this.orbit_angle = ((double)day * 2.0 * pi) / 365.0;
         this.latitude    = ((double)l * pi) / 180.0;
         this.tilt        = ((double)t * pi) / 180.0;
         day2date((long)day);

         // Distance to latitude plane from equatorial plane.
         if ((lat_dist = Math.sin(latitude)) > 1.0) { lat_dist = 1.0; }

         // Radius of planet at given latitude.
         rad_dist = Math.cos(latitude);

         // Determine the projection maxima/minima.
         xmax = rad_dist;
         xmin = -rad_dist;
         ymax = -(lat_dist * Math.sin(tilt)) + (rad_dist * Math.cos(tilt));
         ymin = -(lat_dist * Math.sin(tilt)) - (rad_dist * Math.cos(tilt));
         zmax = (lat_dist * Math.cos(tilt)) + (rad_dist * Math.sin(tilt));
         zmin = (lat_dist * Math.cos(tilt)) - (rad_dist * Math.sin(tilt));

         // Case where sunlight plane is the xz plane:
         // tan(orbit_angle) is zero, y = 0
         if ((orbit_angle == 0.0) || (orbit_angle == pi))
         {
            // Case where tilt == 90 degrees.
            if (tilt == pi2)
            {
               if (lat_dist == 0.0)
               {
                  intType = ALL_POINT;
                  light   = .5;
                  return;
               }
               else
               {
                  intType = NO_POINT;
               }
               if (orbit_angle == 0.0)
               {
                  light = 0.0;                                  // Facing away from sun.
               }
               else
               {
                  light = 1.0;                                  // Facing toward sun.
               }
               return;
            }

            // Tilt < 90 degrees.
            int1.y = 0.0;
            if ((int1.z = lat_dist / Math.cos(tilt)) > 1.0) { int1.z = 1.0; }
            a = int1.z;
            a = 1.0 - (a * a);
            if (a < 0.0)                        // No sqrt, no solution.
            {
               intType = NO_POINT;
               if (orbit_angle == 0.0)
               {
                  light = 0.0;
               }
               else
               {
                  light = 1.0;
               }
               return;
            }
            if ((int1.x = Math.sqrt(a)) > 0.0)
            {
               intType = TWO_POINT;
               int2.x  = -int1.x;
               int2.y  = int1.y;
               int2.z  = int1.z;
               if ((d = dist(int1, int2) / (2.0 * rad_dist)) >= 1.0)
               {
                  light = .5;
                  return;
               }
               if (orbit_angle == 0.0)
               {
                  light = Math.asin(d) / pi;
               }
               else
               {
                  light = (pi - Math.asin(d)) / pi;
               }
               return;
            }
            intType = ONE_POINT;
            if (lat_dist >= 1.0)
            {
               // North pole.
               light = .5;
               return;
            }
            if (orbit_angle == 0.0)
            {
               light = 0.0;
            }
            else
            {
               light = 1.0;
            }
            return;
         }

         // Case where sunlight plane is the yz plane:
         // tan(orbit_angle) is undefined, x = 0
         if ((orbit_angle == pi2) || (orbit_angle == (pi + pi2)))
         {
            // Case where tilt == 90 degrees.
            if (tilt == pi2)
            {
               int1.x = 0.0;
               int1.y = lat_dist;
               if ((int1.z = Math.sqrt(1.0 - (int1.y * int1.y))) > 0.0)
               {
                  intType = TWO_POINT;
                  int2.x  = int1.x;
                  int2.y  = int1.y;
                  int2.z  = -int1.z;
                  light   = .5;
                  return;
               }
               intType = ONE_POINT;
               light   = .5;
               return;
            }

            // Tilt < 90 degrees - use quadratic equation to solve intersection.
            d  = Math.sin(tilt) / Math.cos(tilt);
            d *= d;
            a  = 1.0 + d;
            b  = -2.0 *lat_dist *Math.sin(tilt);

            d  = Math.cos(tilt);
            d *= d;
            b /= d;
            c  = ((lat_dist * lat_dist) / d) - 1.0;
            if ((d = (b * b) - (4.0 * a * c)) < 0.0) { d = 0.0; }
            int1.x = 0.0;
            int1.y = (-b + Math.sqrt(d)) / (2.0 * a);
            int1.z = (lat_dist - (int1.y * Math.sin(tilt))) / Math.cos(tilt);
            if (d > 0.0)
            {
               intType = TWO_POINT;
               int2.x  = 0.0;
               int2.y  = (-b - Math.sqrt(d)) / (2.0 * a);
               int2.z  = (lat_dist - (int2.y * Math.sin(tilt))) / Math.cos(tilt);
               light   = .5;
               return;
            }
            intType = ONE_POINT;
            light   = .5;
            return;
         }
         // End of cases where sunlight plane coincident with xyz planes.

         // Case where tilt == 90 degrees.
         if (tilt == pi2)
         {
            if ((int1.x = lat_dist / Math.tan(orbit_angle)) > 1.0) { int1.x = 1.0; }
            int1.y = lat_dist;
            a      = 1.0 - (int1.x * int1.x) - (int1.y * int1.y);
            if (a < 0.0)
            {
               intType = NO_POINT;
               if ((orbit_angle > pi2) && (orbit_angle < (pi + pi2)))
               {
                  light = 1.0;
               }
               else
               {
                  light = 0.0;
               }
               return;
            }
            if ((int1.z = Math.sqrt(a)) > 0.0)
            {
               intType = TWO_POINT;
               int2.x  = int1.x;
               int2.y  = int1.y;
               int2.z  = -int1.z;
               if ((d = dist(int1, int2) / (2.0 * rad_dist)) >= 1.0)
               {
                  light = .5;
                  return;
               }
               if ((orbit_angle > pi2) && (orbit_angle < (pi + pi2)))
               {
                  light = (pi - Math.asin(d)) / pi;
               }
               else
               {
                  light = Math.asin(d) / pi;
               }
               return;
            }
            intType = ONE_POINT;
            if (lat_dist >= 1.0)
            {
               light = .5;
               return;
            }
            if ((orbit_angle > pi2) && (orbit_angle < (pi + pi2)))
            {
               light = 1.0;
            }
            else
            {
               light = 0.0;
            }
            return;
         }

         // "Main" case - use quadratic equation to solve intersection.
         a  = Math.tan(orbit_angle);
         a *= a;
         a  = 1.0 / a;
         a += 1.0;
         d  = Math.sin(tilt) / Math.cos(tilt);
         d *= d;
         a += d;
         b  = -2.0 *lat_dist *Math.sin(tilt);

         d  = Math.cos(tilt);
         d *= d;
         b /= d;
         c  = ((lat_dist * lat_dist) / d) - 1.0;
         d  = (b * b) - (4.0 * a * c);
         if (d < 0.0)                   // no solution?
         {
            intType = NO_POINT;
            if ((orbit_angle > pi2) && (orbit_angle < (pi + pi2)))
            {
               light = 1.0;
            }
            else
            {
               light = 0.0;
            }
            return;
         }
         int1.y = (-b + Math.sqrt(d)) / (2.0 * a);
         int1.x = int1.y / Math.tan(orbit_angle);
         int1.z = (lat_dist - (int1.y * Math.sin(tilt))) / Math.cos(tilt);
         if (d > 0.0)
         {
            intType = TWO_POINT;
            int2.y  = (-b - Math.sqrt(d)) / (2.0 * a);
            int2.x  = int2.y / Math.tan(orbit_angle);
            int2.z  = (lat_dist - (int2.y * Math.sin(tilt))) / Math.cos(tilt);
            if ((d = dist(int1, int2) / (2.0 * rad_dist)) >= 1.0)
            {
               light = .5;
               return;
            }
            if ((orbit_angle > pi2) && (orbit_angle < (pi + pi2)))
            {
               light = (pi - Math.asin(d)) / pi;
            }
            else
            {
               light = Math.asin(d) / pi;
            }
            return;
         }
         intType = ONE_POINT;
         if (lat_dist >= 1.0)
         {
            light = .5;
            return;
         }
         if ((orbit_angle > pi2) && (orbit_angle < (pi + pi2)))
         {
            light = 1.0;
         }
         else
         {
            light = 0.0;
         }
      }


      // 3-D Euclidean distance
      private double dist(Point int1, Point int2)
      {
         double d, t;

         t  = int1.x - int2.x;
         t *= t;
         d  = t;
         t  = int1.y - int2.y;
         t *= t;
         d += t;
         t  = int1.z - int2.z;
         t *= t;
         d += t;
         return(Math.sqrt(d));
      }


      // Construct date string (MMM dd) for given day.
      // Day 0 is considered the winter solstice, December 21.
      private void day2date(long day)
      {
         Date             d;
         SimpleDateFormat f;

         d    = new Date((((day + 354) % 365) + 1) * 86400 * 1000);
         f    = new SimpleDateFormat("MMM dd");
         date = f.format(d);
      }
   }      // End Calculator class.

   // Main.
   @SuppressWarnings("deprecation")
   public static void main(String[] args)
   {
      // Create applet.
      Sunlight app = new Sunlight();

      // Create frame.
      JFrame frame = new JFrame();

      frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
      frame.setTitle("Sunlight");
      frame.setBounds(0, 0, 500, 649);
      frame.setLayout(new GridLayout(1, 1));
      frame.add(app);
      frame.setVisible(true);

      // Run applet.
      app.init();
      app.start();
      frame.resize(new Dimension(500, 650));
   }
}
