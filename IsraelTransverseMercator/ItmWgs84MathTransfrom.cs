//*****************************************************************************************************
//*                                                                                                   *
//*  This code is free software; you can redistribute it and/or modify it at your will.               *
//*  It is our hope however that if you improve it in any way you will find a way to share it too.     *
//*                                                                                                   *
//*  Original C++ version by jgray77@gmail.com	3/2010								                  *
//*  Ported C# version by mikisiton2@gmail.com	5/2012								                  *
//*  Refactored and created a nuget package    10/2015
//*  Removed ICS, added GeoAPI, .Net Core      03/2017
//*                                                                                                   *
//*  This program is distributed AS-IS in the hope that it will be useful, but WITHOUT ANY WARRANTY;  *
//*  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.        * 
//*                                                                                                   *
//*****************************************************************************************************
//
//
//===================================================================================================
//	Israel Local Grids <==> WGS84 conversions
//===================================================================================================
//
// The Israel New Grid (ITM) is a Transverse Mercator projection of the GRS80 ellipsoid.
//
// To convert from a local grid to WGS84 you first do a "UTM to Lat/Lon" conversion using the 
// known formulas but with the local grid data (Central Meridian, Scale Factor and False 
// Easting and Northing). This results in Lat/Long in the local ellipsoid coordinate system.
// Afterwards you do a Molodensky transformation from this ellipsoid to WGS84.
//
// To convert from WGS84 to a local grid you first do a Molodensky transformation from WGS84
// to the local ellipsoid, after which you do a Lat/Lon to UTM conversion, again with the data of
// the local grid instead of the UTM data.
//
// The UTM to Lat/Lon and Lat/Lon to UTM conversion formulas were taken as-is from the
// excellent article by Prof.Steven Dutch of the University of Wisconsin at Green Bay:
//		http://www.uwgb.edu/dutchs/UsefulData/UTMFormulas.htm
//
// The [abridged] Molodensky transformations were taken from
//		http://home.hiwaay.net/~taylorc/bookshelf/math-science/geodesy/datum/transform/molodensky/
// and can be found in many sources on the net.
// 
// Additional sources:
// ===================
// 1. dX,dY,dZ values:  http://www.geo.hunter.cuny.edu/gis/docs/geographic_transformations.pdf
//
// 2. ITM data:  http://www.mapi.gov.il/geodesy/itm_ftp.txt
//    for the meridional arc false northing, the value is given at
//    http://www.mapi.gov.il/reg_inst/dir2b.doc	
//    (this doc also gives a different formula for Lat/lon -> ITM, but not the reverse)
//
// 3. ICS data:  http://www.mapi.gov.il/geodesy/ics_ftp.txt
//    for the meridional arc false northing, the value is given at several places as the 
//    correction value for Garmin GPS sets, the origin is unknown.
//    e.g. http://www.idobartana.com/etrexkb/etrexisr.htm
//	
// Notes: 
// ======
// 1. The conversions between ICS and ITM are 
//			ITM Lat = ICS Lat - 500000
//			ITM Lon = ICS Lon + 50000
//	  e.g. ITM 678000,230000 <--> ICS 1178000 180000
//
//	  Since the formulas for ITM->WGS84 and ICS->WGS84 are different, the results will differ.
//    For the above coordinates we get the following results (WGS84)
//		ITM->WGS84 32.11'43.945" 35.18'58.782"
//		ICS->WGS84 32.11'43.873" 35.18'58.200"
//      Difference    ~3m            ~15m
//
// 2. If you have, or have seen, formulas that contain the term Sin(1"), I recommend you read 
//    Prof.Dutch's enlightening explanation about it in his link above.
//
//===================================================================================================

using System;
using System.Collections.Generic;
using System.Linq;
using GeoAPI.CoordinateSystems.Transformations;
using GeoAPI.Geometries;

namespace IsraelTransverseMercator
{
    /// <summary>
    /// This class is the main class used to convert between positioning systems.
    /// </summary>
    public class ItmWgs84MathTransfrom : IMathTransform
    {
        private readonly bool _inversed;

        public ItmWgs84MathTransfrom(bool inversed = false)
        {
            _inversed = inversed;
        }

        public bool Identity()
        {
            throw new NotImplementedException();
        }

        public double[,] Derivative(double[] point)
        {
            throw new NotImplementedException();
        }

        public List<double> GetCodomainConvexHull(List<double> points)
        {
            throw new NotImplementedException();
        }

        public DomainFlags GetDomainFlags(List<double> points)
        {
            throw new NotImplementedException();
        }

        public IMathTransform Inverse()
        {
            return new ItmWgs84MathTransfrom(!_inversed);
        }

        public double[] Transform(double[] point)
        {
            var cooridnate = Transform(new Coordinate(point.First(), point.Last()));
            return new[] {cooridnate.X, cooridnate.Y};
        }

        public ICoordinate Transform(ICoordinate coordinate)
        {
            throw new NotImplementedException();
        }

        public Coordinate Transform(Coordinate coordinate)
        {
            if (_inversed)
            {
                // 1. Molodensky WGS84 -> GRS80
                var latLon80 = Molodensky(ToRadians(coordinate), Datum.WGS84, Datum.GRS80);

                // 2. Lat/Lon (GRS80) -> Local Grid (ITM)
                return LatLon2Grid(latLon80, Datum.GRS80, Grid.ITM);
            }
            else
            {
                // 1. Local Grid (ITM) -> GRS80
                var latLon80 = Grid2LatLon(coordinate, Grid.ITM, Datum.GRS80);

                // 2. Molodensky GRS80->WGS84
                var latLon84 = Molodensky(latLon80, Datum.GRS80, Datum.WGS84);

                // final results
                return ToDegrees(latLon84);
            }
        }

        public IList<double[]> TransformList(IList<double[]> points)
        {
            return points.Select(Transform).ToList();
        }

        public IList<Coordinate> TransformList(IList<Coordinate> points)
        {
            return points.Select(Transform).ToList();
        }

        public void Invert()
        {
            throw new NotImplementedException();
        }

        public ICoordinateSequence Transform(ICoordinateSequence coordinateSequence)
        {
            throw new NotImplementedException();
        }

        public int DimSource { get; }
        public int DimTarget { get; }
        public string WKT { get; }
        public string XML { get; }


        private Coordinate Grid2LatLon(Coordinate northEast, Grid from, Datum to)
        {
            double y = northEast.Y + from.FalseNorthing;
            double x = northEast.X - from.FalseEasting;
            double M = y / from.ScaleFactor;

            double a = to.EquatorialEarthRadius;
            double b = to.PolarEarthRadius;
            double e = to.Eccentricity;
            double esq = to.EccentricitySquared;

            double mu = M / (a * (1 - e * e / 4 - 3 * Math.Pow(e, 4) / 64 - 5 * Math.Pow(e, 6) / 256));

            double ee = Math.Sqrt(1 - esq);
            double e1 = (1 - ee) / (1 + ee);
            double j1 = 3 * e1 / 2 - 27 * e1 * e1 * e1 / 32;
            double j2 = 21 * e1 * e1 / 16 - 55 * e1 * e1 * e1 * e1 / 32;
            double j3 = 151 * e1 * e1 * e1 / 96;
            double j4 = 1097 * e1 * e1 * e1 * e1 / 512;

            // Footprint Latitude
            double fp = mu + j1 * Math.Sin(2 * mu) + j2 * Math.Sin(4 * mu) + j3 * Math.Sin(6 * mu) + j4 * Math.Sin(8 * mu);

            double sinfp = Math.Sin(fp);
            double cosfp = Math.Cos(fp);
            double tanfp = sinfp / cosfp;
            double eg = (e * a / b);
            double eg2 = eg * eg;
            double C1 = eg2 * cosfp * cosfp;
            double T1 = tanfp * tanfp;
            double R1 = a * (1 - e * e) / Math.Pow(1 - (e * sinfp) * (e * sinfp), 1.5);
            double N1 = a / Math.Sqrt(1 - (e * sinfp) * (e * sinfp));
            double D = x / (N1 * from.ScaleFactor);

            double Q1 = N1 * tanfp / R1;
            double Q2 = D * D / 2;
            double Q3 = (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * eg2 * eg2) * (D * D * D * D) / 24;
            double Q4 = (61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 3 * C1 * C1 - 252 * eg2 * eg2) * (D * D * D * D * D * D) / 720;
            double Q5 = D;
            double Q6 = (1 + 2 * T1 + C1) * (D * D * D) / 6;
            double Q7 = (5 - 2 * C1 + 28 * T1 - 3 * C1 * C1 + 8 * eg2 * eg2 + 24 * T1 * T1) * (D * D * D * D * D) / 120;
            // result lon
            return new Coordinate
            {
                Y = fp - Q1 * (Q2 - Q3 + Q4),
                X = from.CentralLongitude + (Q5 - Q6 + Q7) / cosfp
            };
        }

        /// <summary>
        /// Lat/Lon to Local Grid conversion
        /// </summary>
        /// <param name="latLon"></param>
        /// <param name="from"></param>
        /// <param name="to"></param>
        /// <returns></returns>
        private Coordinate LatLon2Grid(Coordinate latLon, Datum from, Grid to)
        {
            // Datum data for Lat/Lon to TM conversion
            double a = from.EquatorialEarthRadius;
            double e = from.Eccentricity; 	// sqrt(esq);
            double b = from.PolarEarthRadius;
            double lat = latLon.Y;
            double lon = latLon.X;
            double slat1 = Math.Sin(lat);
            double clat1 = Math.Cos(lat);
            double clat1sq = clat1 * clat1;
            double tanlat1sq = slat1 * slat1 / clat1sq;
            double e2 = e * e;
            double e4 = e2 * e2;
            double e6 = e4 * e2;
            double eg = (e * a / b);
            double eg2 = eg * eg;

            double l1 = 1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256;
            double l2 = 3 * e2 / 8 + 3 * e4 / 32 + 45 * e6 / 1024;
            double l3 = 15 * e4 / 256 + 45 * e6 / 1024;
            double l4 = 35 * e6 / 3072;
            double M = a * (l1 * lat - l2 * Math.Sin(2 * lat) + l3 * Math.Sin(4 * lat) - l4 * Math.Sin(6 * lat));
            //double rho = a*(1-e2) / pow((1-(e*slat1)*(e*slat1)),1.5);
            double nu = a / Math.Sqrt(1 - (e * slat1) * (e * slat1));
            double p = lon - to.CentralLongitude;
            double k0 = to.ScaleFactor;
            // y = northing = K1 + K2p2 + K3p4, where
            double K1 = M * k0;
            double K2 = k0 * nu * slat1 * clat1 / 2;
            double K3 = (k0 * nu * slat1 * clat1 * clat1sq / 24) * (5 - tanlat1sq + 9 * eg2 * clat1sq + 4 * eg2 * eg2 * clat1sq * clat1sq);
            // ING north
            double Y = K1 + K2 * p * p + K3 * p * p * p * p - to.FalseNorthing;

            // x = easting = K4p + K5p3, where
            double K4 = k0 * nu * clat1;
            double K5 = (k0 * nu * clat1 * clat1sq / 6) * (1 - tanlat1sq + eg2 * clat1 * clat1);
            // ING east
            double X = K4 * p + K5 * p * p * p + to.FalseEasting;

            // final rounded results
            return new Coordinate
            {
                Y = Y,
                X = X,
            };
        }

        /// <summary>
        /// Abridged Molodensky transformation between 2 datums
        /// </summary>
        /// <param name="input"></param>
        /// <param name="from"></param>
        /// <param name="to"></param>
        /// <returns></returns>
        private Coordinate Molodensky(Coordinate input, Datum from, Datum to)
        {
            // from->WGS84 - to->WGS84 = from->WGS84 + WGS84->to = from->to
            double dX = from.DeltaX - to.DeltaX;
            double dY = from.DeltaY - to.DeltaY;
            double dZ = from.DeltaZ - to.DeltaZ;

            double slat = Math.Sin(input.Y);
            double clat = Math.Cos(input.Y);
            double slon = Math.Sin(input.X);
            double clon = Math.Cos(input.X);
            double ssqlat = slat * slat;

            //dlat = ((-dx * slat * clon - dy * slat * slon + dz * clat)
            //        + (da * rn * from_esq * slat * clat / from_a)
            //        + (df * (rm * adb + rn / adb )* slat * clat))
            //       / (rm + from.h); 

            double from_f = from.Flatenning;
            double df = to.Flatenning - from_f;
            double from_a = from.EquatorialEarthRadius;
            double da = to.EquatorialEarthRadius - from_a;
            double from_esq = from.EccentricitySquared;
            double adb = 1.0 / (1.0 - from_f);
            double rn = from_a / Math.Sqrt(1 - from_esq * ssqlat);
            double rm = from_a * (1 - from_esq) / Math.Pow((1 - from_esq * ssqlat), 1.5);
            double from_h = 0.0; // we're flat!

            double dlat = (-dX * slat * clon - dY * slat * slon + dZ * clat
                           + da * rn * from_esq * slat * clat / from_a +
                           +df * (rm * adb + rn / adb) * slat * clat) / (rm + from_h);

            // dlon = (-dx * slon + dy * clon) / ((rn + from.h) * clat);
            double dlon = (-dX * slon + dY * clon) / ((rn + from_h) * clat);

            // result lat (radians)
            return new Coordinate
            {
                Y = input.Y + dlat,
                X = input.X + dlon
            };
        }

        /// <summary>
        /// Used to convert latlon from radians to degrees
        /// <param name="latLon">Latitude and Longidute object in radians</param>
        /// <returns>Latitude and Longidute object in degrees</returns>
        /// </summary>
        private Coordinate ToDegrees(Coordinate latLon)
        {
            return new Coordinate
            {
                X = latLon.X * 180 / Math.PI,
                Y = latLon.Y * 180 / Math.PI,
            };
        }

        /// <summary>
        /// Used to convert latlon from degrees to radians
        /// </summary>
        /// <param name="latLon">Latitude and Longidute object in degrees</param>
        /// <returns>Latitude and Longidute object in radians</returns>
        private Coordinate ToRadians(Coordinate latLon)
        {
            return new Coordinate
            {
                X = latLon.X * Math.PI / 180,
                Y = latLon.Y * Math.PI / 180,
            };
        }
    }
}
