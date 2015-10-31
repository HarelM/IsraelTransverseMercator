//*****************************************************************************************************
//*                                                                                                   *
//*  This code is free software; you can redistribute it and/or modify it at your will.               *
//*  It is our hope however that if you improve it in any way you will find a way to share it too.     *
//*                                                                                                   *
//*  Original C++ version by jgray77@gmail.com	3/2010								                  *
//*  Ported C# version by mikisiton2@gmail.com	5/2012								                  *
//*  Refactored and created a nuget package    10/2015
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
// The Israel Old Grid (ICS) is a Cassini-Soldner projection of the modified Clark 1880 ellipsoid.
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

namespace IsraelTransverseMercator
{
    /// <summary>
    /// This class is the main class used to convert between positioning systems.
    /// </summary>
    public class CoordinatesConverter : ICoordinatesConverter
    {
        private const string WGS84 = "WGS84";
        private const string GRS80 = "GRS80";
        private const string CLARK80M = "CLARK80M";
        private const string ICS = "ICS";
        private const string ITM = "ITM";

        private Dictionary<string, Datum> Datums;
        private Dictionary<string, Grid> Grids;

        public CoordinatesConverter()
        {
            Datums = new Dictionary<string, Datum>();
            Grids = new Dictionary<string, Grid>();

            Datums.Add(WGS84, new Datum
            {
                EquatorialEarthRadius = 6378137.0,
                PolarEarthRadius = 6356752.3142,
                DeltaX = 0,
                DeltaY = 0,
                DeltaZ = 0,
            });

            Datums.Add(GRS80, new Datum
            {
                EquatorialEarthRadius = 6378137.0,
                PolarEarthRadius = 6356752.3141,
                // deltas to WGS84
                DeltaX = -48,
                DeltaY = 55,
                DeltaZ = 52
            });

            // Clark 1880 Modified data
            Datums.Add(CLARK80M, new Datum
            {
                EquatorialEarthRadius = 6378300.789,
                PolarEarthRadius = 6356566.4116309,
                // deltas to WGS84
                DeltaX = -235,
                DeltaY = -85,
                DeltaZ = 264,
            });

            Grids.Add(ICS, new Grid
            {
                CentralLongitude = Grid.DegreesStringToRadians("35°12'43.490\""),
                CentralLatitude = Grid.DegreesStringToRadians("31°44'02.749\""),
                ScaleFactor = 1.0,
                FalseEasting = 170251.555,
                FalseNorthing = 2385259.0
            });

            // ITM data
            Grids.Add(ITM, new Grid
            {
                CentralLongitude = Grid.DegreesStringToRadians("35°12'16.261\""),
                CentralLatitude = Grid.DegreesStringToRadians("31°44'03.817\""),
                ScaleFactor = 1.0000067,
                FalseEasting = 219529.584,
                // MAPI says the false northing is 626907.390, and in another place
                // that the meridional arc at the central latitude is 3512424.3388
                FalseNorthing = 2885516.9488
            });
        }

        /// <summary>
        /// Israel New Grid (ITM) to WGS84 conversion
        /// </summary>
        /// <param name="northEast">North East Coordinates in ITM grid</param>
        /// <returns>Latitide and Longitude in WGS84</returns>
        public LatLon ItmToWgs84(NorthEast northEast)
        {
            // 1. Local Grid (ITM) -> GRS80
            var latLon80 = Grid2LatLon(northEast, ITM, GRS80);

            // 2. Molodensky GRS80->WGS84
            var latLon84 = Molodensky(latLon80, GRS80, WGS84);

            // final results
            return ToDegrees(latLon84);
        }

        /// <summary>
        /// WGS84 to Israel New Grid (ITM) conversion
        /// </summary>
        /// <param name="latLon">Latitide and Longitude in WGS84</param>
        /// <returns>North East Coordinates in ITM grid</returns>
        public NorthEast Wgs84ToItm(LatLon latLon)
        {
            // 1. Molodensky WGS84 -> GRS80
            var latLon80 = Molodensky(ToRadians(latLon), WGS84, GRS80);

            // 2. Lat/Lon (GRS80) -> Local Grid (ITM)
            return LatLon2Grid(latLon80, GRS80, ITM);
        }

        /// <summary>
        /// Israel Old Grid (ICS) to WGS84 conversion
        /// </summary>
        /// <param name="northEast">North East Coordinates in ICS grid</param>
        /// <returns>Latitide and Longitude in WGS84</returns>
        public LatLon IcsToWgs84(NorthEast northEast)
        {
            // 1. Local Grid (ICS) -> Clark_1880_modified
            var latLon80 = Grid2LatLon(northEast, ICS, CLARK80M);

            // 2. Molodensky Clark_1880_modified -> WGS84
            var latLon84 = Molodensky(latLon80, CLARK80M, WGS84);

            // final results
            return ToDegrees(latLon84);
        }

        /// <summary>
        /// WGS84 to Israel Old Grid (ICS) conversion
        /// </summary>
        /// <param name="latLon">Latitide and Longitude in WGS84</param>
        /// <returns>North East Coordinates in ICS grid</returns>
        public NorthEast Wgs84ToIcs(LatLon latLon)
        {
            // 1. Molodensky WGS84 -> Clark_1880_modified
            var latLon80 = Molodensky(ToRadians(latLon), WGS84, CLARK80M);

            // 2. Lat/Lon (Clark_1880_modified) -> Local Grid (ICS)
            return LatLon2Grid(latLon80, CLARK80M, ICS);
        }

        private static double sin2(double x)
        {
            return Math.Sin(x) * Math.Sin(x);
        }
        private static double cos2(double x)
        {
            return Math.Cos(x) * Math.Cos(x);
        }
        private static double tan2(double x)
        {
            return Math.Tan(x) * Math.Tan(x);
        }
        private static double tan4(double x)
        {
            return tan2(x) * tan2(x);
        }

        /// <summary>
        /// Local Grid to Lat/Lon conversion
        /// </summary>
        /// <param name="northEast"></param>
        /// <param name="from"></param>
        /// <param name="to"></param>
        /// <returns></returns>
        private LatLon Grid2LatLon(NorthEast northEast, string from, string to)
        {
            double y = northEast.North + Grids[from].FalseNorthing;
            double x = northEast.East - Grids[from].FalseEasting;
            double M = y / Grids[from].ScaleFactor;

            double a = Datums[to].EquatorialEarthRadius;
            double b = Datums[to].PolarEarthRadius;
            double e = Datums[to].Eccentricity;
            double esq = Datums[to].EccentricitySquared;

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
            double D = x / (N1 * Grids[from].ScaleFactor);

            double Q1 = N1 * tanfp / R1;
            double Q2 = D * D / 2;
            double Q3 = (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * eg2 * eg2) * (D * D * D * D) / 24;
            double Q4 = (61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 3 * C1 * C1 - 252 * eg2 * eg2) * (D * D * D * D * D * D) / 720;
            double Q5 = D;
            double Q6 = (1 + 2 * T1 + C1) * (D * D * D) / 6;
            double Q7 = (5 - 2 * C1 + 28 * T1 - 3 * C1 * C1 + 8 * eg2 * eg2 + 24 * T1 * T1) * (D * D * D * D * D) / 120;
            // result lon
            return new LatLon
            {
                Latitude = fp - Q1 * (Q2 - Q3 + Q4),
                Longitude = Grids[from].CentralLongitude + (Q5 - Q6 + Q7) / cosfp
            };
        }

        /// <summary>
        /// Lat/Lon to Local Grid conversion
        /// </summary>
        /// <param name="latLon"></param>
        /// <param name="from"></param>
        /// <param name="to"></param>
        /// <returns></returns>
        private NorthEast LatLon2Grid(LatLon latLon, string from, string to)
        {
            // Datum data for Lat/Lon to TM conversion
            double a = Datums[from].EquatorialEarthRadius;
            double e = Datums[from].Eccentricity; 	// sqrt(esq);
            double b = Datums[from].PolarEarthRadius;
            double lat = latLon.Latitude;
            double lon = latLon.Longitude;
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
            double p = lon - Grids[to].CentralLongitude;
            double k0 = Grids[to].ScaleFactor;
            // y = northing = K1 + K2p2 + K3p4, where
            double K1 = M * k0;
            double K2 = k0 * nu * slat1 * clat1 / 2;
            double K3 = (k0 * nu * slat1 * clat1 * clat1sq / 24) * (5 - tanlat1sq + 9 * eg2 * clat1sq + 4 * eg2 * eg2 * clat1sq * clat1sq);
            // ING north
            double Y = K1 + K2 * p * p + K3 * p * p * p * p - Grids[to].FalseNorthing;

            // x = easting = K4p + K5p3, where
            double K4 = k0 * nu * clat1;
            double K5 = (k0 * nu * clat1 * clat1sq / 6) * (1 - tanlat1sq + eg2 * clat1 * clat1);
            // ING east
            double X = K4 * p + K5 * p * p * p + Grids[to].FalseEasting;

            // final rounded results
            return new NorthEast
            {
                North = (int)(Y + 0.5),
                East = (int)(X + 0.5),
            };
        }

        /// <summary>
        /// Abridged Molodensky transformation between 2 datums
        /// </summary>
        /// <param name="input"></param>
        /// <param name="from"></param>
        /// <param name="to"></param>
        /// <returns></returns>
        private LatLon Molodensky(LatLon input, string from, string to)
        {
            // from->WGS84 - to->WGS84 = from->WGS84 + WGS84->to = from->to
            double dX = Datums[from].DeltaX - Datums[to].DeltaX;
            double dY = Datums[from].DeltaY - Datums[to].DeltaY;
            double dZ = Datums[from].DeltaZ - Datums[to].DeltaZ;

            double slat = Math.Sin(input.Latitude);
            double clat = Math.Cos(input.Latitude);
            double slon = Math.Sin(input.Longitude);
            double clon = Math.Cos(input.Longitude);
            double ssqlat = slat * slat;

            //dlat = ((-dx * slat * clon - dy * slat * slon + dz * clat)
            //        + (da * rn * from_esq * slat * clat / from_a)
            //        + (df * (rm * adb + rn / adb )* slat * clat))
            //       / (rm + from.h); 

            double from_f = Datums[from].Flatenning;
            double df = Datums[to].Flatenning - from_f;
            double from_a = Datums[from].EquatorialEarthRadius;
            double da = Datums[to].EquatorialEarthRadius - from_a;
            double from_esq = Datums[from].EccentricitySquared;
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
            return new LatLon
            {
                Latitude = input.Latitude + dlat,
                Longitude = input.Longitude + dlon
            };
        }

        /// <summary>
        /// Used to convert latlon from radians to degrees
        /// <param name="latLon">Latitude and Longidute object in radians</param>
        /// <returns>Latitude and Longidute object in degrees</returns>
        /// </summary>
        private LatLon ToDegrees(LatLon latLon)
        {
            return new LatLon
            {
                Latitude = latLon.Latitude * 180 / Math.PI,
                Longitude = latLon.Longitude * 180 / Math.PI,
            };
        }

        /// <summary>
        /// Used to convert latlon from degrees to radians
        /// </summary>
        /// <param name="latLon">Latitude and Longidute object in degrees</param>
        /// <returns>Latitude and Longidute object in radians</returns>
        private LatLon ToRadians(LatLon latLon)
        {
            return new LatLon
            {
                Latitude = latLon.Latitude * Math.PI / 180,
                Longitude = latLon.Longitude * Math.PI / 180,
            };
        }
    }
}
