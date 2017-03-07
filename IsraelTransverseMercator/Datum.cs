using System;

namespace IsraelTransverseMercator
{
    internal class Datum
    {
        /// <summary>
        /// Equatorial earth radius (a)
        /// </summary>
        public double EquatorialEarthRadius { get; set; }
        /// <summary>
        /// Polar earth radius (b)
        /// </summary>
        public double PolarEarthRadius { get; set; }
        /// <summary>
        /// Flatenning (f) = (a-b)/a  
        /// </summary>
        public double Flatenning { get { return (EquatorialEarthRadius - PolarEarthRadius) / EquatorialEarthRadius; } }
        /// <summary>
        /// EccentricitySquared (esq) = 1-(b*b)/(a*a)
        /// </summary>
        public double EccentricitySquared { get { return 1 - (PolarEarthRadius * PolarEarthRadius) / (EquatorialEarthRadius * EquatorialEarthRadius); } }
        /// <summary>
        /// Eccentricity (e) = sqrt(esq)
        /// </summary>
        public double Eccentricity { get { return Math.Sqrt(EccentricitySquared); } }
        /// <summary>
        /// Delta X to WGS84
        /// </summary>
        public double DeltaX { get; set; }
        /// <summary>
        /// Delat Y to WGS84
        /// </summary>
        public double DeltaY { get; set; }
        /// <summary>
        /// Delat Z to WGS84
        /// </summary>
        public double DeltaZ { get; set; }

        public static Datum WGS84 = new Datum
        {
            EquatorialEarthRadius = 6378137.0,
            PolarEarthRadius = 6356752.3142,
            DeltaX = 0,
            DeltaY = 0,
            DeltaZ = 0,
        };

        public static Datum GRS80 = new Datum
        {
            EquatorialEarthRadius = 6378137.0,
            PolarEarthRadius = 6356752.3141,
            // deltas to WGS84
            DeltaX = -48,
            DeltaY = 55,
            DeltaZ = 52
        };
    }
}
