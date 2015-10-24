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
    }
}
