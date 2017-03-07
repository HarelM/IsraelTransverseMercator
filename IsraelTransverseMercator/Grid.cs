using System;
using System.Linq;

namespace IsraelTransverseMercator
{
    internal class Grid
    {
        /// <summary>
        /// Central meridian in radians
        /// </summary>
        public double CentralLongitude { get; set; }
        /// <summary>
        /// central latitude in radians
        /// </summary>
        public double CentralLatitude { get; set; }
        /// <summary>
        /// Scale factor
        /// </summary>
        public double ScaleFactor { get; set; }
        /// <summary>
        /// False easting
        /// </summary>
        public double FalseEasting { get; set; }
        /// <summary>
        /// False northing
        /// </summary>
        public double FalseNorthing { get; set; }

        /// <summary>
        /// This method is a helper to translate the angle string into double
        /// </summary>
        /// <param name="angle">The angle string e.g. 35.12'16.261"</param>
        /// <returns></returns>
        public static double DegreesStringToRadians(string angle)
        {
            var splitted = angle.Split(new[] { '°', '\'', '"' }).Where(str => string.IsNullOrWhiteSpace(str) == false).ToArray();
            if (splitted.Length != 3)
            {
                throw new ArgumentException("Angle should look like: 35°12'16.261\"");
            }
            double degrees = double.Parse(splitted[0]);
            double minutes = double.Parse(splitted[1]);
            double seconds = double.Parse(splitted[2]);
            return (degrees + (minutes / 60) + (seconds / 3600.0)) / 180.0 * Math.PI;
        }

        public static Grid ITM = new Grid
        {
            CentralLongitude = DegreesStringToRadians("35°12'16.261\""),
            CentralLatitude = DegreesStringToRadians("31°44'03.817\""),
            ScaleFactor = 1.0000067,
            FalseEasting = 219529.584,
            // MAPI says the false northing is 626907.390, and in another place
            // that the meridional arc at the central latitude is 3512424.3388
            FalseNorthing = 2885516.9488
        };
    };
}
