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
    };
}
