using System;
namespace IsraelTransverseMercator
{
    public class LatLon
    {
        public double Latitude { get; set; }
        public double Longitude { get; set; }

        public void ToDegrees()
        {
            Latitude *= 180 / Math.PI;
            Longitude *= 180 / Math.PI;
        }

        public void ToRadians()
        {
            Latitude *= Math.PI / 180;
            Longitude *= Math.PI / 180;
        }
    }
}
