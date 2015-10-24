using System;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace IsraelTransverseMercator.Tests
{
    [TestClass]
    public class ConverterTests
    {
        [TestMethod]
        public void ConvertWgs84ToItm()
        {
            var converter = new Converter();

            var northEast = converter.Wgs84ToItm(new LatLon { Latitude = 31.99702701, Longitude = 34.9986170 });

            Assert.AreEqual(200000, northEast.East);
            Assert.AreEqual(656000, northEast.North);
        }

        [TestMethod]
        public void ConvertItmToWgs84()
        {
            var converter = new Converter();

            var latlon = converter.ItmToWgs84(new NorthEast { North = 656000, East = 200000 });
            Assert.AreEqual(31.99702701, latlon.Latitude, 1e-7);
            Assert.AreEqual(34.9986170, latlon.Longitude, 1e-7);
        }

    }
}
