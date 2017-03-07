using GeoAPI.Geometries;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace IsraelTransverseMercator.Tests
{
    [TestClass]
    public class ItmWgs84MathTransfromTests
    {
        [TestMethod]
        public void ConvertWgs84ToItm()
        {
            var mathTransform = new ItmWgs84MathTransfrom().Inverse();

            var northEast = mathTransform.Transform(new Coordinate(34.9986170, 31.99702701));

            Assert.AreEqual(200000, northEast.X, 1);
            Assert.AreEqual(656000, northEast.Y, 1);
        }

        [TestMethod]
        public void ConvertItmToWgs84()
        {
            var mathTransform = new ItmWgs84MathTransfrom();

            var latlon = mathTransform.Transform(new Coordinate(200000, 656000));
            Assert.AreEqual(31.99702701, latlon.Y, 1e-7);
            Assert.AreEqual(34.9986170, latlon.X, 1e-7);
        }

        [TestMethod]
        public void ConvertItmToWgs84OfJerusalem()
        {
            var mathTransform = new ItmWgs84MathTransfrom();

            var latlon = mathTransform.Transform(new Coordinate(222286, 631556));
            Assert.AreEqual(31.776747919252124, latlon.Y, 1e-7);
            Assert.AreEqual(35.234383488170444, latlon.X, 1e-7);
        }
    }
}
