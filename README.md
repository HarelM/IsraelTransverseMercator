# Israel Transverse Mercator
Israel Transverse Mercator - Nuget package.

Please note that this package will no longer be maintained as there is a better way to achieve the [relevant results](https://en.wikipedia.org/wiki/Israeli_Transverse_Mercator) using the following code and [Proj4NetGeoAPI](https://www.nuget.org/packages/ProjNET4GeoAPI) NuGet Package:

```
var coordinateTransformFactory = new CoordinateTransformationFactory();
var coordinateSystemFactory = new CoordinateSystemFactory();
var itmParameters = new List<ProjectionParameter>
{
    new ProjectionParameter("latitude_of_origin", 31.734393611111109123611111111111),
    new ProjectionParameter("central_meridian", 35.204516944444442572222222222222),
    new ProjectionParameter("false_northing", 626907.390),
    new ProjectionParameter("false_easting", 219529.584),
    new ProjectionParameter("scale_factor", 1.0000067)
};

var itmDatum = coordinateSystemFactory.CreateHorizontalDatum("Isreal 1993", DatumType.HD_Geocentric,
    Ellipsoid.GRS80, new Wgs84ConversionInfo(-48, 55, 52, 0, 0, 0, 0));

var itmGeo = coordinateSystemFactory.CreateGeographicCoordinateSystem("ITM", AngularUnit.Degrees, itmDatum,
    PrimeMeridian.Greenwich, new AxisInfo("East", AxisOrientationEnum.East), new AxisInfo("North", AxisOrientationEnum.North));

var itmProjection = coordinateSystemFactory.CreateProjection("Transverse_Mercator", "Transverse_Mercator", itmParameters);
var itm = coordinateSystemFactory.CreateProjectedCoordinateSystem("ITM", itmGeo, itmProjection, LinearUnit.Metre,
    new AxisInfo("East", AxisOrientationEnum.East), new AxisInfo("North", AxisOrientationEnum.North));

var wgs84 = ProjectedCoordinateSystem.WGS84_UTM(36, true).GeographicCoordinateSystem;
_inverseTransform = coordinateTransformFactory.CreateFromCoordinateSystems(wgs84, itm);
_transform = coordinateTransformFactory.CreateFromCoordinateSystems(itm, wgs84);
	   
```




Old infromation related to this project:

This Package is taken from the following address:
http://tx.technion.ac.il/~zvikabh/software/ITM/

I refactored the code that was found here:
http://tx.technion.ac.il/~zvikabh/software/ITM/isr84lib.cs

In order to use this transfrom please follow these steps:

```
    Using IsraelTransverseMercator;
	using GeoAPI.Geometries;

    ...
    IMathTransform mathTransform = new ItmWgs84MathTransfrom();
    var latlon = mathTransform.Transform(new Coordinate(200000, 656000));
	var northEast = mathTransform.Inverse().Transform(new Coordinate(34.9986170, 31.99702701));

```

You can download this package from [nuget](https://www.nuget.org/packages/IsraelTransverseMercator/)

Enjoy!

Harel M.
