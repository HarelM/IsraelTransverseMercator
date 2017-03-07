# Israel Transverse Mercator
Israel Transverse Mercator - Nuget package.

This Package is taken from the following address:
http://tx.technion.ac.il/~zvikabh/software/ITM/

I refactored the code that was found here:
http://tx.technion.ac.il/~zvikabh/software/ITM/isr84lib.cs

In order to use this transfrom please follow these steps:
<code>

    Using IsraelTransverseMercator;
	using GeoAPI.Geometries;

    ...
    IMathTransform mathTransform = new ItmWgs84MathTransfrom();
    var latlon = mathTransform.Transform(new Coordinate(200000, 656000));
	var northEast = mathTransform.Inverse().Transform(new Coordinate(34.9986170, 31.99702701));

</code>

You can download this package from [nuget](https://www.nuget.org/packages/IsraelTransverseMercator/)

Enjoy!

Harel M.
