# Israel Transverse Mercator
Israel Transverse Mercator - Nuget package.

This Package is taken from the following address:
http://tx.technion.ac.il/~zvikabh/software/ITM/

I did very little refactoring of the code that was found here:
http://tx.technion.ac.il/~zvikabh/software/ITM/isr84lib.cs

In order to use this converter please follow these steps:
<code>

    Using IsraelTransverseMercator;

    ...
    ICoordinatesConverter converter = new CoordinatesConverter();

    var northEast = converter.Wgs84ToItm(new LatLon { Latitude = 31.99702701, Longitude = 34.9986170 });
    var northEast = converter.Wgs84ToIcs(new LatLon { Latitude = 31.99702701, Longitude = 34.9986170 });
    var latlon = converter.ItmToWgs84(new NorthEast { North = 656000, East = 200000 });
    var latlon = converter.Wgs84ToIcs(new NorthEast { North = 656000, East = 200000 });

</code>

You can download this package from [nuget](https://www.nuget.org/packages/IsraelTransverseMercator/)

Enjoy!

Harel M.
