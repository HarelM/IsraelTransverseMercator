namespace IsraelTransverseMercator
{
    /// <summary>
    /// This is the interface of the coordinates converter.
    /// This interface is used mainly for dependency injection and tests.
    /// </summary>
    public interface ICoordinatesConverter
    {
        /// <summary>
        /// Israel New Grid (ITM) to WGS84 conversion
        /// </summary>
        /// <param name="northEast">North East Coordinates in ITM grid</param>
        /// <returns>Latitide and Longitude in WGS84</returns>
        LatLon IcsToWgs84(NorthEast northEast);

        /// <summary>
        /// WGS84 to Israel New Grid (ITM) conversion
        /// </summary>
        /// <param name="latLon">Latitide and Longitude in WGS84</param>
        /// <returns>North East Coordinates in ITM grid</returns>
        LatLon ItmToWgs84(NorthEast northEast);

        /// <summary>
        /// Israel Old Grid (ICS) to WGS84 conversion
        /// </summary>
        /// <param name="northEast">North East Coordinates in ICS grid</param>
        /// <returns>Latitide and Longitude in WGS84</returns>
        NorthEast Wgs84ToIcs(LatLon latLon);

        /// <summary>
        /// WGS84 to Israel Old Grid (ICS) conversion
        /// </summary>
        /// <param name="latLon">Latitide and Longitude in WGS84</param>
        /// <returns>North East Coordinates in ICS grid</returns>
        NorthEast Wgs84ToItm(LatLon latLon);
    }
}