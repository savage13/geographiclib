GeographicLib
=============

[![Documentation](https://docs.rs/geographiclib/badge.svg)](https://docs.rs/geographiclib)

A Rust interface to [GeographicLib](https://geographiclib.sourceforge.io/html/)

> A library for solving geodesic problems

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
geographiclib = "0.1.0"
```

## Example

```rust
use geographiclib::Geodesic;
let g = Geodesic::wgs84();
let (lat1, lon1) = (37.87622, -122.23558); // Berkeley, California
let (lat2, lon2) = (-9.4047, 147.1597);    // Port Moresby, New Guinea
let (d_deg, d_m, az1, az2) = g.inverse(lat1, lon1, lat2, lon2);

assert_eq!(d_deg, 96.39996198449684); // Distance in degrees
assert_eq!(d_m, 10700471.955233702);  // Distance in meters
assert_eq!(az1, -96.91639942294974);  // Azimuth at (lat1, lon1)
assert_eq!(az2, -127.32548874543627); // Azimuth at (lat2, lon2)
```

## License

This version is released under the same license as GeographicLib; MIT/X11 License


