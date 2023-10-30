# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.3.0] - 2023-10-30

### Added

- Added HOURS_TO_DEG constant in the conversions module.
- JulianDate objects now support rounding of the seconds component in the class string methods date() and time().
- The package now has a proper test library, which provides more reliability when changes are made in future versions.

### Fixed

- Restructured Body inheritance model to allow for better flexibility in adding custom methods
  dependent on the body.

### Changed

- Changed sub-package and module structure of the package. This improves organization and internal importing
  optimization.
- Body names changed to capitalized name of the body, for example EARTH_BODY is now Earth, SUN_BODY is now Sun.

## [0.2.0] - 2023-10-18

### Added

- Greatly enhanced pass finding algorithm.
- This change log file.

### Changed

- Restructured import mechanisms.

### Removed

- Removed attempted importing of private symbols when in debug mode.

### Security

- Bumped requests to v2.31.0 due to potentially leaking Proxy-Authorization headers.
- Bumped Pygments to v2.15.0 due to a ReDoS issue.
- Bumped certifi to 2023.07.22 due to removal of e-Tugra root certificate.

## [0.1.2] - 2023-03-19

### Added

- Added getDetails() method to Orbitables.

## [0.1.1] - 2023-03-09

### Added

- Added toTopocentric*() methods to the public API.
- Added velocity vector method to GeoPosition object.

### Changed

- Improved computations for finding horizon times of overhead passes.
- Better SatellitePass to JSON interface.
- Improved efficiency in some topocentric reference frame conversions.
- Bumped required pyevspace version to 0.0.12.5.

### Fixed

- Get next pass methods will now catch a satellite pass if it is currently overhead.
- Fixed a bug in the SGP4 CPython extension that prevented compilation on Linus environments.
- Catches infinite loops that result from passes not being computable for a given 
  satellite/geo-position combination.

## [0.1.0] - 2023-01-21

### Added

- Stable pass finder for circular, prograde, LEO objects.

[Unreleased]: https://github.com/qbizzle68/sattrack/compare/v0.3.0...HEAD
[0.3.0]: https://github.com/qbizzle68/sattrack/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/qbizzle68/sattrack/compare/v0.1.2...v0.2.0
[0.1.2]: https://github.com/qbizzle68/sattrack/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/qbizzle68/sattrack/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/qbizzle68/sattrack/releases/tag/v0.1.0
