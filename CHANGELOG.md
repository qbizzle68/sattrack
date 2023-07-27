# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- This change log file.

### Changed

- Restructured import mechanisms.

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

[Unreleased]: https://github.com/qbizzle68/sattrack/compare/v0.1.2...HEAD
[0.1.2]: https://github.com/qbizzle68/sattrack/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/qbizzle68/sattrack/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/qbizzle68/sattrack/releases/tag/v0.1.0
