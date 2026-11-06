#!/bin/sh
# Assemble the publishable page in dist/: the web page plus the chemistry core it imports.
set -e
cd "$(dirname "$0")"
rm -rf dist && mkdir -p dist/core dist/vendor
cp web/index.html web/*.js dist/
cp web/vendor/* dist/vendor/
cp core/*.js dist/core/
echo "dist assembled:" && find dist -type f | sort
