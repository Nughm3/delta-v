#!/bin/sh

curl "https://celestrak.org/NORAD/elements/gp.php?GROUP=cosmos-1408-debris&FORMAT=tle" -o data/cosmos-1408.tle
curl "https://celestrak.org/NORAD/elements/gp.php?GROUP=cosmos-2251-debris&FORMAT=tle" -o data/cosmos-2251.tle
curl "https://celestrak.org/NORAD/elements/gp.php?GROUP=fengyun-1c-debris&FORMAT=tle" -o data/fengyun-1c.tle
curl "https://celestrak.org/NORAD/elements/gp.php?GROUP=iridium-33-debris&FORMAT=tle" -o data/iridium-33.tle
