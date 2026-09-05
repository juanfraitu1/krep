#!/bin/bash
awk 'NR==1 || FNR>1' "$@"
