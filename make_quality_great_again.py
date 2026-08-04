#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compatibilité CLI : délègue au package mqga.

Usage inchangé :
    python3 make_quality_great_again.py --help
"""
import sys

from mqga.cli import main

if __name__ == "__main__":
	sys.exit(main())
