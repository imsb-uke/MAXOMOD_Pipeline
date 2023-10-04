#!/usr/bin/env python3
import sys

for line in open(sys.argv[1]):
	ssline = line.strip().split('\t')
	if ssline[0] == sys.argv[2]:
		print(line.strip())
