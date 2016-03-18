#!/usr/bin/env python3

import subprocess
import sys
import re

#########################
## VCF FILE REORDERING ##
#########################

subprocess.call('grep "#" %s > output2.vcf' % sys.argv[1],shell=True)
subprocess.call('grep -v "#" %s | sort -k2 >> output2.vcf' % sys.argv[1],shell=True)

