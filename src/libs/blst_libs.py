# BLAST related methods
# Created by Bryan White, 2016
#
import datetime

def tstamp():
	st = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')
	return st
