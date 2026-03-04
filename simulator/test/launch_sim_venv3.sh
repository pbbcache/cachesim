#!/bin/bash
pbbcache_venv=${PBBCACHE_VENV:-${PBBCACHE}/venv3}

if [ "$PBBCACHE" == "" ]; then
	echo "Undefined $PBBCACHE variable. Please export envvars with \". shrc \""
	exit 1
fi

VIRTUALENV=${pbbcache_venv} ${pbbcache_venv}/bin/python ${PBBCACHE}/test/sim.py $*

