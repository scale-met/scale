#!/bin/sh

for nm in 000*; do rm ${nm}/*nc; done
for nm in 000*; do rm ${nm}/*~; done
for nm in 000*; do rm ${nm}/init_LOG_*; done
for nm in 000*; do rm ${nm}/LOG_*; done
