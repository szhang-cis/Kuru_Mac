#!/bin/bash

abaqus double job=tube input=tube.inp user=umat.f cpus=3

echo "Se ha ejecutado ABAQUS"
