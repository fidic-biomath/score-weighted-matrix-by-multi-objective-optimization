#!/bin/bash

f2py -c -m corr1 corr1.f90 ; mv corr1.cpython-310-x86_64-linux-gnu.so corr1.so
f2py -c -m corr2 corr2.f90 ; mv corr2.cpython-310-x86_64-linux-gnu.so corr2.so
f2py -c -m corr3 corr3.f90 ; mv corr3.cpython-310-x86_64-linux-gnu.so corr3.so

