#!/bin/bash
erl -sasl -kernel error_logger '{file,"log/syserr.log"}' -sname fdsys -pa ebin deps/*/ebin -s fdsys $@

