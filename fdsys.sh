#!/bin/bash
erl -detached -sasl -kernel error_logger '{file,"log/syserr.log"}' -sname fdsys -pa ebin deps/*/ebin -s fdsys $@

