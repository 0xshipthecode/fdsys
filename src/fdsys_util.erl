
-module(fdsys_util).
-author("Martin Vejmelka <vejmelkam@gmail.com>").
-export([fmt_text/2,shift_by_seconds/2]).


fmt_text(Txt,Args) ->
  lists:flatten(io_lib:format(Txt,Args)).


shift_by_seconds(T,S) ->
  TS = calendar:datetime_to_gregorian_seconds(T),
  calendar:gregorian_seconds_to_datetime(TS + S).
