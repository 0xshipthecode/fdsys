
-module(raws_fill).
-author("Martin Vejmelka <vejmelkam@gmail.com>").
-export([raws_fill/0]).

raws_fill() ->
  {ok,Ds} = file:list_dir("inputs"),
  lists:map(fun download_if_needed/1, Ds).


download_if_needed(Dir) ->
  SSel = {station_file, "etc/raws_station_list"},
  [Y1,Y2,Y3,Y4,M1,M2,D1,D2,$-,H1,H2,$0,$0] = Dir,
  Y = list_to_integer([Y1,Y2,Y3,Y4]),
  M = list_to_integer([M1,M2]),
  D = list_to_integer([D1,D2]),
  H = list_to_integer([H1,H2]),

  Path = lists:flatten("inputs/" ++ Dir ++ io_lib:format("/raws_ingest_~4..0B~2..0B~2..0B-~2..0B00.csv", [Y,M,D,H])),
  case filelib:is_regular(Path) of
    true ->
      error_logger:info_report(lists:flatten(io_lib:format("raws_fill: already have data for ~p~n", [[Y,M,D,H]]))),
      ok;
    false ->
      error_logger:info_report(lists:flatten(io_lib:format("raws_fill: retrieving data for ~p~n", [[Y,M,D,H]]))),
      From = fdsys_util:shift_by_seconds({{Y,M,D},{H,0,0}}, -30 * 60),
      To = fdsys_util:shift_by_seconds({{Y,M,D},{H,0,0}}, 30 * 60),
      SSel0 = raws_ingest:resolve_station_selector(SSel),
      io:format("stations selected: ~p~n", [SSel0]),
      raws_ingest:acquire_observations(SSel,[fm10],{From,To},120),
      SSel1 = raws_ingest:resolve_station_selector(SSel),
      Os = raws_ingest:retrieve_observations(SSel1,[fm10],{From,To}),
      raws_export:obs_to_csv(Os,Path),
      error_logger:info_report(lists:flatten(io_lib:format("raws_fill: data for ~p stored.~n", [[Y,M,D,H]]))),
      ok
  end.
