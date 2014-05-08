

-module(fdsys_raws).
-author("Martin Vejmelka <vejmelkam@gmail.com>").
-export([retrieve_and_write_obs/3]).



-spec retrieve_and_write_obs(calendar:datetime(),string(),term()) -> ok.
retrieve_and_write_obs(AtGMT,InDir,SSel) ->
  {{Y,M,D},{H,_,_}} = AtGMT,
  Path = lists:flatten(InDir ++ io_lib:format("/raws_ingest_~4..0B~2..0B~2..0B-~2..0B00.csv", [Y,M,D,H])),
  From = fdsys_util:shift_by_seconds({{Y,M,D},{H,0,0}}, -30 * 60),
  To = fdsys_util:shift_by_seconds({{Y,M,D},{H,0,0}}, 30 * 60),
  SSel1 = raws_ingest:resolve_station_selector(SSel),
  Os = raws_ingest:retrieve_observations(SSel1,[fm10],{From,To}),
  error_logger:info_report(fdsys_util:fmt_text("rtma_server: have ~p observations for ~p, writing to ~p", [length(Os),AtGMT,Path])),
  raws_export:obs_to_csv(Os,Path),
  ok.
