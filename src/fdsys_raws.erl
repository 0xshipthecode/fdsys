

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
  case raws_ingest:retrieve_observations(SSel1,[fm10],{From,To}) of
    Os when is_list(Os) ->
      LogText = fdsys_util:fmt_text("fdsys_raws: have ~p observations for ~p, writing to ~p", [length(Os),AtGMT,Path]),
      error_logger:info_report(LogText),
      raws_export:obs_to_csv(Os,Path),
      ok;
    Err ->
      ErrText = fdsys_util:fmt_text("fdsys_raws: failed to retrieve observations at time ~p with error ~p", [AtGMT, Err]),
      error_logger:error_report(ErrText),
      error
  end.
