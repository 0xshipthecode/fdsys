
-module(fdsys_test).
-author("Martin Vejmelka <vejmelkam@gmail.com>").
-export([execute_history_cycle/2,execute_cycle_range/3]).


-spec execute_cycle_range(calendar:datetime(),pos_integer(),term()) -> [ok|term()].
execute_cycle_range(FromGMT,Count,SSel) ->
  % find the time extent
  Times = lists:map(fun (DH) -> fdsys_util:shift_by_seconds(FromGMT,DH*3600) end, lists:seq(0,Count)),
  % extend both ways by 30 minutes
  Start = fdsys_util:shift_by_seconds(hd(Times),-30*60),
  End = fdsys_util:shift_by_seconds(lists:last(Times),30*60),
  % retrieve all observations of fm-10 in selector during the extended timespan
  raws_ingest:acquire_observations(SSel,[fm10],{Start,End},120),
  % run the data assimilation mechanism
  lists:map(fun (T) -> execute_history_cycle(T,SSel) end, Times).


-spec execute_history_cycle(calendar:datetime(),raws_ingest:station_selector()) -> ok|term().
execute_history_cycle(AtGMT,SSel) ->

  % find current GMT time
  {{Y,M,D},{H,_,_}} = AtGMT,

  % vars are set statically at this time
  Vars = [temp,dewpoint,precip,wind_spd,wind_dir],

  % construct name of input directory for target time
  InDir = lists:flatten(io_lib:format("inputs/~4..0B~2..0B~2..0B-~2..0B00", [Y,M,D,H])),

  % construct name of input directory for origin time [1hr bac]
  AtGMT0={{Y0,M0,D0},{H0,_,_}} = fdsys_util:shift_by_seconds(AtGMT,-3600),
  InDir0 = lists:flatten(io_lib:format("inputs/~4..0B~2..0B~2..0B-~2..0B00", [Y0,M0,D0,H0])),

  % phase 2: asynchronously and concurrently retrieve raws observations & rtma analysis
  Me = self(),
  spawn(fun () -> R = fdsys_raws:retrieve_and_write_obs(AtGMT,InDir,SSel), Me ! {task_done, raws_ingest, R} end),
  DlTasks = lists:flatten(lists:map(fun ({T,Dir}) -> download_if_missing(T,Vars,Dir) end, [{AtGMT,InDir},{AtGMT0,InDir0}])),

  % wait until everything is downloaded and stored
  Result = case wait_for_results([raws_ingest|DlTasks]) of
    ok ->
      % phase 3: run the fmda
      LogPath = lists:flatten(io_lib:format("fmncast-~4..0B~2..0B~2..0B-~2..0B00.log",[Y,M,D,H])),
      error_logger:info_report(fdsys_util:fmt_text("rtma_server: running fmncast.py for GMT time ~p", [AtGMT])),
      os:cmd(io_lib:format("python psrc/fmncast.py ~s ~s state &> outputs/~s", [InDir0,InDir,LogPath])),
      ok;
    Errors ->
      Errors
  end,

  % log result
  error_logger:info_report(fdsys_util:fmt_text("rtma_server: execute_cycle completed at GMT ~p with result ~p", [AtGMT,Result])),
  Result.


-spec download_if_missing(calendar:datetime(),[atom()],string()) -> []|rtma_retr.
download_if_missing(AtGMT,Vars,D) ->
  Me = self(),
  case filelib:is_dir(D) of
    true ->
      [];
    false ->
      spawn(fun () ->
              R = rtma_retr:retrieve_archived(AtGMT,Vars,D),
              Me ! {task_done, rtma_retr, R}
              end),
      rtma_retr
  end.


-spec wait_for_results([atom()],[term()]) -> ok|[term()].
wait_for_results([],[]) ->
  ok;
wait_for_results([],Errs) ->
  Errs;
wait_for_results(Lst,Res) ->
  receive
    {task_done,raws_ingest,ok} ->
      wait_for_results(lists:delete(raws_ingest,Lst),Res);
    {task_done,rtma_retr,[]} ->
      wait_for_results(lists:delete(rtma_retr,Lst),Res);
    {task_done,rtma_retr,Errs} ->
      error_logger:info_report(fdsys_util:fmt_text("fdsys_test: errors encountered during rtma_retr:retrieve, ~p", [Errs])),
      wait_for_results(lists:delete(rtma_retr,Lst),[{rtma_errors,Errs}|Res])
  end.


-spec wait_for_results([atom()]) -> ok|[term()].
wait_for_results(Lst) ->
  wait_for_results(Lst,[]).
