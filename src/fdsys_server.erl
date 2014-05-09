-module(fdsys_server).
-behaviour(gen_server).
-define(SERVER, ?MODULE).

%% ------------------------------------------------------------------
%% API Function Exports
%% ------------------------------------------------------------------

-export([start_link/1]).

-export([test_cycle/1,update_raws_and_store/3]).

%% ------------------------------------------------------------------
%% gen_server Function Exports
%% ------------------------------------------------------------------

-export([init/1, handle_call/3, handle_cast/2, handle_info/2,
         terminate/2, code_change/3]).

%% ------------------------------------------------------------------
%% API Function Definitions
%% ------------------------------------------------------------------

start_link(Cfg) ->
    gen_server:start_link({local, ?SERVER}, ?MODULE, [Cfg], []).

%% ------------------------------------------------------------------
%% gen_server Function Definitions
%% ------------------------------------------------------------------

init(Args) ->
    schedule_next_run(?SERVER),
    {ok, Args}.

handle_call(_Request, _From, State) ->
    {reply, ok, State}.

handle_cast(_Msg, State) ->
    {noreply, State}.

handle_info(_Info, State) ->
    {noreply, State}.

terminate(_Reason, _State) ->
    ok.

code_change(_OldVsn, State, _Extra) ->
    {ok, State}.

%% ------------------------------------------------------------------
%% Internal Function Definitions
%% ------------------------------------------------------------------

schedule_next_run(PID) ->
  {{_,_,_},{_,Min,_}} = calendar:universal_time(),
  TimeoutS = case Min < 20 of
    true ->
      % next run will be H now
      (20 - Min) * 60;
    false ->
      % next run will be NOW
      0
  end,
  timer:send_after(TimeoutS * 1000, PID, cycle_nowcast).


execute_cycle(UrlPfix,AtGMT,SSel) ->

  % find current GMT time
  {{Y,M,D},{H,_,_}} = AtGMT,

  % vars are set statically at this time
  Vars = [temp,dewpoint,precip,wind_spd,wind_dir],

  % phase 1: check if RTMA analysis is already available, if not, wait 60s and retry
  case rtma_retr:is_available(UrlPfix,H,Vars) of
    [] ->
      ok;
    WaitingOnFiles ->
      error_logger:report_info(fdsys_util:fmt_text("fdsys_server: RTMA not ready yet, waiting on ~p", [WaitingOnFiles])),
      timer:sleep(60*1000),
      execute_cycle(UrlPfix,H,SSel)
  end,

  % construct name of input directory for target time
  InDir = lists:flatten(io_lib:format("inputs/~4..0B~2..0B~2..0B-~2..0B00", [Y,M,D,H])),

  % construct name of input directory for origin time [1hr bac]
  {{Y0,M0,D0},{H0,_,_}} = fdsys_util:shift_by_seconds(AtGMT,-3600),
  InDir0 = lists:flatten(io_lib:format("inputs/~4..0B~2..0B~2..0B-~2..0B00", [Y0,M0,D0,H0])),

  % phase 2: asynchronously and concurrently retrieve raws observations & rtma analysis
  Me = self(),
  raws_ingest:update_now(120),
  spawn(fun () -> R = fdsys_raws:retrieve_and_write_obs(AtGMT,InDir,SSel), Me ! {task_done, raws_ingest, R} end),
  DlTasks = lists:flatten(lists:map(fun ({RT,Dir}) -> download_if_missing(UrlPfix,RT,Vars,Dir) end, [{H,InDir},{H0,InDir0}])),

  % wait until everything is downloaded and stored
  Result = case wait_for_results([raws_ingest|DlTasks]) of
    ok ->
      % phase 3: run the fmda
      LogFile = lists:flatten(io_lib:format("fmncast-~4..0B~2..0B~2..0B-~2..0B00.log",[Y,M,D,H])),
      error_logger:info_report(fdsys_util:fmt_text("rtma_server: running fmncast.py for GMT time ~p", [AtGMT])),
      os:cmd(io_lib:format("python psrc/fmncast.py ~s ~s state &> outputs/~s", [InDir0,InDir,LogFile])),
      ok;
    Errors ->
      Errors
  end,

  % log result
  error_logger:info_report(fdsys_util:fmt_text("rtma_server: execute_cycle completed at GMT ~p with result ~p", [AtGMT,Result])),
  Result.


-spec update_raws_and_store(calendar:datetime(),string(),term()) -> ok.
update_raws_and_store(AtGMT,InDir,SSel) ->
  {{Y,M,D},{H,_,_}} = AtGMT,
  Path = lists:flatten(InDir ++ io_lib:format("/raws_ingest_~4..0B~2..0B~2..0B-~2..0B00.csv", [Y,M,D,H])),
  raws_ingest:update_now(120),
  From = fdsys_util:shift_by_seconds({{Y,M,D},{H,0,0}}, -30 * 60),
  To = fdsys_util:shift_by_seconds({{Y,M,D},{H,0,0}}, 30 * 60),
  SSel1 = raws_ingest:resolve_station_selector(SSel),
  Os = raws_ingest:retrieve_observations(SSel1,[fm10],{From,To}),
  error_logger:info_report(fdsys_util:fmt_text("rtma_server: have ~p observations for ~p, writing to ~p", [length(Os),AtGMT,Path])),
  raws_export:obs_to_csv(Os,Path),
  ok.


download_if_missing(UrlPfix,H,Vars,D) ->
  Me = self(),
  case filelib:is_dir(D) of
    true ->
      [];
    false ->
      spawn(fun () ->
              R = rtma_retr:retrieve(UrlPfix,H,Vars,D),
              Me ! {task_done, rtma_retr, R}
              end),
      rtma_retr
  end.


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
      error_logger:info_report(fdsys_util:fmt_text("fdsys_server: errors encountered during rtma_retr:retrieve, ~p", [Errs])),
      wait_for_results(lists:delete(rtma_retr,Lst),[{rtma_errors,Errs}|Res])
  end.


wait_for_results(Lst) ->
  wait_for_results(Lst,[]).


test_cycle(AtGMT) ->
  execute_cycle("http://weather.noaa.gov/pub/SL.us008001/ST.opnl/DF.gr2/DC.ndgd/GT.rtma/AR.conus",AtGMT,{station_file,"etc/raws_station_list"}).
