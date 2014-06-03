-module(fdsys_server).
-behaviour(gen_server).
-define(SERVER, ?MODULE).

%% ------------------------------------------------------------------
%% API Function Exports
%% ------------------------------------------------------------------

-export([start_link/1]).
-export([test_cycle/1]).

%% ------------------------------------------------------------------
%% gen_server Function Exports
%% ------------------------------------------------------------------

-export([init/1, handle_call/3, handle_cast/2, handle_info/2,
         terminate/2, code_change/3]).

%% ------------------------------------------------------------------
%% API Function Definitions
%% ------------------------------------------------------------------

start_link(Cfg) ->
    gen_server:start_link({local, ?SERVER}, ?MODULE, Cfg, []).

%% ------------------------------------------------------------------
%% gen_server Function Definitions
%% ------------------------------------------------------------------

init(Args) ->
    schedule_next_run(?SERVER,0),
    {ok, Args}.

handle_call(_Request, _From, State) ->
    {reply, ok, State}.

handle_cast(_Msg, State) ->
    {noreply, State}.

handle_info(Info, State) ->
  case Info of
    cycle_nowcast ->
      Pfix = proplists:get_value(url_prefix,State),
      AtGMT = calendar:universal_time(),
      SSel = proplists:get_value(station_selector,State),
      execute_cycle(Pfix,AtGMT,SSel),
      schedule_next_run(?SERVER, 40),
      {noreply, State};
    _ ->
      {noreply, State}
  end.

terminate(_Reason, _State) ->
    ok.

code_change(_OldVsn, State, _Extra) ->
    {ok, State}.

%% ------------------------------------------------------------------
%% Internal Function Definitions
%% ------------------------------------------------------------------

-spec schedule_next_run(pid(),non_neg_integer()) -> ok.
schedule_next_run(PID,WaitAtLeastMins) ->
  {{_,_,_},{_,Min,_}} = calendar:universal_time(),
  TimeoutS = case Min < 20 of
    true ->
      % next run will be H now
      error_logger:info_msg("fdsys_server: scheduling next run in ~p minutes.", [20-Min+WaitAtLeastMins]),
      (20 - Min) * 60;
    false ->
      % next run will be NOW
      error_logger:info_msg("fdsys_server: running a cycle immediately (but after ~p mins).", [WaitAtLeastMins]),
      0
  end,
  timer:send_after((WaitAtLeastMins * 60 + TimeoutS) * 1000, PID, cycle_nowcast),
  ok.


-spec execute_cycle(string(),calendar:datetime(),raws_ingest:station_selector()) -> [term()].
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
      error_logger:info_report(fdsys_util:fmt_text("fdsys_server: RTMA not ready yet, waiting on ~p", [WaitingOnFiles])),
      timer:sleep(120*1000),
      execute_cycle(UrlPfix,AtGMT,SSel)
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
  DlTasks = lists:flatten(lists:map(fun ({RT,Dir}) -> download_rtma_if_missing(UrlPfix,RT,Vars,Dir) end, [{H,InDir},{H0,InDir0}])),

  % wait until everything is downloaded and stored
  Result = case wait_for_results([raws_ingest|DlTasks]) of
    ok ->
      % phase 3: run the fmda
      LogFile = lists:flatten(io_lib:format("fmncast-~4..0B~2..0B~2..0B-~2..0B00.log",[Y,M,D,H])),
      error_logger:info_report(fdsys_util:fmt_text("fdsys_server: running fmncast.py for GMT time ~p", [AtGMT])),
      os:cmd(io_lib:format("python psrc/fmncast.py ~s ~s state &> outputs/~s", [InDir0,InDir,LogFile])),
      ok;
    Errors ->
      Errors
  end,

  % execute postprocessing step
  LogFile2 = lists:flatten(io_lib:format("postproc-~4..0B~2..0B~2..0B-~2..0B00.log",[Y,M,D,H])),
  os:cmd(io_lib:format("python psrc/postproc.py ~p ~p ~p ~p &> outputs/~s", [Y,M,D,H,LogFile2])),

  % log result
  error_logger:info_report(fdsys_util:fmt_text("fdsys_server: execute_cycle completed at GMT ~p with result ~p", [AtGMT,Result])),
  Result.


-spec download_rtma_if_missing(string(),non_neg_integer(),[atom()],string()) -> []|rtma_retr.
download_rtma_if_missing(UrlPfix,H,Vars,D) ->
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


-spec wait_for_results([term()],[term()]) -> ok|[term()].
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


-spec wait_for_results([term()]) -> ok|[term()].
wait_for_results(Lst) ->
  wait_for_results(Lst,[]).

-spec test_cycle(calendar:datetime()) -> [term()].
test_cycle(AtGMT) ->
  execute_cycle("http://weather.noaa.gov/pub/SL.us008001/ST.opnl/DF.gr2/DC.ndgd/GT.rtma/AR.conus",AtGMT,{station_file,"etc/raws_station_list"}).
