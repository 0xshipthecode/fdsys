-module(rtma_retr).

%% ------------------------------------------------------------------
%% API Function Exports
%% ------------------------------------------------------------------

-export([retrieve/4,is_available/3,retrieve_archived/3,retrieve_archived/2]).


-spec retrieve_archived(calendar:datetime(),[atom()],string()) -> [{error,term()}].
retrieve_archived({{Y,M,D},{H,_,_}},Vars,Dir) ->
  filelib:ensure_dir(filename:join(Dir,"fake_file")),
  UrlBase = io_lib:format("http://nomads.ncdc.noaa.gov/data/ndgd/~4..0B~2..0B/~4..0B~2..0B~2..0B", [Y,M,Y,M,D]),
  Suffix = io_lib:format("~4..0B~2..0B~2..0B~2..0B00", [Y,M,D,H]),
  Files = lists:map(fun (V) -> [var_to_archive_name(V),"_",Suffix] end, Vars),
  Rs = lists:map(fun(F) -> stream_url_to(UrlBase,Dir,F) end, Files),
  error_logger:info_report(logger_text("rtma_retr: streamed ~p files into directory ~p, converting ...", [length(Rs),Dir])),
  os:cmd("./scripts/grib_to_netcdf_archive.sh " ++ Dir ++ " > " ++ Dir ++ "/grib_to_netcdf.log"),
  error_logger:info_report(logger_text("rtma_retr: conversion process complete (~p files in directory ~p).", [length(Rs),Dir])),
  lists:filter(fun (ok) -> false; (_) -> true end, Rs).


-spec retrieve_archived(calendar:datetime(),[atom()]) -> [{error,term()}].
retrieve_archived(AtGMT={{Y,M,D},{H,_,_}},Vars) ->
  Dir = lists:flatten(io_lib:format("inputs/~4..0B~2..0B~2..0B-~2..0B00", [Y,M,D,H])),
  retrieve_archived(AtGMT,Vars,Dir).


-spec retrieve(string(),non_neg_integer(),[atom()],string()) -> [{error,term()}].
retrieve(UrlPfix,RT,Vars,Dir) ->
  filelib:ensure_dir(filename:join(Dir,"fake_file")),
  UrlBase = lists:flatten(io_lib:format("~s/RT.~2..0B", [UrlPfix, RT])),
  Files = lists:map(fun (V) -> "ds." ++ var_to_rt_name(V) ++ ".bin" end, Vars),
  Rs = lists:map(fun(F) -> stream_url_to(UrlBase,Dir,F) end, Files),
  error_logger:info_report(logger_text("rtma_retr: streamed ~p files into directory ~p, converting ...", [length(Rs),Dir])),
  os:cmd("./scripts/grib_to_netcdf_rt.sh " ++ Dir ++ " > " ++ Dir ++ "/grib_to_netcdf.log"),
  error_logger:info_report(logger_text("rtma_retr: conversion process complete (~p files in directory ~p).", [length(Rs),Dir])),
  lists:filter(fun (ok) -> false; (_) -> true end, Rs).


-spec stream_url_to(string(),string(),string()) -> ok | {error, term()}.
stream_url_to(UrlBase,Dir,F) ->
  Url = lists:flatten([UrlBase,"/",F]),
  error_logger:info_report("rtma_retr: retrieving url " ++ Url),
  case httpc:request(get, {Url, []}, [], [{stream, filename:join(Dir,F)}]) of
    {ok, saved_to_file} ->
      ErrTxt = lists:flatten(io_lib:format("rtma_retr: success retrieving URL ~s", [Url])),
      error_logger:info_report(ErrTxt),
      ok;
    {ok, Other} ->
      ErrTxt = lists:flatten(io_lib:format("rtma_retr: failure ~p retrieving URL ~s", [Other,Url])),
      error_logger:error_report(ErrTxt),
      {error, Other};
    {error, Reason} ->
      ErrTxt = lists:flatten(io_lib:format("rtma_retr: failure ~p retrieving URL ~s", [Reason,Url])),
      error_logger:error_report(ErrTxt),
      {error, Reason}
  end.


% This function checks whether all required files are up-to-date (current as of requested time.)
% Returns a list of files on which we are still waiting
%
-spec is_available(string(),non_neg_integer(),[atom()]) -> [string()].
is_available(UrlPfix,RT,Vars) ->
  error_logger:info_report(logger_text("rtma_retr: checking if RTMA available for RT ~p yet with vars=~p.", [RT,Vars])),
  UrlBase = lists:flatten(io_lib:format("~s/RT.~2..0B", [UrlPfix, RT])),
  Rs = lists:map(fun(V) -> check_file_time(UrlBase,var_to_rt_name(V)) end, Vars),
  case lists:filter(fun ([]) -> false; (_) -> true end, Rs) of
    [] ->
      error_logger:info_report(logger_text("rtma_retr: analysis for RT ~p is ready with vars=~p.",[RT,Vars])),
      [];
    Ws ->
      error_logger:info_report(logger_text("rtma_retr: analysis not ready for RT ~p, waiting on files ~p [vars=~w].", [RT,Ws,Vars])),
      Ws
  end.


-spec check_file_time(string(),string()) -> []|string().
check_file_time(UrlBase,V) ->
  Url = UrlBase ++ "/ds." ++ V ++ ".bin",
  T = calendar:universal_time(),
  Tsec = calendar:datetime_to_gregorian_seconds(T),
  case httpc:request(head, {Url, []}, [], []) of
    {ok, {{_, 200, _}, Hdr, _}} ->
      case proplists:get_value("last-modified", Hdr, undefined) of
        undefined ->
          error_logger:info_report(logger_text("rtma_retr: last-modified header missing for URL ~p.", [Url])),
          V;
        LM ->
          DT = parse_header_datetime(LM),
          Hsec = calendar:datetime_to_gregorian_seconds(DT),
	  error_logger:info_report(logger_text("rtma_retr: now - last-modified is ~p seconds in the past for URL ~p.", [Tsec - Hsec, Url])),
          % check if the last-modified is less than two hours ago [allow for delays]
          case Tsec - Hsec < 2*3600 of
            true ->
              [];
            false ->
              {V,Tsec - Hsec}
          end
      end;
    _ ->
      error_logger:info_report(logger_text("rtma_retr: HTTP request failed for URL ~p.", [Url])),
      V
  end.


% Mon, 05 May 2014 23:51:49 GMT
-spec parse_header_datetime(string()) -> calendar:datetime().
parse_header_datetime([_,_,_,$,,$ ,D1,D2,$ ,M1,M2,M3,$ ,Y1,Y2,Y3,Y4,$ ,H1,H2,$:,Min1,Min2,$:,S1,S2,$ ,$G,$M,$T]) ->
  Day = list_to_integer([D1,D2]),
  Mon = mon_to_int([M1,M2,M3]),
  Year = list_to_integer([Y1,Y2,Y3,Y4]),
  Hour = list_to_integer([H1,H2]),
  Min = list_to_integer([Min1,Min2]),
  Sec = list_to_integer([S1,S2]),
  {{Year,Mon,Day},{Hour,Min,Sec}}.


-spec mon_to_int(string()) -> non_neg_integer().
mon_to_int("Jan") -> 1;
mon_to_int("Feb") -> 2;
mon_to_int("Mar") -> 3;
mon_to_int("Apr") -> 4;
mon_to_int("May") -> 5;
mon_to_int("Jun") -> 6;
mon_to_int("Jul") -> 7;
mon_to_int("Aug") -> 8;
mon_to_int("Sep") -> 9;
mon_to_int("Oct") -> 10;
mon_to_int("Nov") -> 11;
mon_to_int("Dec") -> 12.


-spec var_to_rt_name(atom()) -> string().
var_to_rt_name(dewpoint) -> "td";
var_to_rt_name(temp) -> "temp";
var_to_rt_name(precip) -> "precipa";
var_to_rt_name(wind_spd) -> "wspd";
var_to_rt_name(wind_dir) -> "wdir".

-spec var_to_archive_name(atom()) -> string().
var_to_archive_name(dewpoint) -> "LRIA98_KWBR";
var_to_archive_name(temp) -> "LTIA98_KWBR";
var_to_archive_name(precip) -> "LEIA98_KWBR";
var_to_archive_name(wind_spd) -> "LNIA98_KWBR";
var_to_archive_name(wind_dir) -> "LNIA98_KWBR".


-spec logger_text(string(),[term()]) -> string().
logger_text(Txt,Args) ->
  lists:flatten(io_lib:format(Txt,Args)).
