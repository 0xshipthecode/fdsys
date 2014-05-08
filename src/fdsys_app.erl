-module(fdsys_app).

-behaviour(application).

%% Application callbacks
-export([start/2, stop/1]).

%% ===================================================================
%% Application callbacks
%% ===================================================================

start(_StartType, _Args) ->
    Cfg = load_config("etc/config"),
    fdsys_sup:start_link(Cfg).


stop(_State) ->
    ok.


load_config(Path) ->
  case file:read_file(Path) of
    {ok, Cnt} ->
      case erl_scan:string(binary_to_list(Cnt)) of
        {ok,Toks,_} ->
          case erl_parse:parse_exprs(Toks) of
            {ok,AST} ->
              erl_eval:exprs(AST, erl_eval:new_bindings());
            {error,Info} ->
              ErrTxt = lists:flatten(io_lib:format("Error in eval of configuration file: ~p~n", [Info])),
              error_logger:error_report(ErrTxt),
              throw({config_eval_error, Info})
          end;
        {error,Info,Loc} ->
          ErrTxt = lists:flatten(io_lib:format("Cannot parse config file: info ~p~nlocation: ~p~n",[Info,Loc])),
          error_logger:error_report(ErrTxt),
          throw({config_parse_error,Info,Loc})
      end;
    {error,Reason} ->
      ErrTxt = lists:flatten(io_lib:format("Cannot read config file with error ~p",[Reason])),
      error_logger:error_report(ErrTxt),
      throw({config_read_error,Reason})
  end.
