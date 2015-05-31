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

-spec load_config(string()) -> [term()].
load_config(Path) ->
	{ok,V} = file:consult(Path),
	V.

