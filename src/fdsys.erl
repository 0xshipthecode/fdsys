
-module(fdsys).
-author("Martin Vejmelka <vejmelkam@gmail.com>").
-export([start/0]).


start() ->
    application:start(inets),
    raws_ingest:start(),
    afm_ingest:start(),
    ok = application:start(fdsys).
