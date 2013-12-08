
%
% 2013-12-03 karlll <karl@ninjacontrol.com>
%


-module(util).

-compile([export_all]).

%% -------------------------------------------------------------------------- %%
%% Utility functions                                                          %%
%% -------------------------------------------------------------------------- %%

write_result(Result,Filename) ->
	file:write_file(Filename, io_lib:fwrite("~p.\n", [Result])).

write_result_text(Result,Filename) ->
	file:write_file(Filename, io_lib:fwrite("~s\n", [Result])).

load_string(FileName) ->
	V = load_binary_string(FileName),
	L = binary_to_list(V),
	string:strip(L,both,$\n).

load_binary_string(FileName) ->
	{ok, BinaryString} = file:read_file(FileName),
	BinaryString.
	







