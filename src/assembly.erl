%
% 2013-12-03 karlll <karl@ninjacontrol.com>
%

-module(assembly).
-compile([export_all]).



%% -------------------------------------------------------------------------- %%
%% String composition                                                         %%
%% -------------------------------------------------------------------------- %%

print_str_compo(String,K) ->
	StrCompo = str_compo(String,K),
	lists:foreach(fun(El) -> io:format("~s~n",[El]) end, StrCompo).

get_str_compo(String,K) ->
	StrCompo = str_compo(String,K),
	lists:map(fun(El) -> io_lib:format("~s~n",[El]) end, StrCompo).


str_compo(String,K) ->
	Kmers = motif:get_kmers(String,K),
	lists:sort(Kmers).



%% -------------------------------------------------------------------------- %%
%% Overlapping graph                                                          %%
%% -------------------------------------------------------------------------- %%




write_olgraph_adj_list(Kmers,Filename) ->
	{ok, File} = file:open(Filename,[write]),
	AdjList1 = lists:keysort(1,overlapping_graph(Kmers)),
	AdjList2 = lists:filter(fun(E) -> {_K,L} = E, 
									  case L of 
									  	[] -> false;
									  	_ -> true
									  end
									end, AdjList1),
	AdjList2,
	lists:foreach(fun(F)-> 
		{Kmer,L} = F,
		file:write(File,io_lib:format("~s -> ",[Kmer])),
		lists:foreach(fun(G) -> file:write(File,io_lib:format("~s",[G])) end, L),
		file:write(File,io_lib:format("~n",[])) 
	end, AdjList2).


print_olgraph_adj_list(Kmers) ->
	AdjList1 = lists:keysort(1,overlapping_graph(Kmers)),
	AdjList2 = lists:filter(fun(E) -> {_K,L} = E, 
									  case L of 
									  	[] -> false;
									  	_ -> true
									  end
									end, AdjList1),
	AdjList2,
	lists:foreach(fun(F)-> 
		{Kmer,L} = F,
		io:format("~s -> ",[Kmer]),
		lists:foreach(fun(G) -> io:format("~s ",[G]) end, L),
		io:format("~n")
	end, AdjList2).

overlapping_graph(Kmers) ->
	lists:map(fun(Kmer) -> 

			{Kmer, get_adj_list(Kmer,Kmers)} end,

		Kmers).

get_adj_list(Kmer,Kmers) ->
	Kmers2 = lists:delete(Kmer,Kmers),
	get_adj_list(Kmer,Kmers2,[]).

get_adj_list(_Kmer,[],Acc) ->
	Acc;
get_adj_list(Kmer1,[Kmer2|Tail],Acc) ->
	NewAcc = case is_overlapping(Kmer1,Kmer2) of
		true ->
			[Kmer2|Acc];
		false ->
			Acc
		end,
	get_adj_list(Kmer1,Tail,NewAcc).

is_overlapping(Kmer1,Kmer2) ->
	S = suffix(Kmer1),
	P = prefix(Kmer2),
	case S of 
		P ->
			true;
		_ ->
			false
	end.


prefix(Kmer) ->
	lists:sublist(Kmer,length(Kmer)-1).
suffix([H|Tail]) ->
	Tail.


%% -------------------------------------------------------------------------- %%
%% DeBruijn graph                                                            %%
%% -------------------------------------------------------------------------- %%

get_dbgraph_adj_list(String,K) ->
	AdjList1 = lists:keysort(1,debruijn_graph(String, K)),
	AdjList2 = lists:filter(fun(E) -> {_K,L} = E, 
									  case L of 
									  	[] -> false;
									  	_ -> true
									  end
									end, AdjList1),
	AdjList2,
	lists:map(fun(F)-> 
		{Kmer,L} = F,
		io_lib:format("~s -> ~s~n",[Kmer,string:join(L,",")])
	end, AdjList2).

debruijn_graph(String,K) ->
	Nmers = motif:get_kmers(String,K-1),
	lists:map(fun(E)-> filter_nodes(E,String) end, overlapping_graph(Nmers)).

% remove nodes that does not correspoing to a edge kmer in the original strings
filter_nodes(E,String) ->
	{N1,NL} = E,
	NewNL = lists:filter(fun(N) ->
				EdgeKmer = get_edge_kmer(N1,N),
				case string:str(String,EdgeKmer) of
					0 ->
						false;
					_ -> 
						true
				end
			end,
			NL),
	{N1,lists:sort(NewNL)}.

% ex.: AAGATTGAA -> AGATTGAAG => edge k-mer AAGATTGAAG
get_edge_kmer([H1|T1],N2) ->
	[H1|N2].
	
%% -------------------------------------------------------------------------- %%
%% Util                                                                       %%
%% -------------------------------------------------------------------------- %%


get_input(Filename) ->
	{ok, F} = file:open(Filename,[read]),
	L = read_input(F),
	get_input(L,[]).

get_input(["Input:\n"|T],Acc) ->
	get_input(T,Acc);
get_input(["Output:\n"|T],Acc) ->
	Acc;
get_input([],Acc) ->
	Acc;
get_input([L|T], Acc) ->
	get_input(T,[string:strip(L,right,$\n)|Acc]).


read_input(File) ->
    case file:read_line(File) of
        eof        -> [];
        {ok, Line} -> [Line | read_input(File)]
    end.

