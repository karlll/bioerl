%
% 2013-12-03 karlll <karl@ninjacontrol.com>
%

%
% Genome assembly
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


 test_debruijn_graph() ->
 	L = get_dbgraph_adj_list(util:get_input("in.txt")),
 	util:write_result_text(L,"res.txt").
 

get_dbgraph_adj_list(Kmers) ->
	AdjList = lists:keysort(1,debruijn_graph(Kmers)),
	lists:map(fun(F)-> 
		{Kmer,L} = F,
		io_lib:format("~s -> ~s~n",[Kmer,string:join(L,",")])
	end, AdjList).


get_dbgraph_adj_list(String,K) ->
	AdjList1 = lists:keysort(1,debruijn_graph(String, K)),
	AdjList2 = lists:filter(fun(E) -> {_K,L} = E, 
									  case L of 
									  	[] -> false;
									  	_ -> true
									  end
									end, AdjList1),
	lists:map(fun(F)-> 
		{Kmer,L} = F,
		io_lib:format("~s -> ~s~n",[Kmer,string:join(lists:sort(L),",")])
	end, AdjList2).

% Construct a debruijn graph from a list of K-mers
debruijn_graph(Kmers) ->
	debruijn_graph(Kmers,[]).

debruijn_graph([],Acc) ->
	Acc;

debruijn_graph([Kmer|T],Acc) ->
	P = prefix(Kmer),
	S = suffix(Kmer),
	NewAcc = case lists:keyfind(P,1,Acc) of
			{P,S2} ->
				lists:keyreplace(P,1,Acc,{P,[S|S2]});
			false ->
				[{P,[S]}|Acc]
			end,
	debruijn_graph(T,NewAcc);	

% construct a debruijn graph by dividing the provided String in K-mers
debruijn_graph(String,K) when is_integer(K) ->
	Nmers = motif:get_kmers(String,K-1),
	lists:map(fun(E)-> filter_nodes(E,String) end, overlapping_graph(Nmers)).

% remove nodes that does not correspoing to a edge k-mer in the original strings
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
%% Eulerian path                                                              %%
%% -------------------------------------------------------------------------- %%


% Read an adjecency list from in.txt, output an eulerian path to res.txt
% Precond: Graph has an eulerian path.
test_eulerian_path() ->
	G = build_graph(util:get_input("in.txt")),
	util:write_result_text(string:join(find_path(G,path),"->"),"res.txt").

% Get the directed graph defined by AdjStrings, ["Node1->Node2,Node3"]
% NB:  digraph:delete(G) should be called for the graph G returned by this function.
build_graph(AdjStrings) ->
	AdjList = get_adj_list2(AdjStrings),
	% create empty graph
	G = digraph:new(),
	% create all vertices
	Ns = get_nodes(AdjList),
	VtxList = lists:map(fun(N) -> add_vertex(G,N) end, Ns),
	% add edges from adj list to G
	lists:foreach(fun(NodeEntry) ->
					{N,NL} = NodeEntry,
					lists:foreach(fun(E)-> 
						io:format("Adding edge from node ~s to node ~s~n",[N,E]),
						add_edge(G,VtxList,N,E) 
						end, NL) end,
					AdjList),
	G.


%% Get a adjencency list from a set of strings in the form "Node -> NeighborA[,NeighborB]"
get_adj_list2(AdjStrings) ->
	get_adj_list2(AdjStrings,[]).

get_adj_list2([],Acc) ->
	Acc;
get_adj_list2([AdjStr|T],Acc) ->
	get_adj_list2(T,[get_nbs(AdjStr)|Acc]).

get_nodes(AdjList) ->
	get_nodes(AdjList,[]).
get_nodes([],Acc) ->
	sets:to_list(sets:from_list(Acc));
get_nodes([NodeEntry|T],Acc) ->
	{Node,NeighborList} = NodeEntry,
	get_nodes(T,[Node|NeighborList]++Acc).
	


% Get neighbors for a adjecency string in the format "Node -> Neighbors"
% returns {Node,NeighborList}
get_nbs(AdjStr) ->
	Sep = " -> ",
	NodeA = string:sub_string(AdjStr,1,string:str(AdjStr,Sep)-1),
	NodeBLst = string:sub_string(AdjStr,string:str(AdjStr,Sep)+length(Sep),length(AdjStr)),
	{NodeA, string:tokens(NodeBLst,",")}.


% Add edge to directed graph G
add_edge(G,VertexList,VertexName1,VertexName2) ->
	V1 = lists:keyfind(VertexName1,1,VertexList),
	V2 = lists:keyfind(VertexName2,1,VertexList),
	digraph:add_edge(G,VertexName1,VertexName2),
	ok.

% Add verted to directed graph G
add_vertex(G,VertexName) ->
	io:format("Adding vtx = ~s~n",[VertexName]),
	{VertexName,digraph:add_vertex(G,VertexName)}.


% Find an eulerian path or cycle in digraph G. Type = path|cycle
find_path(G,Type) ->
	S = push(new_stack(),select_startnode(G,Type)),
	io:format("Initial stack = ~p~n",[S]),
	find_path(G,S,[]).
	

find_path(Graph,[],Acc) ->
	digraph:delete(Graph),
	Acc;

find_path(Graph,Stack,Acc) ->
	Node = peek(Stack),
	io:format("Current node = ~p~n",[Node]),
	{NewStack,NewAcc} = case digraph:out_edges(Graph,Node) of
		[E|T] ->
			{E,Node,NewDest,_} = digraph:edge(Graph,E), 
			true = digraph:del_edge(Graph,E),
			io:format("Pushing node = ~p, Deleted edge = ~p~n",[NewDest,E]),
			{push(Stack,NewDest),Acc};
		[] ->
			io:format("Adding/Popping ~p~n",[Node]),
			{Node,NS} = pop(Stack),
			{NS,[Node|Acc]}
	end,
	find_path(Graph,NewStack,NewAcc).

new_stack() -> [].
push(S,E) -> [E|S].
pop([H|S]) -> {H,S}.
peek([H|S]) -> H. 

random_element(List) ->
	lists:nth(random:uniform(length(List)),List).

% G is a digraph which has an Eulerian path
select_startnode(G,path) ->
	Nodes = digraph:vertices(G),
	StartNodes = lists:filter(fun(N) -> D = digraph:out_degree(G,N) - digraph:in_degree(G,N),
						case D of
							1 ->
								true;
							_ -> 
								false
						end
					end,Nodes),
	hd(StartNodes);

% G is a digraph which has an Eulerian cycle
select_startnode(G,cycle) ->
	{X,Y,Z} = now(),
	random:seed(X,Y,Z),
	Nodes = digraph:vertices(G),
	lists:nth(random:uniform(length(Nodes)),Nodes).



