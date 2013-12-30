%
% 2013-12-20 karlll <karl@ninjacontrol.com>
%

-module(compare).
-compile([export_all]).



%% -------------------------------------------------------------------------- %%
%% Dynamic Programming                                                        %%
%% -------------------------------------------------------------------------- %%

dp_change(Money,Coins) ->
	MoneySeq = lists:seq(1,Money),
	CoinSeq = lists:seq(1,length(Coins)),
	Lookup = change_lookup(MoneySeq,CoinSeq,lists:sort(Coins),[{0,0}]),
	proplists:get_value(Money,Lookup).

change_lookup([],_,_,Lookup) ->
	Lookup;

change_lookup([_M|MT],[],Denoms,Lookup) ->
	change_lookup(MT,lists:seq(1,length(Denoms)), Denoms,Lookup);

change_lookup([M|MT],[I|IT],Denoms,Lookup) ->
	Denom = lists:nth(I,Denoms),
	%io:format("M = ~p, I = ~p, Denom = ~p, M-Denom = ~p, Lookup = ~p~n",[M,I,Denom,M-Denom,Lookup]),
	NewLookup = case M of
					M when M >= Denom ->
						MinM_Denom = proplists:get_value(M-Denom,Lookup) + 1,
						MinM = proplists:get_value(M,Lookup),
						case MinM of
							undefined ->
								[{M,MinM_Denom}|Lookup];
							MinM when MinM > MinM_Denom ->
								[{M,MinM_Denom}|proplists:delete(M,Lookup)];
							_ -> Lookup
						end;
					_ ->
						Lookup
				end,
	change_lookup([M|MT],IT,Denoms,NewLookup).



%% -------------------------------------------------------------------------- %%
%% Manhattan tourist problem
%% -------------------------------------------------------------------------- %%

% Find longest path in NxM grid with weighted edges DownW & RightW
manhattan_tourist(N,M,DownW,RightW) ->
	S0 = util:new_mx(N+1,M+1),
	DownCol = [0|init_vec(util:mx_col_at(0,DownW))],
	RightRow = [0|init_vec(util:mx_row_at(0,RightW))],
	%io:format("S0 = ~p~n",[S0]),
	%io:format("D = ~p, R = ~p~n",[DownCol,RightRow]),
	S1 = util:mx_col(0,DownCol,util:mx_row(0,RightRow,S0)),
	S2 = mt(lists:seq(1,N),lists:seq(1,M),M,S1,DownW,RightW),
	util:mx_at(N,M,S2).


mt([I|IT],[],M,S,DownW,RightW) ->
	case IT of
		[] ->
			S;
		_ ->
			mt(IT,lists:seq(1,M),M,S,DownW,RightW)
	end;
	

mt([I|IT],[J|JT],M,S,DownW,RightW) ->
	
	%io:format("I=~p,J=~p~n",[I,J]),
	%io:format("S = ~p~n",[S]),
	A = util:mx_at(I-1,J,S)+util:mx_at(I-1,J,DownW),
	B = util:mx_at(I,J-1,S)+util:mx_at(I,J-1,RightW),
	V = util:v_max(A,B),
	mt([I|IT],JT,M,util:mx(I,J,V,S),DownW,RightW).




init_vec(Vect) ->
	init_vec(Vect,[]).

init_vec([],Acc) ->
	lists:reverse(Acc);
init_vec([VH|VT],[]) ->
	init_vec(VT,[VH]);
init_vec([VH|VT],Acc) ->
	init_vec(VT,[hd(Acc)+VH|Acc]).



%% -------------------------------------------------------------------------- %%
%% Longest common subsequence
%% -------------------------------------------------------------------------- %%


lcs(V,W) when is_list(V), is_list(W) ->
	S0 = util:new_mx(length(V)+1,length(W)+1), % init |V+1| x |W+1| matrix, fill w. zeroes
	B0 = util:new_mx(length(V)+1,length(W)+1,none), % Backtrack matrix
	lcs(lists:seq(1,length(V)),lists:seq(1,length(W)),V,W,S0,B0).

lcs([I|IT],[],V,W,S,B) ->
	lcs(IT,lists:seq(1,length(W)),V,W,S,B);

lcs([],_,_,_,S,B) ->
	{S,B};

lcs([I|IT],[J|JT],V,W,S,B) ->
	{SVal,BVal} = s_value(I,J,V,W,S),
	NewS = util:mx(I,J,SVal,S),
	NewB = util:mx(I,J,BVal,B),
	lcs([I|IT],JT,V,W,NewS,NewB).


s_value(I,J,V,W,S) ->
	A = util:mx_at(I-1,J,S),
	B = util:mx_at(I,J-1,S),
	C = util:mx_at(I-1,J-1,S),
	VI = lists:nth(I,V),
	VJ = lists:nth(J,W),
	SVal = case VI of
		VI when VI =:= VJ ->
			util:v_max(A,B,C+1);
		_ ->
			util:v_max(A,B)
		end,
	case SVal of
		SVal when SVal =:= A ->
			{SVal,deletion};
		SVal when SVal =:= B ->
			{SVal,insertion};
		SVal when SVal =:= C+1 ->
			{SVal,match}
	end.

output_lcs(BackT,V,I,J) ->
	output_lcs(BackT,V,I,J,[]).

output_lcs(BackT,V,I,J,Acc) when I =:= 0; J =:= 0 ->
	Acc;
output_lcs(BackT,V,I,J,Acc) ->
	case util:mx_at(I,J,BackT) of
		deletion ->
			output_lcs(BackT,V,I-1,J,Acc);
		insertion ->
			output_lcs(BackT,V,I,J-1,Acc);
		match ->
			output_lcs(BackT,V,I-1,J-1,[lists:nth(I,V)|Acc])
	end.



%% -------------------------------------------------------------------------- %%
%% Longest path in DAG 
%% -------------------------------------------------------------------------- %%

% "1 -> 2:4" -> {NodeA,NodeB,EdgeWeight}
get_adj_list3(AdjStrings) ->
	get_adj_list3(AdjStrings,[]).
get_adj_list3([],Acc) ->
	Acc;
get_adj_list3([AdjStr|T],Acc) ->
	Sep = " -> ",
	NodeA = string:sub_string(AdjStr,1,string:str(AdjStr,Sep)-1),
	NodeB_W = string:tokens(
				string:sub_string(AdjStr,string:str(AdjStr,Sep)+length(Sep),length(AdjStr)),
				":"),
	get_adj_list3(T,[{NodeA,hd(NodeB_W),hd(tl(NodeB_W))}|Acc]). 


print_longest_path(AdjList,Start,End) ->
	G = build_dag(AdjList),
	{Weight,Path} = get_longest_path(G,Start,End),
	io:format("~p~n~s~n",[Weight,util:print_path(Path)]).
	


% {Weight,NodesInPath}
get_longest_path(G,Start,End) ->
	SinkReach = digraph_utils:reaching([End],G),
	SourceReach =digraph_utils:reachable([Start],G),
	Context = lists:filter(fun(E)-> lists:member(E,SinkReach) end, SourceReach),
	SubG = digraph_utils:subgraph(G,Context,[{keep_labels,true}]),
	SubTS = digraph_utils:topsort(SubG),
	%io:format("SinkReach = ~p~nSourceReach = ~p~nContext = ~p~nSubTS = ~p~n",[SinkReach,SourceReach,Context,SubTS]),
	PathLen = get_lp_and_prune(SubG,SubTS,[]),
	{proplists:get_value(End,PathLen),digraph:get_path(SubG,Start,End)}.
	


get_lp_and_prune(G,[],Acc) ->
	Acc;

get_lp_and_prune(G,[N|T],Acc) ->
	WEdges = weighted_edges(G,Acc,N),
	NodeWeight = case WEdges of
			[] -> % zero in-degree
				{N,0};
		 	_ ->
				{PredMaxWeight,PredMaxEdge} = proplists:lookup(lists:max(proplists:get_keys(WEdges)),WEdges),
				ok = del_edges(G,lists:map(fun(E)->{_K,V} = E, V end, WEdges),[PredMaxEdge]),
				{N,PredMaxWeight}
			end,
	get_lp_and_prune(G,T,[NodeWeight|Acc]).


weighted_edges(G,WL,Node) ->
	Edges = digraph:in_edges(G,Node),
	PredWeights = lists:map(fun(E) ->
		{Edge,Pred,Node,WeightStr} = digraph:edge(G,E),
		PredWeight = proplists:get_value(Pred,WL),
		{Weight,_} = string:to_integer(WeightStr),
		%io:format("Edge = ~p, Pred = ~p, Node = ~p, WeightStr = ~p, PredWeight = ~p, Weight = ~p~n",[Edge,Pred,Node,WeightStr,PredWeight,Weight]),
		{PredWeight+Weight,Edge}
	end, Edges).



del_edges(G,Edges,Exempts) ->
	lists:foreach(fun(E) ->
			case lists:member(E,Exempts) of 
				false ->
					%io:format("Removing edge = ~p~n",[E]),
					digraph:del_edge(G,E),
					ok;
				true ->
					ok
			end
		end,
			Edges).




build_dag(AdjStrings) ->
	AdjList = get_adj_list3(AdjStrings),
	% create empty graph
	G = digraph:new(),
	% create all vertices
	Ns = get_nodes(AdjList),
	lists:foreach(fun(N) -> digraph:add_vertex(G,N) end, Ns),
	lists:foreach(fun(E) -> 
		{NodeA,NodeB,EdgeW} = E,
		io:format("Adding edge ~p -> ~p with weight ~p~n",[NodeA,NodeB,EdgeW]),
		digraph:add_edge(G,NodeA,NodeB,EdgeW)
		end,
		AdjList),
	G.

% Adjlist = [{NodeA,NodeB,EdgeWeight}]
get_nodes(AdjList) ->
	get_nodes(AdjList,[]).
get_nodes([],Acc) ->
	sets:to_list(sets:from_list(Acc));
get_nodes([NodeEntry|T],Acc) ->
	{NodeA,NodeB,_W} = NodeEntry,
	get_nodes(T,[NodeA,NodeB]++Acc).


%% -------------------------------------------------------------------------- %%
%% Global alignment graph
%% -------------------------------------------------------------------------- %%

print_ga(V,W,IndelPenalty) ->
		{S,B} = global_align(V,W,IndelPenalty),
		io:format("~p~n",[util:mx_at(length(V),length(W),S)]),
		{VA,WA} = output_ga(B,V,W,length(V),length(W)),
		io:format("~s~n~s~n",[VA,WA]).

global_align(V,W,IndelPenalty) when is_list(V), is_list(W) ->
	S0 = util:new_mx(length(V)+1,length(W)+1), % init |V+1| x |W+1| matrix,
	B0 = util:new_mx(length(V)+1,length(W)+1,x), % Backtrack matrix
	Init = fun(I,L) -> lists:map(fun(E)->E*I end, lists:seq(0,L)) end,
	S1 = util:mx_col(0,Init(IndelPenalty,length(V)),S0), % init col 0
	S2 = util:mx_row(0,Init(IndelPenalty,length(W)),S1), % init row 0
	B1 = util:mx_col(0,lists:duplicate(length(V)+1,d),B0), % B(i,0) = deletion
	B2 = util:mx_row(0,lists:duplicate(length(W)+1,i),B1),  % B(0,j) = insertion
	B3 = util:mx(0,0,x,B2), % B(0,0) = 0
	global_align(lists:seq(1,length(V)),lists:seq(1,length(W)),V,W,IndelPenalty,S2,B3).

global_align([I|IT],[],V,W,IndelP,S,B) ->
	global_align(IT,lists:seq(1,length(W)),V,W,IndelP,S,B);

global_align([],_,_,_,_,S,B) ->
	{S,B};

global_align([I|IT],[J|JT],V,W,IndelP,S,B) ->
	{SVal,BVal} = g_value(I,J,V,W,IndelP,S),
	NewS = util:mx(I,J,SVal,S),
	NewB = util:mx(I,J,BVal,B),
	global_align([I|IT],JT,V,W,IndelP,NewS,NewB).


g_value(I,J,V,W,IndelP,S) ->
	VI = lists:nth(I,V),
	VJ = lists:nth(J,W),
	A = util:mx_at(I-1,J,S)+IndelP,
	B = util:mx_at(I,J-1,S)+IndelP,
	C = util:mx_at(I-1,J-1,S)+util:mx_at(b62(VI),b62(VJ),blosum62()),
	SVal = util:v_max(A,B,C),
	case SVal of
		SVal when SVal =:= A ->
			{SVal,d};
		SVal when SVal =:= B ->
			{SVal,i};
		SVal when SVal =:= C ->
			{SVal,m}
	end.

output_ga(BackT,V,W,I,J) ->
	output_ga(BackT,V,W,I,J,[],[]).

output_ga(BackT,V,W,I,J,AccV,AccW) when I =:= 0, J =:= 0 ->
	{AccV,AccW};
output_ga(BackT,V,W,I,J,AccV,AccW) ->
	
	case util:mx_at(I,J,BackT) of
		d ->
			output_ga(BackT,V,W,I-1,J,[lists:nth(I,V)|AccV],[$-|AccW]);
		i ->
			output_ga(BackT,V,W,I,J-1,[$-|AccV],[lists:nth(J,W)|AccW]);
		m ->
			output_ga(BackT,V,W,I-1,J-1,[lists:nth(I,V)|AccV],[lists:nth(J,W)|AccW])
	end.





%
% Scoring matrices
%

blosum62() ->[[ 4,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -2, -1, -1, -1,  1,  0,  0, -3, -2],
			  [ 0,  9, -3, -4, -2, -3, -3, -1, -3, -1, -1, -3, -3, -3, -3, -1, -1, -1, -2, -2],
			  [-2, -3,  6,  2, -3, -1, -1, -3, -1, -4, -3,  1, -1,  0, -2,  0, -1, -3, -4, -3],
			  [-1, -4,  2,  5, -3, -2,  0, -3,  1, -3, -2,  0, -1,  2,  0,  0, -1, -2, -3, -2],
			  [-2, -2, -3, -3,  6, -3, -1,  0, -3,  0,  0, -3, -4, -3, -3, -2, -2, -1,  1,  3],
			  [ 0, -3, -1, -2, -3,  6, -2, -4, -2, -4, -3,  0, -2, -2, -2,  0, -2, -3, -2, -3],
			  [-2, -3, -1,  0, -1, -2,  8, -3, -1, -3, -2,  1, -2,  0,  0, -1, -2, -3, -2,  2],
			  [-1, -1, -3, -3,  0, -4, -3,  4, -3,  2,  1, -3, -3, -3, -3, -2, -1,  3, -3, -1],
			  [-1, -3, -1,  1, -3, -2, -1, -3,  5, -2, -1,  0, -1,  1,  2,  0, -1, -2, -3, -2],
			  [-1, -1, -4, -3,  0, -4, -3,  2, -2,  4,  2, -3, -3, -2, -2, -2, -1,  1, -2, -1],
			  [-1, -1, -3, -2,  0, -3, -2,  1, -1,  2,  5, -2, -2,  0, -1, -1, -1,  1, -1, -1],
			  [-2, -3,  1,  0, -3,  0,  1, -3,  0, -3, -2,  6, -2,  0,  0,  1,  0, -3, -4, -2],
			  [-1, -3, -1, -1, -4, -2, -2, -3, -1, -3, -2, -2,  7, -1, -2, -1, -1, -2, -4, -3],
			  [-1, -3,  0,  2, -3, -2,  0, -3,  1, -2,  0,  0, -1,  5,  1,  0, -1, -2, -2, -1],
			  [-1, -3, -2,  0, -3, -2,  0, -3,  2, -2, -1,  0, -2,  1,  5, -1, -1, -3, -3, -2],
			  [ 1, -1,  0,  0, -2,  0, -1, -2,  0, -2, -1,  1, -1,  0, -1,  4,  1, -2, -3, -2],
			  [ 0, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1,  0, -1, -1, -1,  1,  5,  0, -2, -2],
			  [ 0, -1, -3, -2, -1, -3, -3,  3, -2,  1,  1, -3, -2, -2, -3, -2,  0,  4, -3, -1],
			  [-3, -2, -4, -3,  1, -2, -2, -3, -3, -2, -1, -4, -4, -2, -3, -3, -2, -3, 11,  2],
			  [-2, -2, -3, -2,  3, -3,  2, -1, -2, -1, -1, -2, -3, -1, -2, -2, -2, -1,  2,  7]].

% short map
b62(X) -> blosum62_rc_map(X).

% map protein alphabet to blosum62 cost matrix row/col
blosum62_rc_map($A) -> 0;
blosum62_rc_map($C) -> 1;
blosum62_rc_map($D) -> 2;
blosum62_rc_map($E) -> 3;
blosum62_rc_map($F) -> 4;
blosum62_rc_map($G) -> 5;
blosum62_rc_map($H) -> 6;
blosum62_rc_map($I) -> 7;
blosum62_rc_map($K) -> 8;
blosum62_rc_map($L) -> 9;
blosum62_rc_map($M) -> 10;
blosum62_rc_map($N) -> 11;
blosum62_rc_map($P) -> 12;
blosum62_rc_map($Q) -> 13;
blosum62_rc_map($R) -> 14;
blosum62_rc_map($S) -> 15;
blosum62_rc_map($T) -> 16;
blosum62_rc_map($V) -> 17;
blosum62_rc_map($W) -> 18;
blosum62_rc_map($Y) -> 19.










