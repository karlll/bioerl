%
% 2013-12-20 karlll <karl@ninjacontrol.com>
%

%
% Comparing sequences
% 

-module(compare).
-compile([export_all]).


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
		{PredWeight+Weight,Edge}
	end, Edges).



del_edges(G,Edges,Exempts) ->
	lists:foreach(fun(E) ->
			case lists:member(E,Exempts) of 
				false ->
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



%% -------------------------------------------------------------------------- %%
%% Local Alignment
%% -------------------------------------------------------------------------- %%

print_la(V,W,IndelPenalty) ->
		{S,B} = local_align(V,W,IndelPenalty),
		{Max,SI,SJ} = util:mx_max(S),
		io:format("~p~n",[util:mx_at(SI,SJ,S)]),
		{VA,WA} = output_la(B,V,W,SI,SJ),
		io:format("~s~n~s~n",[VA,WA]).

local_align(V,W,IndelPenalty) when is_list(V), is_list(W) ->
	S0 = util:new_mx(length(V)+1,length(W)+1), % init |V+1| x |W+1| matrix,
	B0 = util:new_mx(length(V)+1,length(W)+1,x), % Backtrack matrix
	S1 = util:mx_col(0,lists:duplicate(length(V)+1,0),S0), % init col 0
	S2 = util:mx_row(0,lists:duplicate(length(W)+1,0),S1), % init row 0
	local_align(lists:seq(1,length(V)),lists:seq(1,length(W)),V,W,IndelPenalty,S2,B0).

local_align([I|IT],[],V,W,IndelP,S,B) ->
	local_align(IT,lists:seq(1,length(W)),V,W,IndelP,S,B);

local_align([],_,_,_,_,S,B) ->
	{S,B};

local_align([I|IT],[J|JT],V,W,IndelP,S,B) ->
	{SVal,BVal} = l_value(I,J,V,W,IndelP,S),
	NewS = util:mx(I,J,SVal,S),
	NewB = util:mx(I,J,BVal,B),
	local_align([I|IT],JT,V,W,IndelP,NewS,NewB).


l_value(I,J,V,W,IndelP,S) ->
	VI = lists:nth(I,V),
	VJ = lists:nth(J,W),
	A = util:mx_at(I-1,J,S)+IndelP,
	B = util:mx_at(I,J-1,S)+IndelP,
	C = util:mx_at(I-1,J-1,S)+util:mx_at(p250(VI),p250(VJ),pam250()),
	SVal = lists:max([0,A,B,C]),
	case SVal of
		SVal when SVal =:= 0 ->
			{SVal,x};
		SVal when SVal =:= A ->
			{SVal,d};
		SVal when SVal =:= B ->
			{SVal,i};
		SVal when SVal =:= C ->
			{SVal,m}
	end.

output_la(BackT,V,W,I,J) ->
	output_la(BackT,V,W,I,J,[],[]).

output_la(BackT,V,W,I,J,AccV,AccW) when I =:= 0, J =:= 0 ->
	{AccV,AccW};
output_la(BackT,V,W,I,J,AccV,AccW) ->
	
	case util:mx_at(I,J,BackT) of
		x ->
			{AccV,AccW};
		d ->
			output_la(BackT,V,W,I-1,J,[lists:nth(I,V)|AccV],[$-|AccW]);
		i ->
			output_la(BackT,V,W,I,J-1,[$-|AccV],[lists:nth(J,W)|AccW]);
		m ->
			output_la(BackT,V,W,I-1,J-1,[lists:nth(I,V)|AccV],[lists:nth(J,W)|AccW])
	end.


%% -------------------------------------------------------------------------- %%
%% Levenshtein Distance
%% -------------------------------------------------------------------------- %%

edit_distance(V,W) ->
	S0 = util:new_mx(length(V)+1,length(W)+1), % init |V+1| x |W+1| matrix,
	S1 = util:mx_col(0,lists:seq(0,length(V)),S0), % init col 0
	S2 = util:mx_row(0,lists:seq(0,length(W)),S1),% init row 0
	S3 = ed(lists:seq(1,length(V)),lists:seq(1,length(W)),V,W,S2),
	util:mx_at(length(V),length(W),S3).

ed([I|IT],[],V,W,S) ->
	case IT of
		[] ->
			S;
		_ ->
		 	ed(IT,lists:seq(1,length(W)),V,W,S)
	end;
ed([I|IT],[J|JH],V,W,S) ->
	VI = lists:nth(I,V),
	WJ =  lists:nth(J,W),
	Val = case VI of
			WJ ->
				util:mx_at(I-1,J-1,S);
			_ ->
				A = util:mx_at(I-1,J,S) +1,
				B = util:mx_at(I,J-1,S) +1,
				C = util:mx_at(I-1,J-1,S)+1,
				lists:min([A,B,C])
			end,
	ed([I|IT],JH,V,W,util:mx(I,J,Val,S)).


%% -------------------------------------------------------------------------- %%
%% Fitting Alignment
%% -------------------------------------------------------------------------- %%

print_fa(V,W,IndelPenalty) ->
		{S,B} = fitting_align(V,W,IndelPenalty),
		SCol = length(W), % rightmost column
		{Max,SRow} = util:mx_max_vector(util:mx_col_at(SCol,S)), % Row with maximum value in rightmost col
		io:format("S = ~p ~n",[S]),
		io:format("Max = ~p SCol = ~p SRow = ~p ~n",[Max,SCol,SRow]),
		{VA,WA} = output_fa(B,V,W,SRow,SCol),
		
		io:format("~p~n~s~n~s~n",[Max,VA,WA]).

fitting_align(V,W,IndelPenalty) when is_list(V), is_list(W) ->
	S0 = util:new_mx(length(V)+1,length(W)+1), % init |V+1| x |W+1| matrix,
	B0 = util:new_mx(length(V)+1,length(W)+1,x), % Backtrack matrix
	S1 = util:mx_col(0,lists:duplicate(length(V)+1,0),S0), % init col 0
	S2 = util:mx_row(0,lists:seq(0,-(length(W)),-1),S1), % init row 0
	fitting_align(lists:seq(1,length(V)),lists:seq(1,length(W)),V,W,IndelPenalty,S2,B0).

fitting_align([I|IT],[],V,W,IndelP,S,B) ->
	fitting_align(IT,lists:seq(1,length(W)),V,W,IndelP,S,B);

fitting_align([],_,_,_,_,S,B) ->
	{S,B};

fitting_align([I|IT],[J|JT],V,W,IndelP,S,B) ->
	{SVal,BVal} = f_value(I,J,V,W,IndelP,S),
	NewS = util:mx(I,J,SVal,S),
	NewB = util:mx(I,J,BVal,B),
	fitting_align([I|IT],JT,V,W,IndelP,NewS,NewB).


f_value(I,J,V,W,IndelP,S) ->
	VI = lists:nth(I,V),
	VJ = lists:nth(J,W),
	Score = case VI of
		VJ ->
			1;
		_ ->
			-1
		end,
	A = util:mx_at(I-1,J,S)+IndelP,
	B = util:mx_at(I,J-1,S)+IndelP,
	C = util:mx_at(I-1,J-1,S)+Score,
	SVal = util:v_max(A,B,C),
	case SVal of
		SVal when SVal =:= A ->
			{SVal,d};
		SVal when SVal =:= B ->
			{SVal,i};
		SVal when SVal =:= C ->
			{SVal,m}
	end.

output_fa(BackT,V,W,I,J) ->
	output_fa(BackT,V,W,I,J,[],[]).

output_fa(BackT,V,W,I,J,AccV,AccW) when J =:= 0 ->
	{AccV,AccW};
output_fa(BackT,V,W,I,J,AccV,AccW) ->
	
	case util:mx_at(I,J,BackT) of
		d ->
			output_fa(BackT,V,W,I-1,J,[lists:nth(I,V)|AccV],[$-|AccW]);
		i ->
			output_fa(BackT,V,W,I,J-1,[$-|AccV],[lists:nth(J,W)|AccW]);
		m ->
			output_fa(BackT,V,W,I-1,J-1,[lists:nth(I,V)|AccV],[lists:nth(J,W)|AccW])
	end.


%% -------------------------------------------------------------------------- %%
%% Overlappign align
%% -------------------------------------------------------------------------- %%

print_oa(V,W,IndelPenalty) ->
		{S,B} = overlap_align(V,W,IndelPenalty),
		
		RightmostCol = length(W),
		BottomRow = length(V),
		M1 = util:mx_max_vector(util:mx_col_at(RightmostCol,S)), 
		M2 = util:mx_max_vector(util:mx_row_at(BottomRow,S)), 
		{Max1,I} = M1,
		{Max2,J} = M2,
		case Max1 of
			Max1 when Max1 > Max2 ->
				Max = Max1,
				StartRow = I,
				StartCol = RightmostCol;
			_ ->
				Max = Max2,
				StartRow = BottomRow,
				StartCol = J
		end,
		{VA,WA} = output_oa(B,V,W,StartRow,StartCol),
		
		io:format("~p~n~s~n~s~n",[Max,VA,WA]).

overlap_align(V,W,IndelPenalty) when is_list(V), is_list(W) ->
	S0 = util:new_mx(length(V)+1,length(W)+1), % init |V+1| x |W+1| matrix,
	B0 = util:new_mx(length(V)+1,length(W)+1,x), % Backtrack matrix
	S1 = util:mx_col(0,lists:duplicate(length(V)+1,0),S0), % init col 0
	S2 = util:mx_row(0,lists:duplicate(length(W)+1,0),S1), % init row 0
	overlap_align(lists:seq(1,length(V)),lists:seq(1,length(W)),V,W,IndelPenalty,S2,B0).

overlap_align([I|IT],[],V,W,IndelP,S,B) ->
	overlap_align(IT,lists:seq(1,length(W)),V,W,IndelP,S,B);

overlap_align([],_,_,_,_,S,B) ->
	{S,B};

overlap_align([I|IT],[J|JT],V,W,IndelP,S,B) ->
	{SVal,BVal} = o_value(I,J,V,W,IndelP,S),
	NewS = util:mx(I,J,SVal,S),
	NewB = util:mx(I,J,BVal,B),
	overlap_align([I|IT],JT,V,W,IndelP,NewS,NewB).


o_value(I,J,V,W,IndelP,S) ->
	VI = lists:nth(I,V),
	VJ = lists:nth(J,W),
	Score = case VI of
		VJ ->
			1;
		_ ->
			IndelP
		end,
	A = util:mx_at(I-1,J,S)+IndelP,
	B = util:mx_at(I,J-1,S)+IndelP,
	C = util:mx_at(I-1,J-1,S)+Score,
	SVal = util:v_max(A,B,C),
	case SVal of
		SVal when SVal =:= A ->
			{SVal,d};
		SVal when SVal =:= B ->
			{SVal,i};
		SVal when SVal =:= C ->
			{SVal,m}
	end.

output_oa(BackT,V,W,I,J) ->
	output_oa(BackT,V,W,I,J,[],[]).

output_oa(BackT,V,W,I,J,AccV,AccW) when J =:= 0 ->
	{AccV,AccW};
output_oa(BackT,V,W,I,J,AccV,AccW) ->
	
	case util:mx_at(I,J,BackT) of
		d ->
			output_oa(BackT,V,W,I-1,J,[lists:nth(I,V)|AccV],[$-|AccW]);
		i ->
			output_oa(BackT,V,W,I,J-1,[$-|AccV],[lists:nth(J,W)|AccW]);
		m ->
			output_oa(BackT,V,W,I-1,J-1,[lists:nth(I,V)|AccV],[lists:nth(J,W)|AccW])
	end.


%% -------------------------------------------------------------------------- %%
%% Middle edge in space efficient alignment
%% -------------------------------------------------------------------------- %%

middle_edge(V,W,IndelP) ->
	
	Col0 = lists:seq(IndelP,length(V)*IndelP,IndelP),
	Row0 = lists:seq(IndelP,length(W)*IndelP,IndelP),
	
	MiddleColNum = length(W) div 2,
	NextColNum = MiddleColNum+1,

	FwdMiddleC = get_column(MiddleColNum,V,W,Col0,Row0,0,IndelP),
	RevMiddleC = lists:reverse(get_column(NextColNum,lists:reverse(V),lists:reverse(W),Col0,Row0,0,IndelP)),

	{Type,Max,I} = middle_node(FwdMiddleC,RevMiddleC),

	case Type of
		diag ->
			io:format("(~p,~p) (~p,~p)~n",[I,MiddleColNum,I+1,NextColNum]);
		horiz ->
			io:format("(~p,~p) (~p,~p)~n",[I,MiddleColNum,I,NextColNum])
	end.
	
% ColumnN
get_column(ColN,V,W,Col0,Row0,TopLeftVal,IndelP) ->
	get_col_score(1,ColN,V,W,[TopLeftVal|Col0],Row0,IndelP).

get_col_score(CurrentCol,ColN,V,W,InitCol,Row0,IndelP) ->
	ColScore = col_score(InitCol,lists:nth(CurrentCol,Row0),V,lists:nth(CurrentCol,W),IndelP),
	case CurrentCol of
		CurrentCol when CurrentCol =:= ColN ->
			ColScore;
		_ ->
			get_col_score(CurrentCol+1,ColN,V,W,ColScore,Row0,IndelP)
	end.

middle_node([AH|AColT],[BH|BColT]) ->
	middle_node(AColT,BH,BColT,1,{diag,AH+BH,0}).
middle_node([],_,[],_,Max) ->
	Max;
middle_node([AH|AColT],LastBH,[BH|BColT],Pos,Max) ->
	{_,MaxV,_} = Max,
	TmpMax = lists:last(lists:keysort(2,[{diag,AH+BH,Pos},{horiz,AH+LastBH,Pos}])),
	NewMax = case TmpMax of
					{_,V,_} when V > MaxV ->
						TmpMax;
					_ -> 
						Max
			end,
	middle_node(AColT,BH,BColT,Pos+1,NewMax).

%      Wchar
%      A | B
%     --------
%     TA | TB
% V1  A1 | B1
% V2  A2 | B2
% V3  A3 | ?
%
% Calculating B3 - DiagonalValue = A2, LeftValue = A3, UpValue = B2
% InitialA = [TA,A1,A2,A3,A4...An]
% TopSeed = TB
%
col_score(InitialA,TopSeed,V,Wchar,IndelP) ->
	col_score(hd(InitialA),hd(tl(InitialA)),TopSeed,V,Wchar,IndelP,tl(tl(InitialA)),[TopSeed]).

col_score(DiagValue,LeftValue,UpValue,[VH|VT],Wchar,IndelP,[AH|AT],B) ->
	CurrentB = current_val(DiagValue,LeftValue,UpValue,VH,Wchar,IndelP),
	% Next: DiagValue = OldLeft, Left = Head of A, Up = OldCurrentB
	case AT of
		[] ->
			NextB = current_val(LeftValue,AH,CurrentB,hd(VT),Wchar,IndelP),
			lists:reverse([NextB,CurrentB|B]);
		_ ->
			col_score(LeftValue,AH,CurrentB,VT,Wchar,IndelP,AT,[CurrentB|B])
	end.
	

current_val(DiagValue,LeftValue,UpValue,Vchar,Wchar,IndelP) ->
	X = LeftValue + IndelP,
	Y = UpValue + IndelP,
	Z = DiagValue + util:mx_at(b62(Vchar),b62(Wchar),blosum62()), % todo: use a list since Wchar never changes
	utils:v_max(X,Y,Z).


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
%% Alignment scoring matrices
%% -------------------------------------------------------------------------- %%

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



pam250() ->[[ 2, -2,  0,  0, -3,  1, -1, -1, -1, -2, -1,  0,  1,  0, -2,  1,  1,  0, -6, -3],
            [-2, 12, -5, -5, -4, -3, -3, -2, -5, -6, -5, -4, -3, -5, -4,  0, -2, -2, -8,  0],
            [ 0, -5,  4,  3, -6,  1,  1, -2,  0, -4, -3,  2, -1,  2, -1,  0,  0, -2, -7, -4],
            [ 0, -5,  3,  4, -5,  0,  1, -2,  0, -3, -2,  1, -1,  2, -1,  0,  0, -2, -7, -4],
            [-3, -4, -6, -5,  9, -5, -2,  1, -5,  2,  0, -3, -5, -5, -4, -3, -3, -1,  0,  7],
            [ 1, -3,  1,  0, -5,  5, -2, -3, -2, -4, -3,  0,  0, -1, -3,  1,  0, -1, -7, -5],
            [-1, -3,  1,  1, -2, -2,  6, -2,  0, -2, -2,  2,  0,  3,  2, -1, -1, -2, -3,  0],
            [-1, -2, -2, -2,  1, -3, -2,  5, -2,  2,  2, -2, -2, -2, -2, -1,  0,  4, -5, -1],
            [-1, -5,  0,  0, -5, -2,  0, -2,  5, -3,  0,  1, -1,  1,  3,  0,  0, -2, -3, -4],
            [-2, -6, -4, -3,  2, -4, -2,  2, -3,  6,  4, -3, -3, -2, -3, -3, -2,  2, -2, -1],
            [-1, -5, -3, -2,  0, -3, -2,  2,  0,  4,  6, -2, -2, -1,  0, -2, -1,  2, -4, -2],
            [ 0, -4,  2,  1, -3,  0,  2, -2,  1, -3, -2,  2,  0,  1,  0,  1,  0, -2, -4, -2],
            [ 1, -3, -1, -1, -5,  0,  0, -2, -1, -3, -2,  0,  6,  0,  0,  1,  0, -1, -6, -5],
            [ 0, -5,  2,  2, -5, -1,  3, -2,  1, -2, -1,  1,  0,  4,  1, -1, -1, -2, -5, -4],
            [-2, -4, -1, -1, -4, -3,  2, -2,  3, -3,  0,  0,  0,  1,  6,  0, -1, -2,  2, -4],
            [ 1,  0,  0,  0, -3,  1, -1, -1,  0, -3, -2,  1,  1, -1,  0,  2,  1, -1, -2, -3],
            [ 1, -2,  0,  0, -3,  0, -1,  0,  0, -2, -1,  0,  0, -1, -1,  1,  3,  0, -5, -3],
            [ 0, -2, -2, -2, -1, -1, -2,  4, -2,  2,  2, -2, -1, -2, -2, -1,  0,  4, -6, -2],
            [-6, -8, -7, -7,  0, -7, -3, -5, -3, -2, -4, -4, -6, -5,  2, -2, -5, -6, 17,  0],
            [-3,  0, -4, -4,  7, -5,  0, -1, -4, -1, -2, -2, -5, -4, -4, -3, -3, -2,  0, 10]].

% short map
p250(X) -> pam250_rc_map(X).

% map protein alphabet to pam250 cost matrix row/col
pam250_rc_map($A) -> 0;
pam250_rc_map($C) -> 1;
pam250_rc_map($D) -> 2;
pam250_rc_map($E) -> 3;
pam250_rc_map($F) -> 4;
pam250_rc_map($G) -> 5;
pam250_rc_map($H) -> 6;
pam250_rc_map($I) -> 7;
pam250_rc_map($K) -> 8;
pam250_rc_map($L) -> 9;
pam250_rc_map($M) -> 10;
pam250_rc_map($N) -> 11;
pam250_rc_map($P) -> 12;
pam250_rc_map($Q) -> 13;
pam250_rc_map($R) -> 14;
pam250_rc_map($S) -> 15;
pam250_rc_map($T) -> 16;
pam250_rc_map($V) -> 17;
pam250_rc_map($W) -> 18;
pam250_rc_map($Y) -> 19.



















