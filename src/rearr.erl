%
% 2014-01-07 karlll <karl@ninjacontrol.com>
%

%
% Rearranging sequences
% 

-module(rearr).
-compile([export_all]).


%% -------------------------------------------------------------------------- %%
%% Greeding sorting by reversals
%% -------------------------------------------------------------------------- %%

greedy_sort(P) ->
	greedy_sort(P,0,1,[],[]).

greedy_sort([],Dist,_,_Acc,Perms) ->
	{Dist,Perms};
greedy_sort([PH|PT],Dist,Pos,Acc,Perms)  when PH =:= Pos ->
	greedy_sort(PT,Dist,Pos+1,[PH|Acc],Perms);
greedy_sort([PH|PT],Dist,Pos,Acc,Perms) when -PH =:= Pos ->
	NewAcc = [-PH|Acc],
	Perm = lists:reverse(NewAcc) ++ PT,
	greedy_sort(PT,Dist+1,Pos+1, NewAcc, [Perm|Perms]);
greedy_sort([PH|PT],Dist,Pos,Acc,Perms) ->
	Rev = rev([PH|PT],1,util:idx_abs(Pos,PT)+1),
	Perm = lists:reverse(Acc) ++ Rev,
	greedy_sort(Rev,Dist+1,Pos, Acc, [Perm|Perms]).

%% -------------------------------------------------------------------------- %%
%% Number of breakpoints
%% -------------------------------------------------------------------------- %%

breakpoints(P) ->
	P2 = P ++ [length(P)+1],
	case hd(P) of
		0 ->
			breakpoints(P2,0);
		_ ->
			breakpoints([0|P2],0)
	end.

breakpoints([P|[]],Count) ->
	Count;
breakpoints([Pn,Pn1|PT],Count) when Pn1 =:= Pn+1 ->
	breakpoints([Pn1|PT],Count);
breakpoints([P|PT],Count) ->
	breakpoints(PT,Count+1).


%% -------------------------------------------------------------------------- %%
%% 2-break distance
%% -------------------------------------------------------------------------- %%

twobreak_dist(PermList) ->
	AllNodes = sets:to_list(sets:from_list( % unique
				lists:flatten(
					lists:map(fun(P) -> perm_to_nodelist(P) end,PermList)))
			   ),
	AllEdges = sets:to_list(sets:from_list(lists:flatten(
				lists:map(fun(P) -> perm_to_edgelist(P) end, PermList)
			   ))),
	Blocks = length(AllNodes) div 2,
	Cycles = count_cycles(AllEdges),
	{Blocks, Cycles, Blocks - Cycles}.



% A permutation to a circular edge list; (1,2,-3) -> 
% [{'2h','3h'},{'1h','2t'},{'1t','3t'}]
perm_to_edgelist(Perms) ->
	nodelist_to_edgelist(perm_to_nodelist(Perms)).

perm_to_nodelist(Perms) ->
	
	NL = lists:map(fun(E) ->
				H = list_to_atom(integer_to_list(abs(E))++"h"),
				T = list_to_atom(integer_to_list(abs(E))++"t"),
				case E of
					E when E >= 0 ->
						[H,T];
					E when E < 0 ->
						[T,H]
				end end, Perms
			),
	lists:flatten(NL).

nodelist_to_edgelist(NL) ->
	{L1,L2} = lists:split(1,NL),
	{L3,L4} = lists:split(length(L2)-1,L2),
	nl_to_el(L3,[{hd(L4),hd(L1)},{hd(L1),hd(L4)}]).

nl_to_el([],Acc) ->
	Acc;
nl_to_el([A,B|T],Acc) ->
	nl_to_el(T,[{A,B},{B,A}|Acc]).


count_cycles(Edges) ->
	EFlag = lists:zip(Edges,lists:duplicate(length(Edges),false)),
	count_cycles(EFlag,0).

count_cycles([],Cycles) ->
	Cycles;
count_cycles([{Start,false}|EL],Cycles) ->
	EL2 = follow({Start,false},[{Start,false}|EL]),
	EL3 = clear_visited(EL2),
	count_cycles(EL3,Cycles+1).

clear_visited(Edges) ->
	lists:filter(fun(E) -> case E of
						{_,false} ->
							true;
						_ ->
							false
						end
				end,
				Edges).



% follows an edge until Current edge = start edge. Mark visited along the way

follow({{A,B},Flag},[{{A,B},Flag},{{B,A},Flag}]) -> % only one edge
	[];

follow(Start,EL) ->
	{{A,B},false} = Start,
	StartRev = {{B,A},false},
	EL2 = mark_edge(Start,EL),
	EL3 = mark_edge(StartRev,EL2),
	Next = find_next(Start,EL3),
	follow(Next,Start,EL).

follow({{A,B},_},{{A,B},_},EL) ->
	Start = {{A,B},false},
	StartRev = {{B,A},false},
	EL2 = mark_edge(Start,EL),
	EL3 = mark_edge(StartRev,EL2),
	EL3;

follow(Current,Start,EL) ->
	{{A,B},false} = Current,
	CurrentRev = {{B,A},false},
	EL2 = mark_edge(Current,EL),
	EL3 = mark_edge(CurrentRev,EL2),
	Next = find_next(Current,EL3),
	follow(Next,Start,EL3).


find_next(_,[]) ->
	not_found;
find_next({{A,B},_},[{{B,C},false}|EL]) ->
	{{B,C},false};
find_next({{A,B},F},[{{_C,_D},_}|EL]) ->
	find_next({{A,B},F},EL).

mark_edge(E,EL) ->
	mark_edge(E,EL,[]).
mark_edge(_,[],Acc) ->
	Acc;
mark_edge(E,[E|EL],Acc) ->
	{{A,B},_} = E,
	mark_edge(E,EL,[{{A,B},true}|Acc]);
mark_edge(E,[{H,Flag}|EL],Acc) ->
	mark_edge(E,EL,[{H,Flag}|Acc]).

%% -------------------------------------------------------------------------- %%

rev(El) when is_integer(El) ->
	-El;
rev(List) when is_list(List) ->
	lists:map(fun(El) -> -El end, lists:reverse(List)).

rev(List,Start,End) ->
	rev(List,Start,End,1,[],[],[]).

rev([H|T],Start,End,Pos,A,B,C) when Pos < Start ->
	rev(T,Start,End,Pos+1,[H|A],B,C);
rev([H|T],Start,End,Pos,A,B,C) when Pos >= Start, Pos =< End ->
	rev(T,Start,End,Pos+1,A,[-H|B],C);
rev([H|T],Start,End,Pos,A,B,C) when Pos > End ->
	rev(T,Start,End,Pos+1,A,B,[H|C]);
rev([],_,_,_,A,B,C) ->
	lists:reverse(A) ++ B ++ lists:reverse(C).

% Identity permutation of length L
id_perm(Len) when is_integer(Len) -> lists:seq(1,Len).
is_id_perm(P) when is_list(P) -> 
	IdP = id_perm(length(P)),
	case P of
		IdP ->
			true;
		_ ->
			false
	end.
